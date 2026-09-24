#!/usr/bin/env python3
"""Refresh report presentation without recomputing or rewriting scientific payloads."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from fp_tools.tools.review_multi_comparisons import _decode_payload_b64
from fp_tools.tools.static_comparison_browser import write_embedded_static_browser


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encoded_report_payload(path):
    text = path.read_text(encoding='utf-8')
    match = re.search(r'const\s+reportPayloadB64\s*=\s*"([^"]+)"', text)
    if not match:
        raise ValueError(f'No embedded report payload in {path.name}')
    return match.group(1)


def read_report_payload(path):
    return _decode_payload_b64(encoded_report_payload(path))


def rebuild(source: Path, destination: Path, review_view: str = 'side') -> dict:
    encoded = encoded_report_payload(source)
    payload = _decode_payload_b64(encoded)
    if payload.get('schema') == 'fp-tools.review-multi-comparisons.v1':
        write_embedded_static_browser(payload, destination, default_view=review_view, encoded_payload=encoded)
        comparisons = len(payload['comparisons'])
    elif 'points' in payload and 'conditions' in payload:
        review = dict(schema='fp-tools.review-multi-comparisons.v1', title=payload.get('title',''),
                      comparisons=[dict(label=payload.get('title','Comparison'),payload=payload)])
        write_embedded_static_browser(review, destination, source_payload=payload, encoded_payload=encoded)
        comparisons = 1
    else:
        raise ValueError(f'Not a differential/review report: {source.name}')
    # Keep site-specific branding when refreshing a documentation example.
    with source.open(encoding='utf-8') as handle:
        head = handle.read(8192)
    icons = re.findall(r'''<link\b[^>]*\brel=["']icon["'][^>]*>''', head, re.IGNORECASE)
    if icons:
        document = destination.read_text(encoding='utf-8')
        destination.write_text(document.replace('<head>', '<head>\n' + '\n'.join(icons), 1), encoding='utf-8')
    rebuilt = read_report_payload(destination)
    if rebuilt != payload:
        raise ValueError(f'Scientific payload changed: {source.name}')
    return dict(source_sha256=sha256(source), output_sha256=sha256(destination),
                payload_sha256=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                payload_equal=True,comparisons=comparisons)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--review-view',choices=['single','side'],default='side')
    args=parser.parse_args()
    if args.input_dir.resolve()==args.output_dir.resolve():
        parser.error('Use a separate output directory to preserve the original HTML files.')
    files=sorted(args.input_dir.rglob('*.html'))
    if not files:parser.error('No HTML reports found')
    receipts={}
    for source in files:
        relative=source.relative_to(args.input_dir)
        receipts[relative.as_posix()]=rebuild(source,args.output_dir/relative,args.review_view)
        print(f'Payload unchanged: {relative}',flush=True)
    (args.output_dir/'refresh_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print(f'Verified {len(receipts)} refreshed HTML payloads',flush=True)


if __name__=='__main__':main()
