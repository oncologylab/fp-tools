#!/usr/bin/env bash
set -euo pipefail

checksum() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | cut -d ' ' -f 1
  else
    md5 -q "$1"
  fi
}

download() {
  local expected="$1" destination="$2" accession filename status
  if [[ -f "$destination" ]] && [[ "$(checksum "$destination")" == "$expected" ]]; then
    echo "Verified cached $destination"
    return
  fi
  filename="${destination##*/}"
  accession="${filename%%.*}"
  local url="https://www.encodeproject.org/files/$accession/@@download/$filename"
  # Keep unfinished transfers separate from validated inputs.
  status=0
  curl -fL --retry 3 -C - -o "$destination.part" "$url" || status=$?
  if [[ "$status" == 33 || "$status" == 36 ]]; then
    echo "Server cannot resume $filename; restarting this file."
    curl -fL --retry 3 -o "$destination.part" "$url" || return $?
  elif [[ "$status" != 0 ]]; then
    return "$status"
  fi
  if [[ "$(checksum "$destination.part")" != "$expected" ]]; then
    echo "Checksum mismatch: $destination.part. Remove this partial file before retrying." >&2
    return 1
  fi
  mv "$destination.part" "$destination"
}

main() {
  local tool
  for tool in curl gzip samtools cut mkdir mv; do
    command -v "$tool" >/dev/null 2>&1 || { echo "Required command missing: $tool" >&2; return 1; }
  done
  if ! command -v md5sum >/dev/null 2>&1 && ! command -v md5 >/dev/null 2>&1; then
    echo "Required checksum command missing: md5sum (Linux) or md5 (macOS)" >&2
    return 1
  fi
  mkdir -p encode_data/bams encode_data/peaks
  while read -r expected destination; do
    download "$expected" "$destination" || return $?
  done <<'CHECKSUMS'
4a18e40ee643905cf62ae1fbbd33fb87  encode_data/bams/ENCFF624SON.bam
debb5a616ce1ca26a0d956ce54eae02b  encode_data/bams/ENCFF926KFU.bam
7cc7c6739736b169f56e966784cefcde  encode_data/bams/ENCFF990VCP.bam
135e5a65e5eafac97a8a2242b4c13963  encode_data/bams/ENCFF077FBI.bam
ed3d553646995be008241f6e75f57042  encode_data/bams/ENCFF128WZG.bam
4bd9c1e48617afea2a2cf4834eed961b  encode_data/bams/ENCFF534DCE.bam
4020785aaa90efbebc5a3fc92b5a9cff  encode_data/peaks/ENCFF536RJV.bed.gz
bee81b0efe50b144c24d2060a87aa843  encode_data/peaks/ENCFF855PCP.bed.gz
CHECKSUMS

gzip -dc encode_data/peaks/ENCFF536RJV.bed.gz > encode_data/peaks/ENCFF536RJV.bed
gzip -dc encode_data/peaks/ENCFF855PCP.bed.gz > encode_data/peaks/ENCFF855PCP.bed

  for bam in encode_data/bams/*.bam; do
    samtools index "$bam"
  done
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main
fi
