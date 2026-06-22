"""Bundled motif database registry and path resolution."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable


DEFAULT_MOTIF_DB = "jaspar2026_vertebrates"


@dataclass(frozen=True)
class MotifDatabase:
    key: str
    label: str
    filename: str
    source_url: str
    license: str
    citation: str


_JASPAR_BASE = "https://jaspar.elixir.no/download/data/2026/CORE"
_HOCOMOCO_BASE = "https://hocomoco14.autosome.org/final_bundle/hocomoco14"


MOTIF_DATABASES: dict[str, MotifDatabase] = {
    "jaspar2026_vertebrates": MotifDatabase(
        "jaspar2026_vertebrates",
        "JASPAR 2026 CORE vertebrates non-redundant",
        "JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "jaspar2026_urochordates": MotifDatabase(
        "jaspar2026_urochordates",
        "JASPAR 2026 CORE urochordates non-redundant",
        "JASPAR2026_CORE_urochordates_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_urochordates_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "jaspar2026_plants": MotifDatabase(
        "jaspar2026_plants",
        "JASPAR 2026 CORE plants non-redundant",
        "JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "jaspar2026_nematodes": MotifDatabase(
        "jaspar2026_nematodes",
        "JASPAR 2026 CORE nematodes non-redundant",
        "JASPAR2026_CORE_nematodes_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_nematodes_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "jaspar2026_insects": MotifDatabase(
        "jaspar2026_insects",
        "JASPAR 2026 CORE insects non-redundant",
        "JASPAR2026_CORE_insects_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_insects_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "jaspar2026_fungi": MotifDatabase(
        "jaspar2026_fungi",
        "JASPAR 2026 CORE fungi non-redundant",
        "JASPAR2026_CORE_fungi_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_fungi_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "jaspar2026_core": MotifDatabase(
        "jaspar2026_core",
        "JASPAR 2026 CORE all non-redundant",
        "JASPAR2026_CORE_non-redundant_pfms_jaspar.txt",
        f"{_JASPAR_BASE}/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt",
        "CC BY 4.0",
        "JASPAR 2026",
    ),
    "hocomoco14_core": MotifDatabase(
        "hocomoco14_core",
        "HOCOMOCO v14 CORE",
        "H14CORE_jaspar_format.txt",
        f"{_HOCOMOCO_BASE}/H14CORE/formatted_motifs/H14CORE_jaspar_format.txt",
        "WTFPL; may be treated as CC-BY",
        "HOCOMOCO v14",
    ),
    "hocomoco14_invivo": MotifDatabase(
        "hocomoco14_invivo",
        "HOCOMOCO v14 INVIVO",
        "H14INVIVO_jaspar_format.txt",
        f"{_HOCOMOCO_BASE}/H14INVIVO/formatted_motifs/H14INVIVO_jaspar_format.txt",
        "WTFPL; may be treated as CC-BY",
        "HOCOMOCO v14",
    ),
    "hocomoco14_invitro": MotifDatabase(
        "hocomoco14_invitro",
        "HOCOMOCO v14 INVITRO",
        "H14INVITRO_jaspar_format.txt",
        f"{_HOCOMOCO_BASE}/H14INVITRO/formatted_motifs/H14INVITRO_jaspar_format.txt",
        "WTFPL; may be treated as CC-BY",
        "HOCOMOCO v14",
    ),
    "hocomoco14_rsnp": MotifDatabase(
        "hocomoco14_rsnp",
        "HOCOMOCO v14 RSNP",
        "H14RSNP_jaspar_format.txt",
        f"{_HOCOMOCO_BASE}/H14RSNP/formatted_motifs/H14RSNP_jaspar_format.txt",
        "WTFPL; may be treated as CC-BY",
        "HOCOMOCO v14",
    ),
}


MOTIF_DB_ALIASES = {
    "jaspar2026": "jaspar2026_vertebrates",
    "jaspar2026_all": "jaspar2026_core",
    "jaspar2026_nonredundant": "jaspar2026_core",
    "hocomoco14": "hocomoco14_core",
    "h14core": "hocomoco14_core",
    "h14invivo": "hocomoco14_invivo",
    "h14invitro": "hocomoco14_invitro",
    "h14rsnp": "hocomoco14_rsnp",
}


def normalize_motif_db_key(name: str) -> str:
    key = str(name).strip().lower().replace("-", "_")
    key = MOTIF_DB_ALIASES.get(key, key)
    if key not in MOTIF_DATABASES:
        choices = ", ".join(sorted(MOTIF_DATABASES))
        raise ValueError(f"Unknown built-in motif database '{name}'. Available choices: {choices}")
    return key


def motif_db_path(name: str) -> Path:
    db = MOTIF_DATABASES[normalize_motif_db_key(name)]
    return Path(resources.files("fp_tools.resources.motifs").joinpath(db.filename))


def resolve_motif_inputs(
    motifs: Iterable[str | Path] | str | Path | None,
    motif_db: str | None = None,
    *,
    use_default: bool = True,
) -> list[str]:
    if isinstance(motifs, (str, Path)):
        motif_values = [motifs]
    else:
        motif_values = list(motifs or [])
    paths = [str(path) for path in motif_values if str(path).strip()]
    if motif_db:
        paths.insert(0, str(motif_db_path(motif_db)))
    elif use_default and not paths:
        paths.append(str(motif_db_path(DEFAULT_MOTIF_DB)))
    return paths


def motif_db_table() -> str:
    rows = ["Available built-in motif databases:", ""]
    for key in sorted(MOTIF_DATABASES):
        db = MOTIF_DATABASES[key]
        rows.append(f"{key}\t{db.label}\t{db.license}")
    rows.append("")
    rows.append(f"Default: {DEFAULT_MOTIF_DB}")
    rows.append("Aliases: " + ", ".join(f"{alias}={target}" for alias, target in sorted(MOTIF_DB_ALIASES.items())))
    return "\n".join(rows)
