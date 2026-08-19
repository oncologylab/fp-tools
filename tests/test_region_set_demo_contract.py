import gzip
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/demos/data"
BROWSER = ROOT / "docs/demos/reports/region_set_HepG2_HNF4A_FOXA2"
GROUP_SLUGS = (
    "HNF4A_FOXA2",
    "HNF4A_only",
    "FOXA2_only",
    "No_HNF4A_FOXA2",
)
DEFAULT_MOTIFS = (
    "MA1494.2",
    "MA0484.3",
    "MA0047.4",
    "MA0148.5",
    "MA0046.3",
    "MA0153.2",
    "MA0102.5",
    "MA0466.4",
)


class RegionSetDemoContractTest(unittest.TestCase):
    def test_matched_regions_are_balanced_and_nonoverlapping(self):
        coordinates = set()
        stratum_counts = []
        for slug in GROUP_SLUGS:
            path = DATA / f"region_set_HepG2_{slug}.bed.gz"
            self.assertTrue(path.is_file())
            counts = {}
            rows = 0
            with gzip.open(path, "rt") as handle:
                for line in handle:
                    chrom, start, end, stratum, _signal = line.rstrip("\n").split("\t")
                    coordinate = (chrom, int(start), int(end))
                    self.assertNotIn(coordinate, coordinates)
                    coordinates.add(coordinate)
                    counts[stratum] = counts.get(stratum, 0) + 1
                    rows += 1
            self.assertEqual(rows, 4786)
            self.assertEqual(len(counts), 50)
            stratum_counts.append(counts)
        self.assertTrue(all(counts == stratum_counts[0] for counts in stratum_counts[1:]))

    def test_browser_contains_complete_results_and_curated_defaults(self):
        metadata = json.loads((BROWSER / "data/metadata.json").read_text())
        self.assertEqual(
            metadata["default_comparison"],
            {"condition1": "HNF4A + FOXA2", "condition2": "No HNF4A/FOXA2"},
        )
        self.assertEqual(metadata["default_aggregate_plots"], 8)
        self.assertEqual(len(metadata["comparisons"]), 6)
        self.assertEqual({record["motifs"] for record in metadata["comparisons"]}, {1019})
        self.assertEqual({record["aggregate_motifs"] for record in metadata["comparisons"]}, {8})

        default_record = next(
            record
            for record in metadata["comparisons"]
            if {
                record["condition1"],
                record["condition2"],
            }
            == {"HNF4A + FOXA2", "No HNF4A/FOXA2"}
        )
        with gzip.open(BROWSER / default_record["file"], "rt") as handle:
            payload = json.load(handle)
        motifs = (payload.get("aggregate") or {}).get("motifs") or []
        self.assertEqual([motif["motif_id"] for motif in motifs], list(DEFAULT_MOTIFS))
        self.assertEqual(metadata["default_aggregate_motifs"], [motif["prefix"] for motif in motifs])
        for motif in motifs:
            counts = [int(condition["n_sites"]) for condition in motif["conditions"]]
            self.assertGreaterEqual(min(counts), 150)
            self.assertGreaterEqual(sum(counts), 500)

    def test_source_manifest_records_only_compact_public_inputs(self):
        manifest = (DATA / "region_set_HepG2_source_manifest.tsv").read_text()
        for accession in (
            "ENCFF624SON",
            "ENCFF926KFU",
            "ENCFF990VCP",
            "ENCFF609BSU",
            "ENCFF704BPD",
            "ENCFF656PGC",
        ):
            self.assertIn(accession, manifest)
        self.assertNotIn("FASTQ", manifest)
        self.assertNotIn("local_path", manifest)


if __name__ == "__main__":
    unittest.main()
