import unittest
import random

from fp_tools.utils.regions import OneRegion, RegionCluster, RegionList


class RegionListMergeTest(unittest.TestCase):
    @staticmethod
    def _coordinates(regions):
        return [(region.chrom, region.start, region.end) for region in regions]

    def test_contained_interval_preserves_outer_endpoint(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 10, "outer"]),
                OneRegion(["chr1", 2, 8, "inner"]),
            ]
        )

        regions.merge()

        self.assertEqual(self._coordinates(regions), [("chr1", 0, 10)])

    def test_identical_intervals_merge(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 10, "first"]),
                OneRegion(["chr1", 0, 10, "second"]),
            ]
        )

        regions.merge()

        self.assertEqual(self._coordinates(regions), [("chr1", 0, 10)])

    def test_contained_middle_interval_does_not_break_chained_overlap(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 10]),
                OneRegion(["chr1", 2, 8]),
                OneRegion(["chr1", 9, 12]),
            ]
        )

        regions.merge()

        self.assertEqual(self._coordinates(regions), [("chr1", 0, 12)])

    def test_non_overlapping_and_bookended_intervals_remain_separate(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 5]),
                OneRegion(["chr1", 5, 10]),
                OneRegion(["chr1", 12, 15]),
                OneRegion(["chr2", 0, 5]),
            ]
        )

        regions.merge()

        self.assertEqual(
            self._coordinates(regions),
            [
                ("chr1", 0, 5),
                ("chr1", 5, 10),
                ("chr1", 12, 15),
                ("chr2", 0, 5),
            ],
        )

    def test_name_aware_merge_keeps_different_names_separate(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 10, "A"]),
                OneRegion(["chr1", 2, 8, "B"]),
            ]
        )

        regions.merge(name=True)

        self.assertEqual(self._coordinates(regions), [("chr1", 0, 10), ("chr1", 2, 8)])

    def test_name_aware_merge_handles_interleaved_names(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 10, "A"]),
                OneRegion(["chr1", 2, 8, "B"]),
                OneRegion(["chr1", 4, 6, "A"]),
            ]
        )

        regions.merge(name=True)

        self.assertEqual(
            [(region.chrom, region.start, region.end, region.name) for region in regions],
            [("chr1", 0, 10, "A"), ("chr1", 2, 8, "B")],
        )

    def test_count_overlaps_uses_name_aware_geometric_unions(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 10, "A"]),
                OneRegion(["chr1", 2, 8, "B"]),
                OneRegion(["chr1", 4, 6, "A"]),
            ]
        )

        overlaps = regions.count_overlaps()

        self.assertEqual(overlaps["A"], 10)
        self.assertEqual(overlaps["B"], 6)
        self.assertEqual(overlaps[("A", "B")], 6)
        self.assertEqual(overlaps[("B", "A")], 6)

    def test_count_overlaps_includes_one_base_intersection(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 5, "A"]),
                OneRegion(["chr1", 4, 10, "B"]),
            ]
        )

        overlaps = regions.count_overlaps()

        self.assertEqual(overlaps[("A", "B")], 1)
        self.assertEqual(overlaps[("B", "A")], 1)
        cluster = RegionCluster(overlaps)
        cluster.overlap_to_distance()
        self.assertAlmostEqual(cluster.distance_mat[0, 1], 0.8)
        self.assertAlmostEqual(cluster.distance_mat[1, 0], 0.8)

    def test_count_overlaps_obeys_half_open_bed_boundaries(self):
        regions = RegionList(
            [
                OneRegion(["chr1", 0, 5, "A"]),
                OneRegion(["chr1", 5, 10, "B"]),
                OneRegion(["chr1", 8, 12, "C"]),
                OneRegion(["chr1", 9, 10, "D"]),
                OneRegion(["chr2", 4, 10, "E"]),
            ]
        )

        overlaps = regions.count_overlaps()

        self.assertNotIn(("A", "B"), overlaps)
        self.assertEqual(overlaps[("B", "C")], 2)
        self.assertEqual(overlaps[("B", "D")], 1)
        self.assertEqual(overlaps[("D", "B")], 1)
        self.assertNotIn(("A", "E"), overlaps)

    def test_count_overlaps_matches_randomized_half_open_intersections(self):
        rng = random.Random(20260820)
        for _ in range(2000):
            start_a = rng.randrange(0, 200)
            start_b = rng.randrange(0, 200)
            end_a = start_a + rng.randrange(1, 30)
            end_b = start_b + rng.randrange(1, 30)
            regions = RegionList(
                [
                    OneRegion(["chr1", start_a, end_a, "A"]),
                    OneRegion(["chr1", start_b, end_b, "B"]),
                ]
            )

            overlaps = regions.count_overlaps()
            expected = max(0, min(end_a, end_b) - max(start_a, start_b))
            self.assertEqual(overlaps.get(("A", "B"), 0), expected)
            self.assertEqual(overlaps.get(("B", "A"), 0), expected)


if __name__ == "__main__":
    unittest.main()
