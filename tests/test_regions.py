import unittest

from fp_tools.utils.regions import OneRegion, RegionList


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


if __name__ == "__main__":
    unittest.main()
