from fp_tools.tools.find_signature_fp import DEFAULT_KNN_NEIGHBORS
from fp_tools.tools.pseudobulk_footprints import DEFAULT_KNN_NEIGHBORS as PSEUDOBULK_DEFAULT_KNN_NEIGHBORS


def test_reviewer_validated_knn_default_is_shared() -> None:
    assert DEFAULT_KNN_NEIGHBORS == 150
    assert PSEUDOBULK_DEFAULT_KNN_NEIGHBORS == DEFAULT_KNN_NEIGHBORS
