# TODO: write test code for utils functions
from synbio.utils import *


def test_utils():
    pass


def test_find_subseq():
    seq_to_search = "xxXatATaTxXatatxx"
    subseq_to_find = "atat"  # search is NOT case sensitive

    subseq_ix = find_subseq(seq_to_search, subseq_to_find)
    assert subseq_ix == [
        slice(3, 7, None), slice(5, 9, None), slice(11, 15, None)
    ]
    subseqs = [seq_to_search[ix] for ix in subseq_ix]
    assert subseqs == ['atAT', 'ATaT', 'atat']
