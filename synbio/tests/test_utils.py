# TODO: write test code for utils functions
from synbio.utils import *
import synbio.tests.utils as testutils

__all__ = [
    #############
    # functions #
    #############
    # biology stuff
    "get_codons", "reverse_complement", "is_palindrome", "find_subseq",
    "all_single_mutations", "mutation_pairs",
    # python stuff
    "get_class_name"
]

def test_get_codons():
    raise testutils.TestNotImplemented

def test_reverse_complement():
    raise testutils.TestNotImplemented

def test_is_palindrome():
    raise testutils.TestNotImplemented

def test_find_subseq():
    seq_to_search = "xxXatATaTxXatatxx"
    subseq_to_find = "atat"  # search is NOT case sensitive

    subseq_ix = find_subseq(seq_to_search, subseq_to_find)
    assert subseq_ix == [
        slice(3, 7, None), slice(5, 9, None), slice(11, 15, None)
    ]
    subseqs = [seq_to_search[ix] for ix in subseq_ix]
    assert subseqs == ['atAT', 'ATaT', 'atat']

def test_all_single_mutations():
    raise testutils.TestNotImplemented

def test_mutation_pairs():
    raise testutils.TestNotImplemented

def test_get_class_name():
    raise testutils.TestNotImplemented

