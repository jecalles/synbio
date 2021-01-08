from functools import partial
from typing import Iterable

from synbio import utils
from synbio.wrappers import codesavvy
from synbio.annotations import Part
from synbio.polymers import DNA, Polymer
from synbio.codes import Code
from synbio.codes.functions import get_synonymous_codons


class Library(Part):
    """
    """

    # TODO: write docstring

    def __init__(self, base_seq="", location=None, name=None, SeqType=DNA,
                 variance=lambda seq: seq, variance_params=None, metadata=None,
                 ):
        # Type Check: raise error if seq_type is not a polymer
        if not issubclass(SeqType, Polymer):
            raise TypeError("SeqType must be a subclass of Polymer")

        # convert input seq to Polymer if necessary
        if not isinstance(base_seq, Polymer):
            base_seq = SeqType(base_seq)

        # set default value for variance_params if None
        if variance_params is None:
            variance_params = dict()

        super().__init__(name=name, seq=base_seq, location=location,
                         metadata=metadata)

        # handle variance typing: expect functions -> iterators and iterators
        if callable(variance):
            variance_partial = partial(variance, **variance_params)
            variance = variance_partial(self.seq)
        if not isinstance(variance, Iterable):
            raise TypeError(
                "variance must be an Iterable, or a Callable that accepts a "
                "sequence as its first variable and returns an Iterable"
            )
        self.variance = variance

    def __iter__(self):
        yield from self.variance
#######################
# variance generators #
#######################

# TODO: write docstrings and tests for all of these


@codesavvy
def single_synonymous(seq, code=None):
    # get codons from transcript
    mRNA = DNA(seq).transcribe()
    codon_list = utils.get_codons(mRNA)
    # define a generator that returns all single synonymous_variants
    variants = (
        ''.join(codon_list[:pos] + list(synonym) + codon_list[pos+1:])
        for pos, codon in enumerate(codon_list)
        for synonym in get_synonymous_codons(codon)
    )
    return variants


@codesavvy
def double_synonymous(seq, code=None):
    # get codons from transcript
    mRNA = DNA(seq).transcribe()
    codon_list = utils.get_codons(mRNA)
    # define function that gets all synonymous codons for a given codon
    def get_synonymous_codons(codon): return code.rmap()[code[codon]]
    # define a generator that returns double synonymous_variants
    variants = (
        ''.join(
            codon_list[:ix1] + list(syn1)
            + codon_list[ix1+1:ix1+ix2+1] + list(syn2)
            + codon_list[ix1+ix2+2:]
        )
        for ix1, codon1 in enumerate(codon_list)
        for ix2, codon2 in enumerate(codon_list[ix1+1:])
        for syn1 in get_synonymous_codons(codon1)
        for syn2 in get_synonymous_codons(codon2)
    )
    return variants


@codesavvy
def single_nonsynonymous(seq, code=None, codon_frequencies=None):
    pass
