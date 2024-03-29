from functools import partial
from typing import Any, Callable, Dict, Iterable, Optional

from synbio import utils
from synbio.annotations import Location, Part
from synbio.codes import CodeType
from synbio.codes.functions import get_synonymous_codons
from synbio.codes.wrappers import codesavvy
from synbio.interfaces import SeqType
from synbio.polymers import DNA, Polymer

__all__ = [
    # class
    "Library",
    # variance generators
    "single_synonymous", "single_nonsynonymous",
    "double_synonymous", "double_nonsynonymous"
]


class Library(Part):
    """
    TODO: write docstring
    """

    def __init__(
            self,
            base_seq: Optional[SeqType] = None,
            location: Optional[Location] = None,
            name: Optional[str] = None,
            seq_type: type(Polymer) = DNA,
            variance: Callable[[SeqType], Iterable[SeqType]] = lambda s: [s],
            variance_params: Optional[dict] = None,
            metadata: Optional[Dict[str, Any]] = None
    ):
        # Type Check: raise error if seq_type is not a Polymer
        if not issubclass(seq_type, Polymer):
            raise TypeError("SeqType must be a subclass of Polymer")

        if base_seq is None:
            seq = seq_type()
        elif isinstance(base_seq, Polymer):
            seq = base_seq
        else:
            seq = seq_type(base_seq)

        if variance_params is None:
            variance_params = dict()

        super().__init__(name=name, seq=seq, location=location,
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
def single_synonymous(
        seq: SeqType,
        code: Optional[CodeType] = None) -> Iterable[str]:
    # get codons from transcript
    mRNA = DNA(seq).transcribe()
    codon_list = utils.get_codons(mRNA)
    # define a generator that returns all single synonymous_variants
    variants = (
        ''.join(codon_list[:pos] + list(synonym) + codon_list[pos + 1:])
        for pos, codon in enumerate(codon_list)
        for synonym in get_synonymous_codons(codon, code)
    )
    return variants


@codesavvy
def double_synonymous(
        seq: SeqType,
        code: Optional[CodeType] = None) -> Iterable[str]:
    # get codons from transcript
    mRNA = DNA(seq).transcribe()
    codon_list = utils.get_codons(mRNA)
    # define a generator that returns double synonymous_variants
    variants = (
        ''.join(
            codon_list[:ix1] + list(syn1)
            + codon_list[ix1 + 1:ix1 + ix2 + 1] + list(syn2)
            + codon_list[ix1 + ix2 + 2:]
        )
        for ix1, codon1 in enumerate(codon_list)
        for ix2, codon2 in enumerate(codon_list[ix1 + 1:])
        for syn1 in get_synonymous_codons(codon1, code)
        for syn2 in get_synonymous_codons(codon2, code)
    )
    return variants


@codesavvy
def single_nonsynonymous(
        seq: SeqType,
        code: Optional[CodeType] = None,
        codon_frequencies: Optional[Dict[str, float]] = None):
    # TODO: implement
    raise NotImplementedError("function under construction")


@codesavvy
def double_nonsynonymous():
    # TODO: implement
    raise NotImplementedError("function under construction")
