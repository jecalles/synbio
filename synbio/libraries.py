from typing import Iterable

from synbio.annotations import Part
from synbio.polymers import Polymer, DNA


class Library(Part):
    """
    """
    # TODO: write docstring

    def __init__(self, base_seq="", location=None, name=None, SeqType=DNA,
                 variance=lambda seq: seq, variance_params={}, metadata={},
                 ):
        # Type Check: raise error if seq_type is not a polymer
        if not issubclass(SeqType, Polymer):
            raise TypeError(
                f"SeqType must be a subclass of Polymer (type(SeqType) == {type(SeqType)})")

            # convert input seq to Polymer if necessary
        if not isinstance(base_seq, Polymer):
            base_seq = SeqType(base_seq)

        super().__init__(name=name, seq=base_seq, location=location,
                         metadata=metadata)

        # handle variance typing: expect functions -> iterators and iterators
        if callable(variance):
            variance = variance(self.seq, **variance_params)
        if not isinstance(variance, Iterable):
            raise TypeError(
                "variance must be an Iterable, or a Callable that accepts a sequence as its first variable and returns an Iterable")
        self.variance = variance

    def __iter__(self):
        yield from self.variance
