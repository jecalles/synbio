from collections.abc import Generator
from synbio.annotations import Part
from synbio.polymers import Polymer, DNA


class Library(Part):
    def __init__(self, name=None, base_seq="", SeqType=DNA, variance=None,
                 location=None, metadata={},
                 ):
        # Type Check: raise error if seq_type is not a polymer
        if not issubclass(SeqType, Polymer):
            raise TypeError(
                f"SeqType must be a subclass of Polymer (type(SeqType) == {type(SeqType)})")
        # TODO: handle variance instantiation

            # convert input seq to Polymer if necessary
        if not isinstance(base_seq, Polymer):
            base_seq = SeqType(base_seq)

        super().__init__(name=name, seq=base_seq, location=location,
                         metadata=metadata)
        self.variance = variance
