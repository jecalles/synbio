from synbio.utils import Location
from synbio.polymers import NucleicAcid, DNA


class Part:
    # TODO: write docstring
    '''
    '''

    def __init__(self, name=None, kind=None, seq="", location=None, subparts=None, parent=None, metadata={}):
        # try to cast seq as a DNA, if not already a nucleic acid
        if not isinstance(seq, NucleicAcid):
            seq = DNA(seq)

        # check that parent is of type Part
        if parent is not None and not isinstance(parent, Part):
            raise TypeError("parent must be of type Part")

        # initialize location if appropriate
        if location is None:
            if parent is not None:
                raise ValueError("location must be provided if Part is a subpart")
            else:
                location = Location(0, len(seq))

        # initialize kind if None
        if kind is None:
            kind = self.__class__.__name__

        if subparts is not None:
            # check that subparts don't extend beyond this Part
            location_bool_array = [
                Location.contains(location, sub.location) for sub in subparts
            ]

            if not all(location_bool_array):
                overflowing_subparts = list(
                    itertools.compress(subparts, location_bool_array)
                )
                raise ValueError(
                    f"The following subparts are located outside the range of this part {overflowing_subparts}"
                )

            # update subpart.parent field for each child
            for sub in subparts:
                sub.parent = self

        # TODO: are all these attributes necessary?
        self.name = name
        self.kind = kind
        self._seq_reference = seq
        self.location = location
        self.parent = parent
        self.subparts = subparts
        self.metadata = metadata

    def __add__(self, other):
        # TODO: implement Part concatenation
        raise NotImplementedError

    @property
    def seq(self):
        return self._seq_reference[self.location]

    @seq.setter
    def seq(self, value):
        prev_length = self.location.end - self.location.start
        new_length = len(value)
        self._seq_reference.__setitem__(self.location, value)
        self.location.end += (new_length - prev_length)
