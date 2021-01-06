import itertools
from synbio.utils import compare_same_type


class Location:
    # TODO: write docstring
    '''
    '''

    def __init__(self, start, end, strand="FWD"):
        # check that start and end positions are valid
        if not start <= end:
            raise ValueError(f"start position ({start}) comes after end position ({end})")

        # check that strand input is valid
        strand = strand.upper()
        if strand not in {'FWD', 'REV'}:
            raise ValueError(f"input strand ({strand}) not FWD or REV")

        self.start = start
        self.end = end
        self.strand = strand

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start}, {self.end}, {self.strand})"

    @compare_same_type
    def __eq__(self, other):
        return self.start == other.start and self.end == other.end and self.strand == other.strand

    @staticmethod
    def contains(outer_loc, inner_loc):
        '''
        A static method that compares two locations, determining if the inner_loc is completely contained by the outer location.
        '''
        return (outer_loc.start <= inner_loc.start) and (inner_loc.end <= outer_loc.end)

    @staticmethod
    def overlaps(loc1, loc2):
        '''
        # TODO: write docstring for this method
        '''
        # TODO: implement lol
        return (loc1.start <= loc2.start and loc1.end > loc2.start) or (loc2.start <= loc1.start and loc2.end > loc1.start)

    @ staticmethod
    def find_overlaps(locations):
        '''
        A static method that, given a list of Location objects, constructs a graph representing all Location objects that overlap with each other.

        E.g., suppose our Locations are the following:


        a: ----------
        b:     ---
        c: ----
        d:           ----------
        e:      ----------

        >>> a = Location(0, 10)
        >>> b = Location(4, 7)
        >>> c = Location(0, 4)
        >>> d = Location(10, 20)
        >>> e = Location(5, 15)

        >>> Location.find_overlaps([a, b, c, d, e])
        [
            [1, 1, 1, 0, 1],
            [1, 1, 0, 0, 1],
            [1, 0, 1, 0, 0],
            [0, 0, 0, 1, 1],
            [1, 1, 0, 1, 1],
        ]
        '''
        # calculate overlaps
        loc_pairs = itertools.product(locations, repeat=2)
        bool_array = [
            Location.overlaps(loc1, loc2) for
            (loc1, loc2) in loc_pairs
        ]
        # reshape array to square matrix
        n = len(locations)
        bool_matrix = [
            bool_array[i*n: (i+1)*n] for i in range(n)
        ]
        return bool_matrix

    @classmethod
    def from_slice(cls, slice):
        start = 0 if slice.start is None else slice.start
        end = slice.stop
        return cls(start, end)

    def to_slice(self):
        '''
        '''
        # TODO: write docstring
        return slice(self.start, self.end, 1)


class Part:
    # TODO: write docstring
    '''
    '''

    def __init__(self, name=None, kind=None, seq="", location=None,  metadata={}):
        # assign self as annotation to seq; fails if seq is str (thats okay)
        try:
            seq.annotations.append(self)
        except:
            pass

        # initialize location if appropriate
        if location is None:
            location = Location(0, len(seq))

        # initialize kind if None
        if kind is None:
            kind = self.__class__.__name__

        # TODO: are all these attributes necessary?
        self.name = name
        self.kind = kind
        self._seq_reference = seq
        self.location = location
        self.metadata = metadata

    def __add__(self, other):
        # TODO: implement Part concatenation
        raise NotImplementedError

    @property
    def seq(self):
        try:
            # prefer slicing directly with location, if seq is strand savy
            return self._seq_reference[self.location]
        except:
            # default to using .to_slice()
            return self._seq_reference[self.location.to_slice()]

    @seq.setter
    def seq(self, value):
        length_change = len(value) - len(self.seq)

        try:
            # prefer slicing directly with location, if seq is strand savy
            self._seq_reference.__setitem__(self.location, value)
        except:
            # default to using .to_slice()
            self._seq_reference.__setitem__(self.location.to_slice(),
                                            value)

    def update_location(self, update_loc, length_change, skip_parent=False):
        # if update is fully upstream of Part, update both start and end
        if update_loc.end <= self.location.start:
            self.location.start += length_change
            self.location.end += length_change
        # if update is contained within Part, then only update end
        elif Location.contains(self.location, update_loc):
            self.location.end += length_change
        # if update overlaps with one end of Part, truncate Part
        elif Location.overlaps(self.location, update_loc) and length_change != 0:
            if update_loc.start < self.location.start:
                self.location.end = update_loc.start
            else:
                self.location.start = update_loc.end
