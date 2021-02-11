import itertools
from uuid import uuid4

from synbio.interfaces import ILocation, IPart
from synbio.polymers import DNA
from synbio.utils import ComparableMixin


class Location(ILocation, ComparableMixin):
    # TODO: write docstring
    """
    """

    _comparables = ['start', 'end', 'strand']

    def __init__(self, start, end, strand="FWD"):
        # check that start and end positions are valid
        if not start <= end:
            raise ValueError(
                f"start position ({start}) comes after end position ({end})")

        # check that strand input is valid
        strand = strand.upper()
        if strand not in {'FWD', 'REV'}:
            raise ValueError(f"input strand ({strand}) not FWD or REV")

        self.start = start
        self.end = end
        self.strand = strand

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start}, " \
               f"{self.end}, {self.strand})"

    @staticmethod
    def contains(outer_loc, inner_loc):
        """
        A static method that compares two locations, determining if the
        inner_loc is completely contained by the outer location.
        """
        return (
            outer_loc.start <= inner_loc.start
            and inner_loc.end <= outer_loc.end
        )

    @staticmethod
    def overlaps(loc1, loc2):
        """
        # TODO: write docstring for this method
        """
        return (loc1.start <= loc2.start < loc1.end) or (
            loc2.start <= loc1.start < loc2.end)

    @staticmethod
    def find_overlaps(locations):
        """
        A static method that, given a list of Location objects, constructs a
        graph representing all Location objects that overlap with each other.

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
        """
        # calculate overlaps
        loc_pairs = itertools.product(locations, repeat=2)
        bool_array = [
            Location.overlaps(loc1, loc2) for
            (loc1, loc2) in loc_pairs
        ]
        # reshape array to square matrix
        n = len(locations)
        bool_matrix = [
            bool_array[i * n: (i + 1) * n] for i in range(n)
        ]
        return bool_matrix

    @classmethod
    def from_slice(cls, slice_):
        start = 0 if slice_.start is None else slice_.start
        end = slice_.stop
        return cls(start, end)

    def to_slice(self):
        """
        """
        # TODO: write docstring
        return slice(self.start, self.end, 1)


class Part(IPart, ComparableMixin):
    # TODO: write docstring
    """
    """

    _comparables = [
        '_seq_id', 'location', 'name', 'kind', 'metadata'
    ]

    def __init__(self, seq=None, location=None, name=None, kind=None,
                 metadata=None):

        # default parameter initializations
        if seq is None:
            seq = DNA()

        if location is None:
            location = Location(0, len(seq))

        if name is None:
            name = str(uuid4())

        if kind is None:
            kind = self.__class__.__name__

        if metadata is None:
            metadata = dict()

        self._seq_reference = seq
        self._seq_id = id(seq)
        self.location = location
        self.name = str(name)
        self.kind = kind
        self.metadata = metadata

        # assign self as annotation to seq; fails if seq is str (that's okay)
        try:
            seq.annotations[self.name] = self
        except AttributeError:
            pass

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name}, " \
               f"{self.kind}, {self.location})"

    def __len__(self):
        return len(self.seq)

    @property
    def seq(self):
        try:
            # prefer slicing directly with location, if seq is strand savvy
            return self._seq_reference[self.location]
        except TypeError:
            # default to using .to_slice()
            return self._seq_reference[self.location.to_slice()]

    @seq.setter
    def seq(self, value):
        try:
            # prefer slicing directly with location, if seq is strand savvy
            self._seq_reference.__setitem__(self.location, value)
        except AttributeError:
            # default to using .to_slice()
            self._seq_reference.__setitem__(self.location.to_slice(),
                                            value)

    def update_location(self, key, length_change):
        # handle parsing of key into Location object
        if isinstance(key, ILocation):
            update_loc = key

        elif isinstance(key, slice):
            update_loc = Location.from_slice(key)

        elif isinstance(key, int):
            update_loc = Location(key, key + 1)

        else:
            raise TypeError("could not convert key into Location")

        # if update is fully upstream of Part, update both start and end
        if update_loc.end <= self.location.start:
            self.location.start += length_change
            self.location.end += length_change
        # if update is contained within Part, then only update end
        elif Location.contains(self.location, update_loc):
            self.location.end += length_change
        # if update overlaps with one end of Part, truncate Part
        elif Location.overlaps(self.location, update_loc) \
                and length_change != 0:
            if update_loc.start < self.location.start:
                self.location.end = update_loc.start
            else:
                self.location.start = update_loc.end
