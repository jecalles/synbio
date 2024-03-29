from __future__ import annotations

import itertools
from functools import reduce
from typing import Any, Dict, List, Optional, Sequence
from uuid import uuid4

from synbio.interfaces import ILocation, IPart, LocationType, SeqType
from synbio.polymers import DNA

__all__ = [
    "Location", "Part"
]


class Location(ILocation):
    """
    A class used to specify a location on a SeqType, as well as a
    strand--FWD or REV for NucleicAcids, FWD by default for all other types.
    """

    def __init__(self, start: int, end: int, strand: str = "FWD") -> None:
        # check that start and end positions are valid
        if not start <= end:
            raise ValueError(
                f"start position ({start}) comes after end position ({end})")

        # check that strand input is valid
        strand = strand.upper()
        if strand not in {'FWD', 'REV'}:
            raise ValueError(f"input strand ({strand}) not FWD or REV")

        self.start = int(start)
        self.end = int(end)
        self.strand = str(strand)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.start}, " \
               f"{self.end}, {self.strand})"

    def _comparables(self) -> List[str]:
        return ['start', 'end', 'strand']

    @staticmethod
    def contains(outer_loc: Location, inner_loc: Location) -> bool:
        """
        A static method that compares two locations, determining if the
        inner_loc is completely contained by the outer location.

        E.g.,
        >>> a = Location(0, 10) # a: ----------
        >>> b = Location(4, 7)  # b:     ----
        >>> c = Location(5, 15) # c:      ----------
        >>> Location.contains(a, a)
        True
        >>> Location.contains(a, b)
        True
        >>> Location.contains(b, a)
        False
        >>> Location.contains(a, c)
        False
        """
        return (
                outer_loc.start <= inner_loc.start
                and inner_loc.end <= outer_loc.end
        )

    @staticmethod
    def overlaps(loc1: Location, loc2: Location) -> bool:
        """
        A static method that, given two Location objects, returns True if the
        two locations overlap and False otherwise

        E.g.,

        >>> a = Location(0, 10) # a: ----------
        >>> b = Location(4, 7)  # b:     ----
        >>> c = Location(0, 4)  # c: ----
        >>> Location.overlaps(a, a)
        True
        >>> Location.overlaps(a, b)
        True
        >>> Location.overlaps(a, c)
        True
        >>> Location.overlaps(b, c)
        False
        """
        return (loc1.start <= loc2.start < loc1.end) or (
                loc2.start <= loc1.start < loc2.end)

    @staticmethod
    def find_overlaps(locations: Sequence[Location]) -> List[List[bool]]:
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

    def offset(self, offset: int) -> Location:
        """
        A method that returns a new Location object whose start and end
        positions are offset by an integer value

        >>> loc1 = Location(10, 20, "FWD")
        >>> loc2 = loc1.offset(5)
        >>> print(loc1, loc2)
        Location(10, 20, FWD) Location(15, 25, FWD)
        >>> loc3 = Location(15, 25, "REV")
        >>> loc4 = loc3.offset(5)
        >>> print(loc3, loc4)
        Location(15, 25, REV) Location(20, 30, REV)
        """
        # check type of offset
        if not isinstance(offset, int):
            raise TypeError("offset must be of type int")

        return Location(
            start=self.start + offset,
            end=self.end + offset,
            strand=self.strand
        )

    @classmethod
    def from_slice(cls, slice_: slice) -> Location:
        """
        A method that returns a new Location object from a slice object
        """
        start = 0 if slice_.start is None else slice_.start
        end = slice_.stop
        return cls(start, end)

    def to_slice(self) -> slice:
        """
        A method that returns a slice object representing the Location obj
        """
        return slice(self.start, self.end, 1)


class Part(IPart):
    """
    A class used to represent an annotation of or an abstraction over
    NucleicAcids (e.g., promoters, genes, terminators, restriction sites, etc.)
    """

    def __init__(
            self,
            seq: Optional[SeqType] = None,
            location: Optional[LocationType] = None,
            name: Optional[str] = None,
            kind: Optional[str] = None,
            metadata: Optional[Dict[str, Any]] = None
    ) -> None:

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

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.name}, " \
               f"{self.kind}, {self.location})"

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, key: LocationType) -> Part:
        if isinstance(key, ILocation):
            loc = key
        elif isinstance(key, slice):
            loc = Location.from_slice(key)
        else:
            raise TypeError("key must be of type ILocation or slice")

        loc.strand = self.location.strand

        return Part(
            seq=self._seq_reference,
            location=loc.offset(self.location.start),
            name=f"{self.name}_subset",
            kind=self.kind,
            metadata=self.metadata
        )

    def _comparables(self) -> List[str]:
        return [
            '_seq_id', 'location', 'name', 'kind', 'metadata'
        ]

    def _seq_index(self, loc: LocationType):
        try:
            return self._seq_reference[loc]
        except TypeError:
            return self._seq_reference[loc.to_slice()]

    @property
    def seq(self) -> SeqType:
        if isinstance(self.location, list):
            return reduce(
                lambda x, y: x + y,
                [self._seq_index(loc) for loc in self.location]
            )
        else:
            return self._seq_index(self.location)

    @seq.setter
    def seq(self, value: SeqType) -> None:
        try:
            # prefer slicing directly with location, if seq is strand savvy
            self._seq_reference.__setitem__(self.location, value)
        except AttributeError:
            # default to using .to_slice()
            self._seq_reference.__setitem__(self.location.to_slice(),
                                            value)

    def update_location(self, key: LocationType, length_change: int) -> None:
        # handle parsing of key into Location object
        if isinstance(key, Location):
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
