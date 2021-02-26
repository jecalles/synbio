from abc import ABC
from typing import NewType, Union


class ILocation(ABC):
    def __init__(self):
        start = None
        end = None
        strand = None

    @staticmethod
    def contains(outer_loc, inner_loc):
        raise NotImplementedError

    @staticmethod
    def overlaps(loc1, loc2):
        raise NotImplementedError

    @staticmethod
    def find_overlaps(locations):
        raise NotImplementedError

    @classmethod
    def from_slice(cls, slice_):
        raise NotImplementedError

    def to_slice(self):
        raise NotImplementedError

LocationType = NewType("LocationType", Union[int, slice, ILocation])

class IPart(ABC):
    def __init__(self):
        self._seq_reference = None
        self._seq_id = None
        self.location = None
        self.name = None
        self.kind = None
        self.metadata = None

    @property
    def seq(self):
        raise NotImplementedError

    @seq.setter
    def seq(self, value):
        raise NotImplementedError

    def update_location(self, key, length_change):
        raise NotImplementedError
