from __future__ import annotations

from abc import ABC
from collections import abc
from typing import NewType, Union, List


__all__ = [
    # interface classes
    "ILocation", "IPart", "IPolymer",
    # types
    "LocationType", "SeqType"
]


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


class IPolymer(abc.MutableSequence):
    def __init__(self, seq: SeqType = '') -> None:
        self.seq = self._seq_check(seq)

    def _seq_check(self, value: SeqType) -> str:
        raise NotImplementedError

    def alphabet(self) -> List[str]:
        raise NotImplementedError


SeqType = NewType("SeqType", Union[str, IPolymer])
