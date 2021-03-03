from __future__ import annotations

from abc import ABC
from collections import abc
from typing import NewType, Union, List
from functools import wraps

__all__ = [
    # wrappers
    "compare_same_type",
    # mixins
    "ComparableMixin",
    # interface classes
    "ILocation", "IPart", "IPolymer",
    # types
    "LocationType", "SeqType"
]


#############
# wrappers  #
#############
def compare_same_type(func):
    @wraps(func)
    def wrapper(self, other):
        if not issubclass(type(self), type(other)):
            raise TypeError(
                f"Cannot compare {self.__class__.__name__} with {type(other)}")
        else:
            return func(self, other)

    return wrapper


#############
#  Mixins   #
#############
class ComparableMixin:
    """
    A Mixin class that automatically provides some comparison dunder methods,
    given a class attribute _comparables
    """
    _comparables = list()

    def __hash__(self) -> int:
        return hash(
            (getattr(self, comp) for comp in self._comparables)
        )

    @compare_same_type
    def __eq__(self, other: "Same Type") -> bool:
        return all(
            getattr(self, comp) == getattr(other, comp)
            for comp in self._comparables
        )


################
#  Interfaces  #
################
class ILocation(ABC, ComparableMixin):
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


class IPart(ABC, ComparableMixin):
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


class IPolymer(abc.MutableSequence, ComparableMixin):
    def __init__(self, seq: SeqType = '') -> None:
        self.seq = self._seq_check(seq)

    def _seq_check(self, value: SeqType) -> str:
        raise NotImplementedError

    def alphabet(self) -> List[str]:
        raise NotImplementedError


SeqType = NewType("SeqType", Union[str, IPolymer])

