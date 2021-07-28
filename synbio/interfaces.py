from __future__ import annotations

from abc import ABC, abstractmethod
from collections import abc
from functools import wraps
from pathlib import Path
from typing import List, TypeVar

__all__ = [
    # wrappers
    "compare_same_type",
    # mixins
    "ComparableMixin",
    # annotations
    "ILocation", "IPart", "IPolymer", "LocationType", "IndexType", "SeqType",
    # # parsers
    # "IParser", "FileType"
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
    A Mixin class that automatically provides an equality comparison method.
    Requires concrete implementation of _comparables method
    """
    @abstractmethod
    def _comparables(self) -> List[str]:
        """
        Abstract method that must be implemented to subclass from
        ComparableMixin or subclasses. Should return a list of strings
        representing instance attribute names to be compared.

        E.g.,

        >>> class Foo(ComparableMixin):
        >>>     def __init__(self, a, b, c):
        >>>         self.val_a = a
        >>>         self.val_b = b
        >>>         self.val_c = c
        >>>     def _comparables(self):
        >>>         return ['val_a', 'val_b']

        >>> alice = Foo(1, 2, "not_compared")
        >>> bobby = Foo(1, 2, "just_ignored")
        >>> alice == bobby
        True
        """
        raise NotImplementedError

    @compare_same_type
    def __eq__(self, other: "Same Type") -> bool:
        return all(
            getattr(self, comp) == getattr(other, comp)
            for comp in self._comparables()
        )

class HashableMixin(ComparableMixin):
    """
    A Mixin class that extends ComparableMixin to provide both equality and
    hash methods, given a class attribute _comparables
    """
    def __hash__(self) -> int:
        return hash(
            (getattr(self, comp) for comp in self._comparables())
        )


################
#  Interfaces  #
################


# polymers
class IPolymer(abc.MutableSequence, ComparableMixin):
    def __init__(self, seq: SeqType = '') -> None:
        self.seq = self._seq_check(seq)

    @abstractmethod
    def _seq_check(self, value: SeqType) -> str:
        raise NotImplementedError

    @abstractmethod
    def alphabet(self) -> List[str]:
        raise NotImplementedError


SeqType = TypeVar("SeqType", str, IPolymer)


# annotations
class ILocation(ABC, ComparableMixin):
    @abstractmethod
    def __init__(self):
        start = None
        end = None
        strand = None

    @staticmethod
    @abstractmethod
    def contains(outer_loc, inner_loc):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def overlaps(loc1, loc2):
        raise NotImplementedError

    @staticmethod
    @abstractmethod
    def find_overlaps(locations):
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def from_slice(cls, slice_):
        raise NotImplementedError

    @abstractmethod
    def to_slice(self):
        raise NotImplementedError


LocationType = TypeVar("LocationType", int, slice, ILocation)
IndexType = TypeVar("IndexType", str, LocationType)


class IPart(ABC, ComparableMixin):
    @abstractmethod
    def __init__(self):
        self._seq_reference = None
        self._seq_id = None
        self.location = None
        self.name = None
        self.kind = None
        self.metadata = None

    @abstractmethod
    def __len__(self) -> int:
        raise NotImplementedError

    @abstractmethod
    def __getitem__(self, key: LocationType) -> "Own Type":
        raise NotImplementedError

    @property
    @abstractmethod
    def seq(self):
        raise NotImplementedError

    @seq.setter
    @abstractmethod
    def seq(self, value):
        raise NotImplementedError

    @abstractmethod
    def update_location(self, key, length_change):
        raise NotImplementedError


# parsers


class IParser(ABC):
    @abstractmethod
    def read(self, filename: FileType) -> str:
        pass

    @abstractmethod
    def write(self, obj: IWriteable, filename: FileType) -> None:
        pass

    @abstractmethod
    def write_string(self, obj: IWriteable) -> str:
        pass


FileType = TypeVar("FileType", str, Path)


class IWriteable:
    pass
