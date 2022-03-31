from __future__ import annotations

from abc import ABC, abstractmethod
from collections import abc
from functools import wraps
from typing import List, Union

__all__ = [
    # wrappers
    "compare_same_type",
    # mixins
    "ComparableMixin",
    # annotations
    "ILocation", "IPart", "IPolymer", "LocationType", "IndexType", "SeqType",
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


################
#  Interfaces  #
################
class IPolymer(abc.MutableSequence):
    @abstractmethod
    def __hash__(self):
        raise NotImplementedError

    @abstractmethod
    def _seq_check(self, value: SeqType) -> str:
        """
        A private method used to check that the characters of a given input
        sequence "value" all belong to a given Polymer's alphabet (defined in
        the alphabet() method).

        Parameters
        ----------
        value: sequence to check

        Returns
        -------
        True
            if all characters are part of self.alphabet(), and
        False
            otherwise
        """
        raise NotImplementedError

    @abstractmethod
    def alphabet(self) -> List[str]:
        """
        A function prototype that returns a list of characters
        representing valid inputs for a given Polymer subclass. Implement this
        method in order to inherit from Polymer.

        E.g.,

        >>> class XNA(Polymer):
        >>>     def alphabet(self):
        >>>         return ["X", "Y", "W", "Z"]

        >>> XNA("XXYYZZ")
        XNA(XXYYZZ)
        >>> XNA("AATTCCGG")#raises ValueError("input value not in XNA alphabet")
        """
        raise NotImplementedError


SeqType = Union[str, IPolymer]


# annotations
class ILocation(ABC, ComparableMixin):
    @abstractmethod
    def __init__(self):
        self.start = None
        self.end = None
        self.strand = None

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


LocationType = Union[int, slice, ILocation, List[ILocation]]
IndexType = Union[str, LocationType]


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
