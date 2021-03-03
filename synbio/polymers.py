from __future__ import annotations
from typing import Dict, List, Optional

from synbio import utils
from synbio.codes import Code, CodeType
from synbio.interfaces import ILocation, IPart, IPolymer, LocationType, SeqType

__all__ = [
    "Polymer", "NucleicAcid", "DNA", "RNA", "Protein"
]


class Polymer(IPolymer, utils.ComparableMixin):
    """
    An abstract base class from which NucleicAcid and Protein inherit.

    In order to create custom classes that inherit from Polymer, implement
    self.alphabet() (see Polymer.alphabet() docstring for details).
    """
    _comparables = ['seq']

    def __init__(self, seq: SeqType = '') -> None:
        self.seq = str(self._seq_check(seq))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.seq})"

    def __str__(self) -> str:
        return self.seq

    def __eq__(self, other: SeqType) -> bool:
        return self.seq.upper() == other

    def __len__(self) -> int:
        return len(self.seq)

    def __hash__(self) -> int:
        comparables = (
            val if "__hash__" in dir(val) else str(val)
            for val in vars(self).values()
        )
        return hash(comparables)

    def __getitem__(self, key: LocationType) -> "Own Type":
        return self.__class__(self.seq[key])

    def __setitem__(self, key: LocationType, value: SeqType) -> None:
        seqlist = list(self.seq)
        seqlist.__setitem__(key, self._seq_check(value))
        self.seq = ''.join(seqlist)

    def __delitem__(self, key: LocationType) -> None:
        seqlist = list(self.seq)
        seqlist.__delitem__(key)
        self.seq = ''.join(seqlist)

    def insert(self, key: int, value: SeqType) -> None:
        seqlist = list(self.seq)
        seqlist.insert(key, self._seq_check(value))
        self.seq = ''.join(seqlist)

    def _seq_check(self, value: SeqType) -> str:
        """
        A private method used to check that the characters of a given input
        sequence "value" all belong to a given Polymer's alphabet (defined in
        the alphbabet() method).

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
        # strip seq of whitespace
        seq = ''.join(str(value).split())
        # check if input sequence has appropriate alphabet
        bool_array = [item in self.alphabet() for item in seq.upper()]
        if not all(bool_array):
            raise ValueError(
                f"input value not in {self.__class__.__name__} alphabet")

        return seq

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
        >>> XNA("AATTCCGG") # raises ValueError: input value not in XNA alphabet

        Returns
        -------
        alphabet: list of valid characters
        """
        raise NotImplementedError


class NucleicAcid(Polymer):
    """
    """
    # TODO: write docstring
    _comparables = ['seq', 'annotations']
    basepairing = None

    def __init__(self,
                 seq: SeqType = "",
                 annotations: Optional[Dict[str, IPart]] = None) -> None:
        if annotations is None:
            annotations = {}
        else:
            try:
                annotations = dict(annotations)
            except TypeError:
                raise TypeError("annotations must be castable to dict")

            # raise TypeError if not all annotations are Parts
            if not all(
                    isinstance(part, IPart)
                    for part in annotations.values()
            ):
                raise TypeError("annotations must be of type Part")

        super().__init__(seq)
        self.annotations = annotations

    def __getitem__(self, key: str | LocationType) -> "Own Type":
        if isinstance(key, str):
            slice_ = self.annotations[key].location.to_slice()

        elif isinstance(key, ILocation):
            slice_ = key.to_slice()

            # shortcircuit - if REV strand, return rev comp
            if key.strand == "REV":
                return utils.reverse_complement(
                    self.seq.__getitem__(slice_), self.basepairing
                )
        else:
            slice_ = key

        return super().__getitem__(slice_)

    def __setitem__(self, key: str | LocationType, value: SeqType) -> None:
        length_change = len(value) - len(self[key])
        value = self._seq_check(value)

        if isinstance(key, str):
            slice_ = self.annotations[key].location.to_slice()
        elif isinstance(key, ILocation):
            slice_ = key.to_slice()

            if key.strand == "REV":
                value = utils.reverse_complement(
                    value, self.basepairing
                )
        else:
            slice_ = key

        super().__setitem__(slice_, value)
        self.update_annotations(slice_, length_change)

    def __delitem__(self, key: str | LocationType) -> None:
        length_change = -len(self.__getitem__(key))

        if isinstance(key, str):
            slice_ = self.annotations[key].location.to_slice()
        elif isinstance(key, ILocation):
            slice_ = key.to_slice()
        else:
            slice_ = key

        super().__delitem__(slice_)
        self.update_annotations(slice_, length_change)

    def insert(self, key: int, value: SeqType) -> None:
        """
        Please, don't use this method. Mkay? There are more idiomatic ways to
        work with NucleicAcids.

        Note: key may only be an int representing the index at which value
        will be added
        """
        if not isinstance(key, int):
            raise TypeError("key must be of type int")

        length_change = len(value)
        value = self._seq_check(value)

        super().insert(key, value)
        self.update_annotations(key, length_change)

    def update_annotations(self, key: LocationType, length_change: int) -> None:
        if self.annotations is not None:
            for part in self.annotations.values():
                part.update_location(key, length_change)

    def transcribe(self) -> RNA:
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError

    def reverse_transcribe(self) -> DNA:
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError

    def translate(self, code: Optional[CodeType] = None) -> Protein:
        """
        # TODO: write docstring
        Returns
        -------

        """
        if code is None:
            code = Code()
        elif isinstance(code, dict):
            code = Code(code)
        else:
            raise TypeError("code must be a dict or dict-like obj")

        mRNA = self.transcribe()
        prot_seq = code.translate(mRNA.seq)
        return Protein(prot_seq)

    def reverse_complement(self) -> "Own Type":
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError


class DNA(NucleicAcid):
    basepairing = utils.dna_basepair_WC

    def alphabet(self) -> List[str]:
        return utils.dNTPs

    def transcribe(self) -> RNA:
        return RNA(self.seq.replace('T', 'U').replace('t', 'u'))

    def reverse_transcribe(self) -> DNA:
        return self

    def reverse_complement(self) -> DNA:
        return DNA(
            utils.reverse_complement(self, utils.dna_basepair_WC)
        )


class RNA(NucleicAcid):
    basepairing = utils.rna_basepair_WC

    def alphabet(self) -> List[str]:
        return utils.rNTPs

    def transcribe(self) -> RNA:
        return self

    def reverse_transcribe(self) -> DNA:
        return DNA(self.seq.replace('U', 'T').replace('u', 't'))

    def reverse_complement(self) -> RNA:
        return RNA(
            utils.reverse_complement(self, utils.rna_basepair_WC)
        )


class Protein(Polymer):
    def alphabet(self) -> List[str]:
        return utils.aminoacids


if __name__ == "__main__":
    pass
