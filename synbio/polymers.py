from abc import ABCMeta
from collections import abc
from copy import copy

from synbio import utils
from synbio.annotations import Location, Part
from synbio.codes import Code


class Polymer(abc.MutableSequence, utils.ComparableMixin):
    """
    # TODO: write docstring
    """
    _comparables = ['seq']

    def __init__(self, seq=''):
        self.seq = self._seq_check(seq)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.seq})"

    def __str__(self):
        return self.seq

    def __eq__(self, other):
        return self.seq == other

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.__class__(self.seq[key])

    def __setitem__(self, key, value):
        seqlist = list(self.seq)
        seqlist.__setitem__(key, self._seq_check(value))
        self.seq = ''.join(seqlist)

    def __delitem__(self, key):
        seqlist = list(self.seq)
        seqlist.__delitem__(key)
        self.seq = ''.join(seqlist)

    def insert(self, key, value):
        seqlist = list(self.seq)
        seqlist.insert(key, self._seq_check(value))
        self.seq = ''.join(seqlist)

    def _seq_check(self, value):
        # handle different input types
        if isinstance(value, Polymer):
            seq = value.seq
        else:
            seq = str(value)
        # check if input sequence has appropriate alphabet
        bool_array = [item in self.alphabet() for item in seq.upper()]
        if not all(bool_array):
            raise ValueError(
                f"input value not in {self.__class__.__name__} alphabet")
        return seq

    def alphabet(self):
        raise NotImplementedError


class NucleicAcid(Polymer, metaclass=ABCMeta):
    """
    """
    # TODO: write docstring
    _comparables = ['seq', 'annotations']
    basepairing = None

    def __init__(self, seq, annotations=None):
        """
        Parameters
        ----------
        seq
        annotations
        """
        # TODO: write docstring

        if annotations is None:
            # default value is empty dict
            annotations = {}
        else:
            # try to cast annotations as dict
            if isinstance(annotations, dict):
                pass
            else:
                try:
                    annotations = {
                        part.name: part
                        for part in annotations
                    }
                except TypeError:
                    raise TypeError("annotations must be castable to dict")

            # raise TypeError if not all annotations are Parts
            if not all(
                isinstance(part, Part)
                for part in annotations.values()
            ):
                raise TypeError("annotation must be a Part")

        super().__init__(seq)
        self.annotations = annotations

    def __getitem__(self, key):
        if isinstance(key, Location):
            slice_ = key.to_slice()

            # shortcircuit - if REV strand, return rev comp
            if key.strand == "REV":
                return utils.reverse_complement(
                    self.seq.__getitem__(slice_), self.basepairing
                )
        else:
            slice_ = key

        return super().__getitem__(slice_)

    def __setitem__(self, key, value):
        length_change = len(value) - len(self[key])

        value = self._seq_check(value)
        if isinstance(key, Location):
            loc = copy(key)
            key = key.to_slice()

            if loc.strand == "REV":
                value = utils.reverse_complement(
                    value, self.basepairing
                )

        elif isinstance(key, slice):
            loc = Location.from_slice(key)

        elif isinstance(key, int):
            loc = Location(key, key + 1)

        else:
            raise TypeError("could not convert key into Location")

        super().__setitem__(key, value)

        if self.annotations is not None:
            for part in self.annotations.values():
                part.update_location(loc, length_change)

    def __delitem__(self, key):
        length_change = -len(self.__getitem__(key))

        if isinstance(key, Location):
            loc = copy(key)
            key = key.to_slice()

        elif isinstance(key, slice):
            loc = Location.from_slice(key)

        elif isinstance(key, int):
            loc = Location(key, key + 1)

        else:
            raise TypeError("could not convert key into Location")

        super().__delitem__(key)

        if self.annotations is not None:
            for part in self.annotations.values():
                part.update_location(loc, length_change)

    def insert(self, key, value):
        """
        Please, don't use this method. Mkay? There are more idiomatic ways to
        work with NucleicAcids.

        Note: key may only be an int representing the index at which value
        will be added
        """
        # TODO: write a better docstring
        if not isinstance(key, int):
            raise TypeError("key must be of type int")

        length_change = len(value)
        value = self._seq_check(value)
        loc = Location(key, key + 1)

        super().insert(key, value)

        if self.annotations is not None:
            for part in self.annotations.values():
                part.update_location(loc, length_change)

    def transcribe(self):
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError

    def reverse_transcribe(self):
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError

    def translate(self):
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError

    def reverse_complement(self):
        """
        # TODO: write docstring
        Returns
        -------

        """
        raise NotImplementedError


class DNA(NucleicAcid):
    basepairing = utils.dna_basepair_WC

    def alphabet(self):
        return utils.dNTPs

    def transcribe(self):
        return RNA(self.seq.replace('T', 'U'))

    def reverse_transcribe(self):
        return self

    def translate(self, code=None):
        if code is None:
            code = Code()
        elif isinstance(code, dict):
            code = Code(code)
        else:
            raise TypeError("code must be a dict or dict-like obj")
        mRNA = self.transcribe()
        return mRNA.translate(code)

    def reverse_complement(self):
        return DNA(
            utils.reverse_complement(self, utils.dna_basepair_WC)
        )


class RNA(NucleicAcid):
    basepairing = utils.rna_basepair_WC

    def alphabet(self):
        return utils.rNTPs

    def transcribe(self):
        return self

    def reverse_transcribe(self):
        return DNA(self.seq.replace('U', 'T'))

    def translate(self, code=None):
        if code is None:
            code = Code()
        elif isinstance(code, dict):
            code = Code(code)
        else:
            raise TypeError("code must be a dict or dict-like obj")
        prot_seq = code.translate(self.seq)
        return Protein(prot_seq)

    def reverse_complement(self):
        return RNA(
            utils.reverse_complement(self, utils.rna_basepair_WC)
        )


class Protein(Polymer):
    def alphabet(self):
        return utils.aminoacids
