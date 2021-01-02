from collections import abc
from synbio import utils
from synbio.codes import Code
from synbio.annotations import Location


class Polymer(abc.MutableSequence):

    def __init__(self, seq=''):
        self.seq = self._seq_check(seq)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.seq})"

    def __eq__(self, other):
        return self.seq == other

    def __len__(self):
        return self.seq.__len__()

    def __getitem__(self, key):
        return self.seq.__getitem__(key)

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
            raise ValueError(f"input value not in {self.__class__.__name__} alphabet")
        return seq

    def alphabet(self):
        raise NotImplementedError


class NucleicAcid(Polymer):

    def __init__(self, seq, annotations=None):
        super().__init__(seq)
        self.annotations = annotations

    def __getitem__(self, key):
        # TODO: implement compatibility with Location objects
        if type(key) == Location:
            slice = key.to_slice()

            # shortcircuit - if REV strand, return rev comp
            if key.strand == "REV":
                return utils.reverse_complement(
                    self.seq.__getitem__(slice), self.basepairing
                )
        else:
            slice = key

        return super().__getitem__(slice)

    def __setitem__(self, key, value):
        value = self._seq_check(value)

        if type(key) == Location:
            slice = key.to_slice()

            if key.strand == "REV":
                value = utils.reverse_complement(
                    value, self.basepairing
                )

        else:
            slice = key

        return super().__setitem__(slice, value)

    def __delitem__(self, key):
        if type(key) == Location:
            slice = key.to_slice()

        return super().__delitem__(slice, value)

    def insert(self, key, value):
        value = self._seq_check(value)

        if type(key) == Location:
            slice = key.to_slice()

            if key.strand == "REV":
                value = utils.reverse_complement(
                    value, self.basepairing
                )

        else:
            slice = key

        return super().insert(slice, value)

    def transcribe(self):
        raise NotImplementedError

    def reverse_transcribe(self):
        raise NotImplementedError

    def translate(self):
        raise NotImplementedError

    def reverse_complement(self):
        raise NotImplementedError


class DNA(NucleicAcid):
    basepairing = utils.dna_basepair_WC

    def alphabet(self):
        return utils.dNTPs

    def transcribe(self):
        return RNA(self.seq.replace('T', 'U'))

    def reverse_transcribe(self):
        return self

    def translate(self, code=Code()):
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

    def translate(self, code=Code()):
        return Code(code).translate(self)

    def reverse_complement(self):
        return RNA(
            utils.reverse_complement(self, utils.rna_basepair_WC)
        )


class Protein(Polymer):
    def alphabet(self):
        return utils.aminoacids
