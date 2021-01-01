from collections import abc
from synbio import utils
from synbio.codes import Code


class Polymer(abc.MutableSequence):

    def __init__(self, seq=''):
        self.seq = self._seq_check(seq)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.seq})"

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

    def transcribe(self):
        raise NotImplementedError

    def reverse_transcribe(self):
        raise NotImplementedError

    def translate(self):
        raise NotImplementedError


class DNA(NucleicAcid):
    def alphabet(self):
        return utils.dNTPs

    def transcribe(self):
        return RNA(self.seq.replace('T', 'U'))

    def reverse_transcribe(self):
        return self

    def translate(self, code=Code()):
        mRNA = self.transcribe()
        return mRNA.translate(code)


class RNA(NucleicAcid):
    def alphabet(self):
        return utils.rNTPs

    def transcribe(self):
        return self

    def reverse_transcribe(self):
        return DNA(self.seq.replace('U', 'T'))

    def translate(self, code=Code()):
        return Code(code).translate(self)


class Protein(Polymer):
    def alphabet(self):
        return utils.aminoacids
