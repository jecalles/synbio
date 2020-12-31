from synbio import utils
from synbio.codes import Code


class Polymer(str):
    def __new__(cls, value):
        if not cls._alphabet_check(value):
            raise ValueError(f"input value not in {cls.__name__} alphabet")
        else:
            return super().__new__(cls, value.upper())

    def __repr__(self):
        # truncate representation if too long
        if len(self) > 160:
            seq = f"{self[:20]}...{str(self[-20:])}"
        else:
            seq = str(self)
        return f"{self.__class__.__name__}({seq})"

    @classmethod
    def _alphabet_check(cls, value):
        str_value = str(value).upper()
        return all([item in cls.alphabet() for item in str_value])

    @classmethod
    def alphabet(cls):
        raise NotImplementedError(f"{cls.__name__} is an abstract class!")


class NucleicAcid(Polymer):
    def __new__(cls, seq, annotation=None):
        return super().__new__(cls, seq)

    def __init__(self, seq, annotation=None):
        self.annotation = annotation

    def transcribe(self):
        raise NotImplementedError

    def reverse_transcribe(self):
        raise NotImplementedError

    def translate(self):
        raise NotImplementedError


class DNA(NucleicAcid):
    @classmethod
    def alphabet(cls):
        return utils.dNTPs

    def transcribe(self):
        return RNA(self.replace('T', 'U'))

    def reverse_transcribe(self):
        return self

    def translate(self, code=Code()):
        mRNA = self.transcribe()
        return mRNA.translate(code)


class RNA(NucleicAcid):
    @ classmethod
    def alphabet(cls):
        return utils.rNTPs

    def transcribe(self):
        return self

    def reverse_transcribe(self):
        return DNA(self.replace('U', 'T'))

    def translate(self, code=Code()):
        return Code(code).translate(self)


class Protein(Polymer):
    @ classmethod
    def alphabet(cls):
        return utils.aminoacids
