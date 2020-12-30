class Polymer(str):
    def __new__(cls, value):
        if not cls._alphabet_check(value):
            raise ValueError(f"input value not in {cls.__name__} alphabet")
        else:
            return super().__new__(cls, value.upper())

    @classmethod
    def _alphabet_check(cls, value):
        str_value = str(value).upper()
        return all([item in cls.alphabet() for item in str_value])

    @classmethod
    @abstractmethod
    def alphabet(cls):
        raise NotImplementedError(f"{cls.__name__} is an abstract class!")

class DNA(Polymer):
    @classmethod
    def alphabet(cls):
        # TODO: implement using utils.dNTPs, but first you gotta refactor, oof
        pass

class RNA(Polymer):
    @classmethod
    def alphabet(cls):
        # TODO: implement using utils.rNTPs, but first you gotta refactor, oof
        pass

class Protein(Polymer):
    @classmethod
    def alphabet(cls):
        # TODO: implement, and also update utils to have amino acid library
        pass
