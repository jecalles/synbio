import pickle
from pathlib import Path
from functools import wraps

# define scope of package
__all__ = [
    # definitions
    "dNTPs", "rNTPs", "aminoacids", "triplet_codons", "triplet_mut_pairs",
    "dna_basepair_WC", "rna_basepair_WC",
    "quadruplet_codons", "quadruplet_mut_pairs", "PRS", "kdHydrophobicity",
    # functions
    "get_codons", "reverse_complement",
    # wrappers
    "compare_same_type",
    # Mixins
    "ComparableMixin"
]

#############
# definitions#
#############
module_path = Path(__file__).absolute().parent
definitions_path = Path(module_path, "./res/utils_definitions.pickle").absolute()
with open(definitions_path, 'rb') as handle:
    un_pickled = pickle.load(handle)
    [
        dNTPs, rNTPs, aminoacids, triplet_codons, triplet_mut_pairs,
        quadruplet_codons, quadruplet_mut_pairs,
        dna_basepair_WC, rna_basepair_WC,
        PRS, kdHydrophobicity
    ] = un_pickled


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
# functions #
#############
def get_codons(seq, n=3):
    """
    A function that takes an input DNA/RNA sequence representing an open
    reading frame (ORF)and splits that sequence into codons of length n
    (default=3)

    Parameters
    ----------
        str seq: string representing an ORF (DNA or RNA)
        int n: length of codons (default=3)

    Returns
    -------
        list<str> codons: input sequence split into codons
    """
    # check that sequence is divisible by n
    if len(seq) % n != 0:
        raise ValueError(f"seq is not divisible by n ({n})")

    num_codons = int(len(seq) / n)
    codons = [
        seq[n * i:n * (i + 1)] for i in range(num_codons)
    ]

    return codons


def reverse_complement(seq, complement=dna_basepair_WC):
    return ''.join(
        complement[nt] for nt in seq[::-1]
    )


#############
#  Mixins   #
#############
class ComparableMixin:
    """
    A Mixin class that automatically provides some comparison dunder methods,
    given a class attribute _comparables
    """
    _comparables = list()

    def __hash__(self):
        return hash(
            (getattr(self, comp) for comp in self._comparables)
        )

    @compare_same_type
    def __eq__(self, other):
        return all(
            getattr(self, comp) == getattr(other, comp)
            for comp in self._comparables
        )
