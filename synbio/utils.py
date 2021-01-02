import os
import pickle

# define scope of package
__all__ = [
    # definitions
    "dNTPs", "rNTPs", "aminoacids", "triplet_codons", "triplet_mut_pairs",
    "dna_basepair_WC", "rna_basepair_WC",
    "quadruplet_codons", "quadruplet_mut_pairs", "PRS", "kdHydrophobicity",
    # functions
    "get_codons", "reverse_complement"
]

path = os.path.dirname(os.path.abspath(__file__))
with open(path + '/res/utils_definitions.pickle', 'rb') as handle:
    un_pickled = pickle.load(handle)
    [
        dNTPs, rNTPs, aminoacids, triplet_codons, triplet_mut_pairs,
        quadruplet_codons, quadruplet_mut_pairs,
        dna_basepair_WC, rna_basepair_WC,
        PRS, kdHydrophobicity
    ] = un_pickled

#############
# functions #
#############


def get_codons(seq, n=3):
    '''
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
    '''
    # check that sequence is divisible by n
    if len(seq) % n != 0:
        raise ValueError(f"seq is not divisible by n ({n})")

    num_codons = int(len(seq) / n)
    codons = [
        seq[n*i:n*(i+1)] for i in range(num_codons)
    ]

    return codons


def reverse_complement(seq, complement=dna_basepair_WC):
    return ''.join(
        [complement[nt] for nt in seq[::-1]]
    )
