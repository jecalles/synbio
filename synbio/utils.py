import itertools
from typing import Dict, List

from synbio.interfaces import SeqType, LocationType

# define scope of package
__all__ = [
    ###############
    # definitions #
    ###############
    # monomers
    "dNTPs", "rNTPs", "nonstandard_NTPs", "extended_dNTPs", "extended_rNTPs",
    "aminoacids", "nonstandard_aminoacids", "extended_aminoacids",
    # codons
    "triplet_rna_codons", "triplet_rna_mut_pairs",
    "quadruplet_rna_codons", "quadruplet_rna_mut_pairs",
    "triplet_dna_codons", "triplet_dna_mut_pairs",
    "quadruplet_dna_codons", "quadruplet_dna_mut_pairs",
    # basepairing
    "dna_basepairing", "rna_basepairing", "nonstandard_basepairing",
    "extended_dna_basepairing", "extended_rna_basepairing",
    # aminoacid physicochemistry
    "PRS", "kdHydrophobicity",
    #############
    # functions #
    #############
    "get_codons", "reverse_complement", "is_palindrome", "find_subseq",
    "all_single_mutations", "mutation_pairs"
]

###############
# definitions #
###############
# define NTPs list
dNTPs = ['T', 'C', 'A', 'G']
rNTPs = ['U', 'C', 'A', 'G']
nonstandard_NTPs = [
    'I',  # Inosine
    'R',  # A/G/I (puRine)
    'Y',  # C/T/U (pYrimidine)
    'K',  # G/T/U (have Ketone groups)
    'M',  # A/C (have aMino groups)
    'S',  # C/G (Strong interaction)
    'W',  # A/T (Weak interaction)
    'B',  # NOT A (e.g., C, G, T, or U)
    'D',  # NOT C (e.g., A, G, T, or U)
    'H',  # NOT G (e.g., A, C, T, or U)
    'V',  # NOT T/U (e.g., A, C, or G)
    'N',  # any (e.g., A, C, T, G, or U)
    '-',  # gap / no nucleic acid
]
extended_dNTPs = dNTPs + nonstandard_NTPs
extended_rNTPs = rNTPs + nonstandard_NTPs
# define list of amino acids, including stop
aminoacids = [
    'G', 'A', 'V', 'L', 'I', 'P', 'M', 'C', 'S', 'F', 'Y',
    'W', 'T', 'N', 'Q', 'D', 'E', 'R', 'H', 'K', '*'
]
nonstandard_aminoacids = [
    'B',  # Asp or Asn
    'J',  # Leu or Ile
    'O',  # Pyr / Pyrrolysine
    'U',  # Sel / Selenocysteine
    'Z',  # Glu or Gln
    'X',  # any amino acid
    '-',  # gap / no amino acid
]
extended_aminoacids = aminoacids + nonstandard_aminoacids
# define polar requirement scale
PRS = {
    'F': 5.0,
    'L': 4.9,
    'I': 4.9,
    'M': 5.3,
    'V': 5.6,
    'S': 7.5,
    'P': 6.6,
    'T': 6.6,
    'A': 7.0,
    'Y': 5.7,
    'H': 8.4,
    'Q': 8.6,
    'N': 10.0,
    'K': 10.1,
    'D': 13.0,
    'E': 12.5,
    'C': 11.5,
    'W': 5.3,
    'R': 9.1,
    'G': 7.9,
}
# define default hydrophobicity metric
kdHydrophobicity = {
    'I': 4.5,
    'V': 4.2,
    'L': 3.8,
    'F': 2.8,
    'C': 2.5,
    'M': 1.9,
    'A': 1.8,
    '*': 0,  # stop codon's given neutral value for visualization
    '0': 0,  # null codon's given neutral value for visualization
    'G': -0.4,
    'T': -0.7,
    'W': -0.9,
    'S': -0.8,
    'Y': -1.3,
    'P': -1.6,
    'H': -3.2,
    'E': -3.5,
    'Q': -3.5,
    'D': -3.5,
    'N': -3.5,
    'K': -3.9,
    'R': -4.5,
}

# define codon sets
triplet_rna_codons = [
    ''.join(tup)
    for tup in itertools.product(rNTPs, repeat=3)
]
quadruplet_rna_codons = [
    ''.join(tup)
    for tup in itertools.product(rNTPs, repeat=4)
]
triplet_dna_codons = [
    ''.join(tup)
    for tup in itertools.product(dNTPs, repeat=3)
]
quadruplet_dna_codons = [
    ''.join(tup)
    for tup in itertools.product(dNTPs, repeat=4)
]

# define Watson Crick Wobbling Rules
dna_basepairing = {
    'T': 'A',
    'C': 'G',
    'A': 'T',
    'G': 'C'
}
rna_basepairing = {
    'U': 'A',
    'C': 'G',
    'A': 'U',
    'G': 'C'
}
nonstandard_basepairing = {
    'I': 'H',
    'R': 'Y',
    'Y': 'R',
    'K': 'M',
    'M': 'K',
    'S': 'S',
    'W': 'W',
    'B': 'V',
    'V': 'B',
    'D': 'H',
    'H': 'D',
    'N': 'N',
    '-': '-',
}
extended_dna_basepairing = {**dna_basepairing, **nonstandard_basepairing}
extended_rna_basepairing = {**rna_basepairing, **nonstandard_basepairing}


# define all pairs of codons 1 mutation away
def all_single_mutations(codon, NTPs):
    return set(
        codon[:i] + nt + codon[i + 1:]
        for i, base in enumerate(codon)
        for nt in NTPs if not nt == base
    )


def mutation_pairs(codons, NTPs):
    return {
        codon: all_single_mutations(codon, NTPs)
        for codon in codons
    }


triplet_rna_mut_pairs = mutation_pairs(triplet_rna_codons, rNTPs)
quadruplet_rna_mut_pairs = mutation_pairs(quadruplet_rna_codons, rNTPs)
triplet_dna_mut_pairs = mutation_pairs(triplet_dna_codons, dNTPs)
quadruplet_dna_mut_pairs = mutation_pairs(quadruplet_dna_codons, dNTPs)


#############
# functions #
#############
def get_codons(seq: SeqType, n=3) -> List[SeqType]:
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


def reverse_complement(seq: SeqType,
                       complement: Dict[str, str] = dna_basepairing) -> str:
    return ''.join(
        complement[nt.upper()] for nt in seq[::-1]
    )


def is_palindrome(seq: SeqType) -> bool:
    comp_seq = seq.casefold()
    return comp_seq == comp_seq[::-1]


def find_subseq(seq: SeqType, subseq: SeqType) -> List[LocationType]:
    """
    >>> from synbio.utils import find_subseq
    >>>
    >>> seq_to_search = "xxXatATaTxXatatxx"
    >>> subseq_to_find = "atat"     # search is NOT case sensitive
    >>>
    >>> subseq_ix = find_subseq(seq_to_search, subseq_to_find)
    >>> subseq_ix
    [slice(3, 7, None), slice(5, 9, None), slice(11, 15, None)]
    >>> subseqs = [seq_to_search[ix] for ix in subseq_ix]
    ['atAT', 'ATaT', 'atat']
    """

    seq_len = len(seq)
    subseq_len = len(subseq)
    all_substrings = (seq[ix:ix+subseq_len] for ix in range(seq_len-subseq_len+1))
    matches = (slice(ix, ix + subseq_len) for ix, sub in enumerate(
        all_substrings) if str(sub).casefold() == str(subseq).casefold())
    return list(matches)