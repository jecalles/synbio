import os
import pickle
import itertools

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
#  classes  #
#############


class Location:
    # TODO: write docstring
    '''
    '''

    def __init__(self, start, end, strand='FWD'):
        # check that start and stop are valid
        if not start <= end:
            raise ValueError(f"start position ({start}) comes after end position ({end})")

        # check that strand input is valid
        if strand not in {'FWD', 'REV'}:
            raise ValueError(f"input strand ({strand}) not FWD or REV")

        self.start = start
        self.end = end
        self.strand = strand

    def __repr__(self):
        return f"{self.__class__.__name__}({self.start}, {self.end}, {self.strand})"

    def __eq__(self, other):
        if not isinstance(other, Location):
            raise TypeError(f"cannot compare Location to {type(other)}")

        return self.start == other.start and self.end == other.end and self.strand == other.strand

    @staticmethod
    def contains(outer_loc, inner_loc):
        '''
        A static method that compares two locations, determining if the inner_loc is completely contained by the outer location.
        '''
        return (outer_loc.start <= inner_loc.start) and (inner_loc.end <= outer_loc.end)

    @staticmethod
    def overlaps(loc1, loc2):
        '''
        # TODO: write docstring for this method
        '''
        # TODO: implement lol
        return (loc1.start <= loc2.start and loc1.end > loc2.start) or (loc2.start <= loc1.start and loc2.end > loc1.start)

    @ staticmethod
    def find_overlaps(locations):
        '''
        A static method that, given a list of Location objects, constructs a graph representing all Location objects that overlap with each other.

        E.g., suppose our Locations are the following:


        a: ----------
        b:     ---
        c: ----
        d:           ----------
        e:      ----------

        >>> a = Location(0, 10)
        >>> b = Location(4, 7)
        >>> c = Location(0, 4)
        >>> d = Location(10, 20)
        >>> e = Location(5, 15)

        >>> Location.find_overlaps([a, b, c, d, e])
        [
            [1, 1, 1, 0, 1],
            [1, 1, 0, 0, 1],
            [1, 0, 1, 0, 0],
            [0, 0, 0, 1, 1],
            [1, 1, 0, 1, 1],
        ]
        '''
        # calculate overlaps
        loc_pairs = itertools.product(locations, repeat=2)
        bool_array = [
            Location.overlaps(loc1, loc2) for
            (loc1, loc2) in loc_pairs
        ]
        # reshape array to square matrix
        n = len(locations)
        bool_matrix = [
            bool_array[i*n: (i+1)*n] for i in range(n)
        ]
        return bool_matrix

    def to_slice(self):
        '''
        '''
        # TODO: write docstring
        return slice(self.start, self.end, 1)


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
