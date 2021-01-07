import os
import pickle

# define NTPs list
dNTPs = ['T', 'C', 'A', 'G']
rNTPs = ['U', 'C', 'A', 'G']
# define list of natural amino acids, including stop
aminoacids = [
    'G', 'A', 'V', 'L', 'I', 'P', 'M', 'C', 'S', 'F', 'Y',
    'W', 'T', 'N', 'Q', 'D', 'E', 'R', 'H', 'K', '*'
]
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
triplet_codons = []
quadruplet_codons = []
for nt1 in rNTPs:
    for nt2 in rNTPs:
        for nt3 in rNTPs:
            triplet_codons.append(nt1 + nt2 + nt3)
            for nt4 in rNTPs:
                quadruplet_codons.append(nt1 + nt2 + nt3 + nt4)

# define Watson Crick Wobbling Rules
dna_basepair_WC = {
    'T': 'A',
    'C': 'G',
    'A': 'T',
    'G': 'C'
}
rna_basepair_WC = {
    'U': 'A',
    'C': 'G',
    'A': 'U',
    'G': 'C'
}

# define all pairs of codons 1 mutation away
triplet_mut_pairs = set()
for codon in triplet_codons:
    for i, base in enumerate(codon):
        for nt in rNTPs:
            # handle if nt is the same as base
            if nt == base:
                continue
            # if not, generate new codon
            c_new = codon[:i] + nt + codon[i + 1:]
            # add to set
            triplet_mut_pairs.add((codon, c_new))

# define all pairs of quadruplet codons 1 mutation away
quadruplet_mut_pairs = set()
for codon in quadruplet_codons:
    for i, base in enumerate(codon):
        for nt in rNTPs:
            # handle if nt is the same as base
            if nt == base:
                continue
            # if not, generate new codon
            c_new = codon[:i] + nt + codon[i + 1:]
            # add to set
            quadruplet_mut_pairs.add((codon, c_new))

if __name__ == '__main__':
    # time to pickle!
    toDump = [
        dNTPs, rNTPs, aminoacids, triplet_codons, triplet_mut_pairs,
        quadruplet_codons, quadruplet_mut_pairs,
        dna_basepair_WC, rna_basepair_WC,
        PRS, kdHydrophobicity
    ]
    path = os.path.dirname(os.path.abspath(__file__))
    with open(path + '/utils_definitions.pickle', 'wb') as handle:
        pickle.dump(toDump, handle)
    # test the pickle
    with open(path + '/utils_definitions.pickle', 'rb') as handle:
        unDumped = pickle.load(handle)
    # taste the pickle
    if (toDump == unDumped):
        print("Mmmm, that's a tasty pickle")
    else:
        print("Hmmm, something's up with this pickle")
