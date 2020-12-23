import pickle

# define NTPs list
dNTPs = ['T', 'C', 'A', 'G']
rNTPs = ['U', 'C', 'A', 'G']
# define list of natural amino acids, including stop
residues = [
    'G', 'A', 'V', 'L', 'I', 'P', 'M', 'C', 'S', 'F', 'Y',
    'W', 'T', 'N', 'Q', 'D', 'E', 'R', 'H', 'K', '*'
];
# define polar requirement scale
PRS = {
    'F' : 5.0,
    'L' : 4.9,
    'I' : 4.9,
    'M' : 5.3,
    'V' : 5.6,
    'S' : 7.5,
    'P' : 6.6,
    'T' : 6.6,
    'A' : 7.0,
    'Y' : 5.7,
    'H' : 8.4,
    'Q' : 8.6,
    'N' : 10.0,
    'K' : 10.1,
    'D' : 13.0,
    'E' : 12.5,
    'C' : 11.5,
    'W' : 5.3,
    'R' : 9.1,
    'G' : 7.9,
}

# define default hydrophobicity metric
kdHydrophobicity = {
    'I' : 4.5,
    'V' : 4.2,
    'L' : 3.8,
    'F' : 2.8,
    'C' : 2.5,
    'M' : 1.9,
    'A' : 1.8,
    '*' : 0, #stop codon's given neutral value for visualization
    '0' : 0, #null codon's given neutral value for visualization
    'G' : -0.4,
    'T' : -0.7,
    'W' : -0.9,
    'S' : -0.8,
    'Y' : -1.3,
    'P' : -1.6,
    'H' : -3.2,
    'E' : -3.5,
    'Q' : -3.5,
    'D' : -3.5,
    'N' : -3.5,
    'K' : -3.9,
    'R' : -4.5,
}

# define block structures and codon set
unrestricted_block = {}
triplet_codons = []
count = 0;
for c1 in rNTPs:
    for c2 in rNTPs:
        for c3 in rNTPs:
            triplet_codons.append(c1+c2+c3)
            unrestricted_block[count] = [c1+c2+c3]
            count+=1

standard_block = {}
for i in range(48):
    standard_block[i] = []

count = 0
for c1 in rNTPs:
    for c2 in rNTPs:
        for c3 in rNTPs:
            standard_block[count].append(c1+c2+c3)
            if (c3 != 'U'):
                count+=1

# define quadruplet codon set
quadruplet_codons = []
for nt1 in rNTPs:
    for nt2 in rNTPs:
        for nt3 in rNTPs:
            for nt4 in rNTPs:
                quadruplet_codons.append(nt1+nt2+nt3+nt4)

natural_block = {
    0 : ['UUU', 'UUC'],
    1 : ['UUA', 'UUG'],
    2 : ['CUU', 'CUC', 'CUA', 'CUG'],
    3 : ['AUU', 'AUC', 'AUA'],
    4 : ['AUG'],
    5 : ['GUU', 'GUC', 'GUA', 'GUG'],
    6 : ['UCU', 'UCC', 'UCA', 'UCG'],
    7 : ['CCU', 'CCC', 'CCA', 'CCG'],
    8 : ['ACU', 'ACC', 'ACA', 'ACG'],
    9 : ['GCU', 'GCC', 'GCA', 'GCG'],
    10: ['UAU', 'UAC'],
    11: ['UAA', 'UAG'],
    12: ['CAU', 'CAC'],
    13: ['CAA', 'CAG'],
    14: ['AAU', 'AAC'],
    15: ['AAA', 'AAG'],
    16: ['GAU', 'GAC'],
    17: ['GAA', 'GAG'],
    18: ['UGU', 'UGC'],
    19: ['UGA'],
    20: ['UGG'],
    21: ['CGU', 'CGC', 'CGA', 'CGG'],
    22: ['AGU', 'AGC'],
    23: ['AGA', 'AGG'],
    24: ['GGU', 'GGC', 'GGA', 'GGG']
}

# define Watson Crick Wobbling Rules
basepair_WC = {
    'U' : 'A',
    'C' : 'G',
    'A' : 'U',
    'G' : 'C'
}
wobble_WC = {
    'U' : ['A', 'G'],
    'C' : ['G'],
    'A' : ['U', 'C'],
    'G' : ['U', 'C'],
    'I' : ['A', 'C', 'U']
}

#define all pairs of codons 1 mutation away
triplet_mut_pairs = set()
for codon in triplet_codons:
    for i, base in enumerate(codon):
        for nt in rNTPs:
            # handle if nt is the same as base
            if nt == base:
                continue
            # if not, generate new codon
            c_new = codon[:i] + nt + codon[i+1:]
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
            c_new = codon[:i] + nt + codon[i+1:]
            # add to set
            quadruplet_mut_pairs.add((codon, c_new))

# define standard code
standard_code = {
    'UUU':'F',
    'UUC':'F',
    'UUA':'L',
    'UUG':'L',
    'UCU':'S',
    'UCC':'S',
    'UCA':'S',
    'UCG':'S',
    'UAU':'Y',
    'UAC':'Y',
    'UAA':'*',
    'UAG':'*',
    'UGU':'C',
    'UGC':'C',
    'UGA':'*',
    'UGG':'W',
    'CUU':'L',
    'CUC':'L',
    'CUA':'L',
    'CUG':'L',
    'CCU':'P',
    'CCC':'P',
    'CCA':'P',
    'CCG':'P',
    'CAU':'H',
    'CAC':'H',
    'CAA':'Q',
    'CAG':'Q',
    'CGU':'R',
    'CGC':'R',
    'CGA':'R',
    'CGG':'R',
    'AUU':'I',
    'AUC':'I',
    'AUA':'I',
    'AUG':'M',
    'ACU':'T',
    'ACC':'T',
    'ACA':'T',
    'ACG':'T',
    'AAU':'N',
    'AAC':'N',
    'AAA':'K',
    'AAG':'K',
    'AGU':'S',
    'AGC':'S',
    'AGA':'R',
    'AGG':'R',
    'GUU':'V',
    'GUC':'V',
    'GUA':'V',
    'GUG':'V',
    'GCU':'A',
    'GCC':'A',
    'GCA':'A',
    'GCG':'A',
    'GAU':'D',
    'GAC':'D',
    'GAA':'E',
    'GAG':'E',
    'GGU':'G',
    'GGC':'G',
    'GGA':'G',
    'GGG':'G',
}
# define refactored [sic] code from Pines et al 2017 (aka Colorado code)
colorado_code = {
    'GAA' : 'V',
    'UCG' : 'V',
    'CGU' : 'V',
    'UGA' : 'L',
    'AAU' : 'L',
    'CUC' : 'L',
    'CCA' : 'I',
    'GGG' : 'I',
    'UUU' : 'I',
    'UAC' : 'I',
    'CAG' : 'A',
    'AUA' : 'A',
    'GCU' : 'A',
    'AGC' : 'A',
    'GAU' : 'E',
    'ACA' : 'E',
    'UUC' : 'E',
    'CGG' : 'E',
    'UGU' : 'D',
    'AAC' : 'D',
    'GUG' : 'D',
    'UAA' : '*',
    'UCU' : 'P',
    'AUG' : 'P',
    'GUC' : 'P',
    'CAA' : 'P',
    'GAC' : 'T',
    'UCA' : 'T',
    'CCC' : 'S',
    'AGG' : 'S',
    'AUU' : 'Q',
    'GGA' : 'Q',
    'UGC' : 'N',
    'CAU' : 'N',
    'GCG' : 'M',
    'CUA' : 'M',
    'AAA' : 'C',
    'UUG' : 'C',
    'GGU' : 'C',
    'CUU' : 'G',
    'AGU' : 'G',
    'ACC' : 'G',
    'UAG' : 'G',
    'UGG' : 'R',
    'GCA' : 'R',
    'CAC' : 'R',
    'GGC' : 'H',
    'CCG' : 'H',
    'UUA' : 'H',
    'ACU' : 'H',
    'CGA' : 'K',
    'UCC' : 'K',
    'GUU' : 'K',
    'AAG' : 'K',
    'CCU' : 'Y',
    'GAG' : 'Y',
    'AUC' : 'Y',
    'CGC' : 'W',
    'ACG' : 'W',
    'GUA' : 'W',
    'UAU' : 'W',
    'GCC' : 'F',
    'CUG' : 'F',
    'AGA' : 'F',
}

if __name__ == '__main__':
    # time to pickle!
    toDump = [dNTPs, rNTPs, residues, triplet_codons, triplet_mut_pairs,
                quadruplet_codons, quadruplet_mut_pairs,
                PRS, kdHydrophobicity, 
                unrestricted_block, standard_block, natural_block,
                basepair_WC, wobble_WC,
                standard_code, colorado_code]
    with open('./utils_definitions.pickle', 'wb') as handle:
        pickle.dump(toDump, handle)

    # test the pickle
    with open('./utils_definitions.pickle', 'rb') as handle:
        unDumped = pickle.load(handle)
    # taste the pickle
    if (toDump == unDumped) :
        print("Mmmm, that's a tasty pickle")
    else:
        print("Hmmm, something's up with this pickle")
