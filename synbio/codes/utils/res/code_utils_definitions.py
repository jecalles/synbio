import os
import pickle
from synbio.utils import rNTPs

# define block structures
unrestricted_block = {}
count = 0
for c1 in rNTPs:
    for c2 in rNTPs:
        for c3 in rNTPs:
            unrestricted_block[count] = [c1+c2+c3]
            count += 1

standard_block = {}
for i in range(48):
    standard_block[i] = []

count = 0
for c1 in rNTPs:
    for c2 in rNTPs:
        for c3 in rNTPs:
            standard_block[count].append(c1+c2+c3)
            if (c3 != 'U'):
                count += 1

natural_block = {
    0: ['UUU', 'UUC'],
    1: ['UUA', 'UUG'],
    2: ['CUU', 'CUC', 'CUA', 'CUG'],
    3: ['AUU', 'AUC', 'AUA'],
    4: ['AUG'],
    5: ['GUU', 'GUC', 'GUA', 'GUG'],
    6: ['UCU', 'UCC', 'UCA', 'UCG'],
    7: ['CCU', 'CCC', 'CCA', 'CCG'],
    8: ['ACU', 'ACC', 'ACA', 'ACG'],
    9: ['GCU', 'GCC', 'GCA', 'GCG'],
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
    'U': 'A',
    'C': 'G',
    'A': 'U',
    'G': 'C'
}
wobble_WC = {
    'U': ['A', 'G'],
    'C': ['G'],
    'A': ['U', 'C'],
    'G': ['U', 'C'],
    'I': ['A', 'C', 'U']
}

# define standard code
standard_code = {
    'UUU': 'F',
    'UUC': 'F',
    'UUA': 'L',
    'UUG': 'L',
    'UCU': 'S',
    'UCC': 'S',
    'UCA': 'S',
    'UCG': 'S',
    'UAU': 'Y',
    'UAC': 'Y',
    'UAA': '*',
    'UAG': '*',
    'UGU': 'C',
    'UGC': 'C',
    'UGA': '*',
    'UGG': 'W',
    'CUU': 'L',
    'CUC': 'L',
    'CUA': 'L',
    'CUG': 'L',
    'CCU': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',
    'CAU': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGU': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AUU': 'I',
    'AUC': 'I',
    'AUA': 'I',
    'AUG': 'M',
    'ACU': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',
    'AAU': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'AGU': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GUU': 'V',
    'GUC': 'V',
    'GUA': 'V',
    'GUG': 'V',
    'GCU': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',
    'GAU': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'GGU': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',
}
# define refactored [sic] code from Pines et al 2017 (aka Colorado code)
colorado_code = {
    'GAA': 'V',
    'UCG': 'V',
    'CGU': 'V',
    'UGA': 'L',
    'AAU': 'L',
    'CUC': 'L',
    'CCA': 'I',
    'GGG': 'I',
    'UUU': 'I',
    'UAC': 'I',
    'CAG': 'A',
    'AUA': 'A',
    'GCU': 'A',
    'AGC': 'A',
    'GAU': 'E',
    'ACA': 'E',
    'UUC': 'E',
    'CGG': 'E',
    'UGU': 'D',
    'AAC': 'D',
    'GUG': 'D',
    'UAA': '*',
    'UCU': 'P',
    'AUG': 'P',
    'GUC': 'P',
    'CAA': 'P',
    'GAC': 'T',
    'UCA': 'T',
    'CCC': 'S',
    'AGG': 'S',
    'AUU': 'Q',
    'GGA': 'Q',
    'UGC': 'N',
    'CAU': 'N',
    'GCG': 'M',
    'CUA': 'M',
    'AAA': 'C',
    'UUG': 'C',
    'GGU': 'C',
    'CUU': 'G',
    'AGU': 'G',
    'ACC': 'G',
    'UAG': 'G',
    'UGG': 'R',
    'GCA': 'R',
    'CAC': 'R',
    'GGC': 'H',
    'CCG': 'H',
    'UUA': 'H',
    'ACU': 'H',
    'CGA': 'K',
    'UCC': 'K',
    'GUU': 'K',
    'AAG': 'K',
    'CCU': 'Y',
    'GAG': 'Y',
    'AUC': 'Y',
    'CGC': 'W',
    'ACG': 'W',
    'GUA': 'W',
    'UAU': 'W',
    'GCC': 'F',
    'CUG': 'F',
    'AGA': 'F',
}

if __name__ == '__main__':
    # time to pickle!
    toDump = [
        unrestricted_block, standard_block, natural_block,
        basepair_WC, wobble_WC,
        standard_code, colorado_code
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
