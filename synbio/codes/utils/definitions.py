# import necessary modules
import os
from math import comb as binomial
import pickle
import random
from collections import deque
from copy import copy

# define scope of package
__all__ = [
    "dNTPs", "rNTPs", "residues", "triplet_codons", "triplet_mut_pairs",
    "quadruplet_codons", "quadruplet_mut_pairs", "PRS", "kdHydrophobicity",
    "unrestricted_block", "standard_block", "natural_block", "basepair_WC",
    "wobble_WC", "standard_code", "colorado_code", "RED20", "RED15"
]
def __dir__():
    default = [key for key in globals().keys() if key[:2] == '__']
    return default + __all__

#####################
# define properties #
#####################
path = os.path.dirname(os.path.abspath(__file__))
with open(path + '/res/utils_definitions.pickle', 'rb') as handle:
    un_pickled = pickle.load(handle)
    [
        dNTPs, rNTPs, residues, triplet_codons, triplet_mut_pairs,
        quadruplet_codons, quadruplet_mut_pairs,
        PRS, kdHydrophobicity,
        unrestricted_block, standard_block, natural_block,
        basepair_WC, wobble_WC,
        standard_code, colorado_code
    ] = un_pickled

with open(path + '/res/RED20.pickle', 'rb') as handle:
    RED20 = pickle.load(handle)
with open(path + '/res/RED15.pickle', 'rb') as handle:
    RED15 = pickle.load(handle)
