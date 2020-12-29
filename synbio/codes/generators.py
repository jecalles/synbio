from random import choice

from .code import Code
from .utils import definitions

def red20():
    '''
    A generator that randomly generates a stream of RED20 codes as python dictionaries without yielding repeated codes

    Parameters
    ----------
        None

    Returns
    -------
        dict r20: a python dict representation of a genetic code
    '''
    rmap = Code().rmap()
    cache = set()

    while True:
        code = {codon:"0" for codon in definitions.triplet_codons}
        choices = [(choice(c_set), aa) for aa, c_set in rmap.items()]
        code.update(choices)
        hashed = hash(str(code))

        if hashed not in cache:
            cache.add(hashed)
            yield code
