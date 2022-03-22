from random import choice
from typing import Iterable

from synbio.codes import utils as codeutils
from .code import Code


def red20() -> Iterable[Code]:
    """
    A generator that randomly generates a stream of RED20 codes as python
    dictionaries without yielding repeated codes
    """
    rmap = Code().rmap
    cache = set()

    while True:
        code = {codon: "0" for codon in codeutils.triplet_rna_codons}
        choices = [(choice(c_set), aa) for aa, c_set in rmap.items()]
        code.update(choices)
        hashed = hash(str(code))

        if hashed not in cache:
            cache.add(hashed)
            yield code
