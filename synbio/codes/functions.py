from typing import Dict, Set

from synbio.codes import Code


def get_synonymous_codons(codon: str, code: Code) -> Set[str]:
    return code.rmap()[code[codon]]


def get_nonsynonymous_codons(
        codon: str,
        code: Code,
        codon_frequency: Dict[str, float]) -> Set[str]:
    """
    this code doesn't do what it thinks it does! whaat is this supposed to do?
    """
    raise NotImplementedError("This function under construction")

    def max_codon(codon_list, freq):
        # finds the codon in the list with the highest frequency
        codon_and_frequency_pairs = (
            (c, freq[c])
            for c in codon_list
        )
        return max(
            codon_and_frequency_pairs,
            key=lambda tup: tup[1]
        )[0]

    return [
        max_codon(codon_list, codon_frequency)
        for aa, codon_list in code.rmap().items()
        if aa != code[codon]
    ]
