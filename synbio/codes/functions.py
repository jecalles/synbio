
# TODO: write tests and docstrings for the functions below


def get_synonymous_codons(codon, code):
    return code.rmap()[code[codon]]


def get_nonsynonymous_codons(codon, code, codon_frequency):
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
