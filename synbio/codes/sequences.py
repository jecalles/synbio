def transcribe(dna_seq):
    '''
    A function that converts a DNA sequence into an RNA sequence

    Parameters
    ----------
        str dna_seq: string representing DNA sequence

    Returns
    -------
        str rna_seq: string representting RNA sequence
    '''
    return dna_seq.upper().replace('T', 'U')

def reverse_transcribe(rna_seq):
    '''
    A function that converts an RNA sequence into a DNA sequence

    Parameters
    ----------
        str rna_seq: string representting RNA sequence

    Returns
    -------
        str dna_seq: string representing DNA sequence
    '''
    return rna_seq.upper().replace('U', 'T')

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
