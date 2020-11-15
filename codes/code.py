from collections import UserDict

from .utils import definitions
from .utils import functions

from . import sequences

class Code(UserDict):
    ''' A class used to represent genetic codes. '''

    def __init__(self, table=None):
        '''Automatically loads object with a Code and a
        comparison function between amino acids "ordering". bool norm is used
        to tell dict_to_graph whether or not to set node values based on
        absolute value of metric or just the ordering. This nomenclature is
        consistent throughout the object.

        Parameters
        ----------
            dict table=None: a python dict representing a genetic code, optionally takes a string 'random' to assign a random table
        Returns
        -------
            Code obj: returns a Code object '''
        # optionally generate a random table, or load a preset table
        if type(table) == str:
            table_options = {
                'STANDARD': definitions.standard_code,
                'COLORADO': definitions.colorado_code,
                'RED20': definitions.RED20,
                'RED15': definitions.RED15,
                'RANDOM': functions.random_code()
            }
            try:
                table = table_options[table.upper()]
            except:
                raise ValueError(
                    'Table string not recognized. Use one of the following options: {0}'.format(
                        set(table_options.keys())
                    )
                )
        # default to standard table if not specified or wrong datatype
        elif table == None or not type(table) == dict:
            table = definitions.standard_code

        # determine table ambiguity
        ambiguous = functions.is_promiscuous(table)

        # determine if table is one-to=one
        one_to_one = functions.is_one_to_one(table)

        # Assign assign instance attributes
        self.data = table
        self.ambiguous = ambiguous
        self.one_to_one = one_to_one
        self.codon_length = len(table['AUG'])

    def __repr__(self):
        table = self.table()

        crossline = '-'*25 + '\n'

        out = object.__repr__(self) + '\n'
        for col in table:
            out += crossline
            for row in col:
                out += '|'
                for elem in row:
                    out += (elem + '|')
                out +='\n'
        out += crossline[:-1]
        return out

    def table(self):
        '''a method used to represent a genetic code as a 4x4x4 array
        '''
        rNTPs = definitions.rNTPs
        out = [ [ [c1+c2+c3 + ':' + self.data[c1+c2+c3] for c2 in rNTPs]
                        for c3 in rNTPs]
                            for c1 in rNTPs]

        return out

    def translate(self, gene):
        '''
        A method used to translate a DNA/RNA sequence into its corresponding
        protein sequence. Raises an error if the input sequence length is not
        divisible by the codon length.

        Parameters
        ----------
            str gene: string representing translated amino acid sequence

        Returns
        -------
            str prot_seq: string representing translated amino acid sequence
        '''
        mRNA = sequences.transcribe(gene)
        codons = sequences.get_codons(mRNA)
        return ''.join(
            [self.data[c] for c in codons]
        )

    def reverse_translate(self, prot_seq, stop_codon='UGA'):
        '''
        A method used to reverse translate a given protein sequence, assuming
        the genetic code is one-to-one. Raises an error otherwise

        Parameters
        ----------
            str prot_seq: string representing amino acid sequence

        Returns
        -------
            str gene: string representing reverse translated gene sequence
        '''
        # handle error if table is not one-to-one
        if not self.one_to_one:
            raise TypeError(
                'cannot reverse translate sequence. genetic code is not one-to-one')
        # otherwise, create reverse translation dictionary
        rev_dict = {aa: codon for codon, aa in self.data.items()}
        rev_dict['*'] = stop_codon
        # translate gene and return
        gene = ''
        for aa in prot_seq:
            gene += rev_dict[aa]
        return gene

    def recode(self, gene, encoding=None):
        '''
        A method used to recode an input sequence, given an initial genetic code, into an RNA sequence in this genetic code.

        Parameters
        ----------
            str gene: input gene sequence (DNA or RNA)
            Code encoding: genetic code used to encode input gene (default:
                Standard Code). Note: input must be valid input to Code()
                constructor

        Returns
        -------
            str out: output gene sequence (RNA) in this genetic code
        '''
        in_code = type(self)(encoding)
        mRNA = sequences.transcribe(gene)
        protein = in_code.translate(mRNA)
        return self.reverse_translate(protein)
