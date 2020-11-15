from .utils import definitions
from .utils import functions

class Code:
    ''' A class used to represent genetic codes. '''

    def __init__(self, table=None):
        '''Automatically loads object with a Code and a
        comparison function between amino acids "ordering". bool norm is used
        to tell dict_to_graph whether or not to set node values based on
        absolute value of metric or just the ordering. This nomenclature is
        consistent throughout the object.

        Parameters
        ----------
        - dict table=None: a python dict representing a genetic code, optionally
            takes a string 'random' to assign a random table
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
        self.dict = table
        self.ambiguous = ambiguous
        self.one_to_one = one_to_one

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
        out = [ [ [c1+c2+c3 + ':' + self.dict[c1+c2+c3] for c2 in rNTPs]
                        for c3 in rNTPs]
                            for c1 in rNTPs]

        return out

    def reverse_translate(self, prot_seq, stop_codon='uga'):
        '''a method used to translate a given protein sequence, assuming the
        genetic code is one-to-one. raises an error otherwise

        parameters
        ----------
        - str prot_seq: string representing amino acid sequence to translate
        - str stop_codon: stop codon to use for goi (default to uga)

        returns
        -------
        str gene: string representing translated amino acid sequence
        '''
        # handle error if table is not one-to-one
        if not self.one_to_one:
            raise TypeError(
                'cannot translate sequence. genetic code is not one-to-one')
        # otherwise, create reverse translation dictionary
        rev_dict = {aa: codon for codon, aa in self.dict.items()}
        rev_dict['*'] = stop_codon
        # translate gene and return
        gene = ''
        for aa in prot_seq:
            gene += rev_dict[aa]
        return gene
