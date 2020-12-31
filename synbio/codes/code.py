# TODO: refactor code, now that utils has been split up
from synbio.codes import utils as codeutils
from synbio import utils


class Code(dict):
    ''' A class used to represent genetic codes. '''

    def __init__(self, code=None):
        '''Automatically loads object with a Code and a
        comparison function between amino acids "ordering". bool norm is used
        to tell dict_to_graph whether or not to set node values based on
        absolute value of metric or just the ordering. This nomenclature is
        consistent throughout the object.

        Parameters
        ----------
            dict code=None: a python dict representing a genetic code, optionally takes a string 'random' to assign a random code
        Returns
        -------
            Code obj: returns a Code object '''
        # optionally generate a random code, or load a preset code
        if type(code) == str:
            code_options = {
                'STANDARD': codeutils.definitions.standard_code,
                'COLORADO': codeutils.definitions.colorado_code,
                'RED20': codeutils.definitions.RED20,
                'RED15': codeutils.definitions.RED15,
                'RANDOM': codeutils.functions.random_code()
            }
            try:
                code = code_options[code.upper()]
            except:
                raise ValueError(
                    'Code string not recognized. Use one of the following options: {0}'.format(
                        set(code_options.keys())
                    )
                )
        # default to standard code
        elif code is None:
            code = codeutils.definitions.standard_code
        # else, try and cast input to a dict
        else:
            try:
                code = dict(code)
            except:
                raise ValueError(
                    'Code input parameter not recognized. Pass in a dictionary, Code object, or one of the following strings: {0}'.format(
                        set(code_options.keys())
                    )
                )

        # finally, call super() for init
        super().__init__(code)

        # Assign additional attributes
        self.ambiguous = codeutils.functions.is_promiscuous(code)
        self.one_to_one = codeutils.functions.is_one_to_one(code)
        self.codon_length = len(next(iter(self)))

    def __repr__(self):
        code = self.table()

        crossline = '-'*25 + '\n'

        out = object.__repr__(self) + '\n'
        for col in code:
            out += crossline
            for row in col:
                out += '|'
                for elem in row:
                    out += (elem + '|')
                out += '\n'
        out += crossline[:-1]
        return out

    def table(self):
        '''a method used to represent a genetic code as a 4x4x4 array
        '''
        rNTPs = utils.rNTPs
        out = [[[c1+c2+c3 + ':' + self[c1+c2+c3] for c2 in rNTPs]
                for c3 in rNTPs]
               for c1 in rNTPs]

        return out

    def rmap(self):
        '''
        A method used to generate the reverse map of a genetic code. Returns an amino acid -> set(codon) dictionary.

        Parameters
        ----------
            None

        Returns
        -------
            dict rmap
        '''
        rmap = {}
        for c, aa in self.items():
            codons = rmap.get(aa, [])
            codons.append(c)
            rmap.update({aa: codons})

        return rmap

    def translate(self, seq):
        '''
        A method used to translate a RNA sequence into its corresponding
        Protein sequence. Raises an error if the input sequence length is not
        divisible by the codon length, or if the sequence is not RNA.

        Parameters
        ----------
            str seq: string representing gene to translate

        Returns
        -------
            Protein prot_seq: Protein obj representing translated input sequence
        '''
        codons = utils.get_codons(seq, self.codon_length)
        return ''.join(
            [self[c] for c in codons]
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
        # handle error if code is not one-to-one
        if not self.one_to_one:
            raise TypeError(
                'cannot reverse translate sequence. genetic code is not one-to-one')
        # otherwise, create reverse translation dictionary
        rev_dict = {aa: codon for codon, aa in self.items()}
        rev_dict['*'] = stop_codon
        # translate gene and return
        return ''.join(
            [rev_dict[aa] for aa in prot_seq]
        )

    def recode(self, gene, encoding=None):
        '''
        A method used to recode an input sequence, given an initial genetic code, into an RNA sequence in this genetic code. Note: requires input sequence to be in RNA

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
        in_code = Code(encoding)
        protein = in_code.translate(gene)
        return self.reverse_translate(protein)
