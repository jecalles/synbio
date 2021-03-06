from synbio.libraries import Library
from synbio.polymers import DNA, RNA
from synbio.tests import utils as testutils


class TestLibrary:
    dna_str = "ATCGAATTCCGG"
    dna_obj = DNA(dna_str)
    rna_str = "AUCGAAUUCCGG"
    rna_obj = RNA(rna_str)

    def test_init(self):
        dna_str = self.dna_str
        dna_obj = self.dna_obj
        rna_str = self.rna_str
        rna_obj = self.rna_obj
        # should pass
        Library(base_seq=dna_str)  # from raw dna string
        Library(base_seq=dna_str,
                seq_type=DNA)  # from dna string w/ proper SeqType
        Library(base_seq=dna_obj)  # from DNA object

        Library(base_seq=rna_str,
                seq_type=RNA)  # from rna string w/ proper SeqType
        Library(base_seq=rna_obj)  # from RNA object

        # should fail
        assert isinstance(
            testutils.raises(
                Library, [], {'base_seq': dna_str, 'seq_type': RNA}
            ), ValueError
        )
        assert isinstance(
            testutils.raises(
                Library, [], {'base_seq': rna_str, 'seq_type': DNA}
            ), ValueError
        )
        assert isinstance(
            testutils.raises(
                Library, [], {'base_seq': dna_str, 'seq_type': str}
            ), TypeError
        )

    def test_variance(self):
        # TODO: write this test
        pass
