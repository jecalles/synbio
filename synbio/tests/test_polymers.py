from synbio.polymers import *
from synbio.annotations import Location


class TestDNA:
    seq = DNA("AATTCCGG")
    loc_fwd = Location(2, 7, "FWD")  # 5' --TTCCG- 3'
    loc_rev = Location(1, 5, "REV")  # 3' -TAAG--- 5'

    def test_slicing(self):
        seq = self.seq
        loc_fwd = self.loc_fwd
        loc_rev = self.loc_rev

        # sense strand slicing
        fwd_expect = "TTCCG"
        fwd_slice = seq[2:7]
        assert fwd_slice == fwd_expect
        fwd_slice = seq[loc_fwd]
        assert fwd_slice == fwd_expect

        # reverse strand slicing
        rev_expect = "GAAT"
        rev_slice = seq[loc_rev]
        assert rev_slice == rev_expect

    def test_reassignment(self):
        loc_fwd = self.loc_fwd
        loc_rev = self.loc_rev

        # sense strand reassignment
        fwd_seq = DNA("AAA")
        fwd_seq[1] = "T"
        fwd_expect = DNA("ATA")
        assert fwd_seq == fwd_expect

        fwd_loc_seq = DNA("AAAAAAAA")
        fwd_loc_seq[loc_fwd] = "T"
        fwd_loc_expect = DNA("AATA")
        assert fwd_loc_seq == fwd_loc_expect

        # reverse strand reassignment
        rev_seq = DNA("AAAAAAAA")
        rev_seq[loc_rev] = "CCCC"
        rev_expect = DNA("AGGGGAAA")
        assert rev_seq == rev_expect

    def test_delete(self):
        # TODO: test DNA.__delitem__
        assert 1 == 2

    def test_insert(self):
        # TODO: test DNA.insert()
        assert 1 == 2

    def test_central_dogma(self):
        gfp_str = "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAG"

        GFP = DNA(gfp_str)
        GFP_transcript = GFP.transcribe()
        GFP_prot1 = GFP.translate()
        GFP_prot2 = GFP_transcript.translate()
        assert GFP_prot1 == GFP_prot2
