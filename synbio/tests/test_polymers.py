from synbio.polymers import *


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
        seq = DNA("AATTCCGG")
        seq.__delitem__(slice(2, 4))
        assert seq == "AACCGG"
        seq.__delitem__(Location(2, 4))
        assert seq == "AAGG"
        seq.__delitem__(0)
        assert seq == "AGG"

    def test_insert(self):
        seq = DNA("ATCG")
        seq.insert(2, "T")

    def test_part_integration(self):
        dna = DNA("ATCGAATTCCGG")
        part1 = Part(seq=dna, location=Location(0, 4))
        part2 = Part(seq=dna, location=Location(4, 8))
        part3 = Part(seq=dna, location=Location(8, 12))
        part4 = Part(seq=dna, location=Location(2, 10))

        # test indexing
        assert part1.seq == "ATCG"
        assert part2.seq == "AATT"
        assert part3.seq == "CCGG"
        assert part4.seq == "CGAATTCC"

        # test reassignment from DNA obj
        dna[4:8] = "AT"
        assert part1.seq == "ATCG"
        assert part1.location == Location(0, 4)

        assert part2.seq == "AT"
        assert part2.location == Location(4, 6)

        assert part3.seq == "CCGG"
        assert part3.location == Location(6, 10)

        assert part4.seq == "CGATCC"
        assert part4.location == Location(2, 8)

        dna[Location(4, 6)] = "AATT"
        assert part1.seq == "ATCG"
        assert part1.location == Location(0, 4)

        assert part2.seq == "AATT"
        assert part2.location == Location(4, 8)

        assert part3.seq == "CCGG"
        assert part3.location == Location(8, 12)

        assert part4.seq == "CGAATTCC"
        assert part4.location == Location(2, 10)

    def test_duplicate_parts(self):
        # TODO: write test that ensures that only unique parts are found in DNA.annotations
        x = DNA("AAAA")
        x_annotation = Part(seq=x, name='x', kind='source')
        y = DNA("TTTT")
        y_annotation = Part(seq=y, name='y', kind='source')

        assert len(x.annotations) == 1
        assert len(y.annotations) == 1
        assert x.annotations != y.annotations

    def test_central_dogma(self):
        gfp_str = "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAG"

        GFP = DNA(gfp_str)
        GFP_transcript = GFP.transcribe()
        GFP_prot1 = GFP.translate()
        GFP_prot2 = GFP_transcript.translate()
        assert GFP_prot1 == GFP_prot2


if __name__ == '__main__':
    TestDNA().test_duplicate_parts()
