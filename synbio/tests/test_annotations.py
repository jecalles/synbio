from synbio.utils import Location
from synbio.annotations import *
from synbio.polymers import DNA


class TestPart:
    def test_DNA_integration(self):
        dna = DNA("ATCGAATTCCGG")
        part1 = Part(seq=dna, location=Location(0, 4))
        part2 = Part(seq=dna, location=Location(4, 8))
        part3 = Part(seq=dna, location=Location(8, 12))

        # test indexing
        assert part1.seq == "ATCG"
        assert part2.seq == "AATT"
        assert part3.seq == "CCGG"

        # test reassignment from DNA obj
        dna[4:8] = "AT"
        assert part1.seq == "ATCG"
        assert part1.location == Location(0, 4)

        assert part2.seq == "AT"
        assert part2.location == Location(4, 6)

        assert part3.seq == "CCGG"
        assert part3.location == Location(6, 10)

        # test reassignment from Part
        part2.seq[4:6] = "AATT"
        assert part1.seq == "ATCG"
        assert part1.location == Location(0, 4)

        assert part2.seq == "AATT"
        assert part2.location == Location(4, 8)

        assert part3.seq == "CCGG"
        assert part3.location == Location(8, 12)

    def test_initialization(self):
        assert 1 == 2

    def test_concatenation(self):
        assert 1 == 2

    def test_subparts(self):
        assert 1 == 2

    def test_parent(self):
        assert 1 == 2
