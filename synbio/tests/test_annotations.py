from synbio.annotations import *
from synbio.polymers import DNA


class TestLocation:
    ##################
    # Test locations #
    ##################
    # a: ----------
    # b:     ----
    # c: ----
    # d:           ----------
    # e:      ----------
    a = Location(0, 10)
    b = Location(4, 7)
    c = Location(0, 4)
    d = Location(10, 20)
    e = Location(5, 15)

    def test_eq(self):
        w = Location(0, 4, "FWD")
        x = Location(0, 4, "FWD")
        y = Location(0, 4, "REV")
        z = Location(3, 7, "FWD")

        assert x == x
        assert w == x
        assert x != y
        assert x != z

    def test_contains(self):
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e

        # testing contains on self should return True
        assert Location.contains(a, a) == True

        # Truly contained should return True ...
        assert Location.contains(a, b) == True
        # ... even with the same start index
        assert Location.contains(a, c) == True
        # ... but not vice versa (order matters!!)
        assert Location.contains(b, a) == False
        assert Location.contains(c, a) == False

        # Overlaped, but not contained, should return False
        assert Location.contains(a, e) == False
        assert Location.contains(e, a) == False

        # No overlap at all should return False
        assert Location.contains(a, d) == False
        assert Location.contains(d, a) == False

    def test_overlaps(self):
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        e = self.e

        # same location should overlap
        assert Location.overlaps(a, a) == True

        # sequential locations should not overlap
        assert Location.overlaps(a, d) == False
        assert Location.overlaps(d, a) == False

        # overlapping locations should return True
        assert Location.overlaps(a, e) == True
        assert Location.overlaps(e, a) == True

        # contained location should return True ...
        assert Location.overlaps(a, b) == True
        assert Location.overlaps(b, a) == True
        # ... even with the same start index
        assert Location.overlaps(a, c) == True
        assert Location.overlaps(c, a) == True

    def test_find_overlaps(self):
        locations = [self.a, self.b, self.c, self.d, self.e]

        assert Location.find_overlaps(locations) == [
            [1, 1, 1, 0, 1],
            [1, 1, 0, 0, 1],
            [1, 0, 1, 0, 0],
            [0, 0, 0, 1, 1],
            [1, 1, 0, 1, 1],
        ]

    # TODO: write test for to_slice()
    def test_to_slice(self):
        loc = self.b
        assert loc.to_slice() == slice(4, 7, 1)


class TestPart:
    def test_eq(self):
        # TODO: write test
        pass

    def test_slice(self):
        dna = DNA("ATCGAATTCCGG")
        part1 = Part(seq=dna, name="part1", location=Location(2, 8, "FWD"))
        part2 = Part(seq=dna, name="part2", location=Location(4, 10, "REV"))

        assert part1[2:4].seq == DNA("AA")
        assert part1[2:4] == Part(seq=dna, name="part1_subset",
                                  location=Location(4, 6))
        assert part1[2:4] == part1[Location(2, 4, "FWD")]

        assert part2[2:5].seq == DNA("GAA")
        assert part2[2:5] == Part(seq=dna, name="part2_subset",
                                  location=Location(6, 9, "REV"))
        assert part2[2:5] == part2[Location(2, 5, "REV")]

    def test_DNA_integration(self):
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

        # test reassignment from Part
        part2.seq = "AT"
        assert part1.seq == "ATCG"
        assert part1.location == Location(0, 4)

        assert part2.seq == "AT"
        assert part2.location == Location(4, 6)

        assert part3.seq == "CCGG"
        assert part3.location == Location(6, 10)

        assert part4.seq == "CGATCC"
        assert part4.location == Location(2, 8)

        part2.seq = "AATT"
        assert part1.seq == "ATCG"
        assert part1.location == Location(0, 4)

        assert part2.seq == "AATT"
        assert part2.location == Location(4, 8)

        assert part3.seq == "CCGG"
        assert part3.location == Location(8, 12)

        assert part4.seq == "CGAATTCC"
        assert part4.location == Location(2, 10)


if __name__ == '__main__':
    TestPart().test_DNA_integration()
