from synbio.utils import *
from synbio.utils import Location


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

    def test_lt(self):
        a = self.a
        b = self.b
        c = self.c

        assert (b < c) == False
        assert (c < b) == True

        assert (a < b) == False
        assert (a < c) == False

    def test_gt(self):
        a = self.a
        b = self.b
        c = self.c

        assert (b > c) == True
        assert (c > b) == False

        assert (a > b) == False
        assert (a > c) == False

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
