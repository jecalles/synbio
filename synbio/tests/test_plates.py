import synbio.tests.utils as testutils
from synbio.plates import *
from synbio.reagents import Recipe, reagent_registry as r
from synbio.units import unit_registry as u


class TestPlates:
    def test_init(self):
        name = "test_plate"
        shape = (3, 4)  # not a real plate
        max_vol = 30 * u.uL
        dead_vol = 5 * u.uL

        Plate(
            name=name,
            shape=shape,
            max_vol=max_vol,
            dead_vol=dead_vol
        )

    def test___getitem__(self):
        plate = make_96_well()
        assert plate["E3"] == plate[4, 2]

    def test___setitem__(self):
        plate = make_96_well()
        new_well = Well(
            name="test_well",
            dead_vol=0 * u.uL, max_vol=2000 * u.uL,
            content=r["PURE"], vol=10 * u.uL
        )
        plate["C4"] = new_well
        assert plate["C4"].name == "test_well"
        assert plate["C4"].content == r["PURE"]
        assert plate["C4"].volume == 10 * u.uL
        assert plate["C4"].location == "C4"

        plate[4, 2] = new_well
        assert plate[4, 2].name == "test_well"
        assert plate[4, 2].content == r["PURE"]
        assert plate[4, 2].volume == 10 * u.uL
        assert plate[4, 2].location == "E3"

    def test___add__(self):
        """
        test1 = [
            [H2O, xxx, xxx, xxx],
            [xxx, xxx, H2O, xxx],
            [xxx, xxx, xxx, H2O],
        ]
        test2 = [
            [xxx,  xxx, xxx,  PURE],
            [xxx,  xxx, PURE, xxx],
            [PURE, xxx, xxx,  xxx],
        ]

        test1 + test 2 = [
            [H2O,  xxx, xxx,      PURE],
            [xxx,  xxx, PURE+H2O, xxx],
            [PURE, xxx, xxx,      H2O],
        ]

        """
        shape = (3, 4)  # not a real plate
        max_vol = 30 * u.uL
        dead_vol = 5 * u.uL

        test1 = Plate(
            name="test_plate1",
            shape=shape,
            max_vol=max_vol,
            dead_vol=dead_vol
        )
        test1.fill_wells(r["H2O"], 10 * u.uL, ["A1", (1, 2), "C4"])
        test2 = Plate(
            name="test_plate2",
            shape=shape,
            max_vol=max_vol,
            dead_vol=dead_vol
        )
        test2.fill_wells(r["PURE"], 10 * u.uL, ["C1", (1, 2), "A4"])

        res = test1 + test2

        pure_H2O = (r["PURE"].recipe * 10 * u.uL) + Recipe({r["H2O"]: 10 * u.uL})
        assert res["B3"].content.recipe == pure_H2O
        assert res["B3"].volume == 20 * u.uL

        assert res["A1"].content == r["H2O"]
        assert res["A1"].volume == 10 * u.uL

        assert res["C4"].content == r["H2O"]
        assert res["C4"].volume == 10 * u.uL

        assert res["A4"].content == r["PURE"]
        assert res["A4"].volume == 10 * u.uL

        assert res["C1"].content == r["PURE"]
        assert res["C1"].volume == 10 * u.uL

    def test___sub__(self):
        raise testutils.TestNotImplemented

    def test_well_volumes(self):
        plate = make_384_ldv_well()

        # getter
        assert plate.well_volumes == {
            'max_vol': 14 * u.uL,
            'dead_vol': 6 * u.uL,
            'working_vol': 8 * u.uL
        }

        # setter
        new_vols = {
            'max_vol': 200 * u.uL,
            'dead_vol': 50 * u.uL,
            'working_vol': 150 * u.uL
        }
        plate.well_volumes = new_vols

        test_well = plate.array[4][2]
        assert test_well.max_vol == 200 * u.uL
        assert test_well.dead_vol == 50 * u.uL
        assert test_well.working_vol == 150 * u.uL

    def test_contents(self):
        raise testutils.TestNotImplemented

    def test_wells_by_content(self):
        raise testutils.TestNotImplemented

    def test_content_volumes(self):
        raise testutils.TestNotImplemented

    def test_reagents(self):
        raise testutils.TestNotImplemented

    def test_reagent_volumes(self):
        raise testutils.TestNotImplemented

    def test_empty_wells(self):
        raise testutils.TestNotImplemented

    def test_full_wells(self):
        raise testutils.TestNotImplemented
