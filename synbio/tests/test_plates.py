from synbio.plates import *
from synbio.units import unit_registry as u


class TestWell:
    def test_init(self):
        pass

    def test_volume(self):
        # getter

        # setter
        pass


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

    def test_dict(self):
        # TODO: write
        # getter

        # setter
        pass

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
