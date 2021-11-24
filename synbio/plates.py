import itertools
from dataclasses import dataclass
from datetime import date
from math import prod
from typing import Dict, List, Optional, Tuple, TypeVar

import numpy as np
import pandas as pd
import pint

from synbio.reagents import Mixture, Reagent
from synbio.units import unit_registry as u

__all__ = [
    # Dataclasses
    "Well", "Plate",
    "ExperimentalCondition", "Experiment",
    # Functions
    "make_96_well",
    "make_384_well", "make_384_ldv_well",
    "make_1536_ldv_well"
]


@dataclass
class Well:
    name: str
    content: Reagent
    max_vol: pint.Quantity
    dead_vol: pint.Quantity

    vol: pint.Quantity = 0 * u.uL

    def __str__(self):
        return self.name

    @property
    def volume(self) -> pint.Quantity:
        return np.round(self.vol, decimals=4)

    @volume.setter
    def volume(self, v: pint.Quantity) -> None:
        vol = v.magnitude * u.Unit(v.units)
        if vol > self.max_vol:
            raise ValueError(
                f"set volume exceeds max well volume ({self.max_vol})"
            )
        elif vol < 0:
            raise ValueError(f"set volume is less than zero {vol}")

        self.vol = vol

    @property
    def working_vol(self) -> pint.Quantity:
        """
        A property that returns the maximum transferable volume in a Well,
        given its max_vol and dead_vol. Working volume is a property of the
        geometry of the Well
        """
        return max(
            np.round(self.max_vol - self.dead_vol, decimals=4),
            0 * u.uL
        )

    @property
    def available_vol(self) -> pint.Quantity:
        """
        A property that returns the total volume in the well that can be
        transferred, given that wells have dead volume.
        """
        return max(
            np.round(self.volume - self.dead_vol, decimals=4),
            0 * u.uL
        )

    @property
    def recipe(self) -> Dict[Reagent, pint.Quantity]:
        rec = self.content.recipe

        total = sum(rec.values())
        vol_factor = self.volume / total

        return {
            reagent: np.round(num * vol_factor, decimals=4)
            for reagent, num in rec.items()
        }


class Plate:
    def __init__(
            self, name: str, shape: Tuple[int, int],
            max_vol: pint.Quantity, dead_vol: pint.Quantity
    ):
        self.name = name

        rows = self.get_row_names(shape[0])
        cols = self.get_col_names(shape[1])

        self.name_map = {
            f"{r}{c}": (i, j)
            for i, r in enumerate(rows)
            for j, c in enumerate(cols)
        }

        array = np.empty(shape, Well)
        for well_name, (i, j) in self.name_map.items():
            array[i][j] = Well(
                name=well_name,
                content=Reagent("empty"),
                max_vol=max_vol,
                dead_vol=dead_vol
            )
        self.array = array

    def __repr__(self) -> str:
        cls_ = self.__class__.__name__
        name = self.name
        shp = self.shape
        n_full = len(self.full_wells)

        vols = self.well_volumes

        return f'{cls_}(name="{name}", shape={shp}, {n_full} full wells, ' \
               f'max_vol={vols["max_vol"]}), dead_vol={vols["dead_vol"]})'

    @property
    def shape(self) -> Tuple[int, int]:
        return self.array.shape

    @property
    def dict(self) -> Dict[str, Well]:
        return {
            well_name: self.array[i][j]
            for well_name, (i, j) in self.name_map.items()
        }

    @dict.setter
    def dict(self, well_dict) -> None:
        for well_name, well in well_dict.items():
            i, j = self.name_map[well_name]
            self.array[i][j] = well

    @property
    def well_volumes(self) -> Dict[str, pint.Quantity]:
        volumes_set = {
            (well.max_vol, well.dead_vol)
            for well in self.dict.values()
        }
        if len(volumes_set) != 1:
            raise AttributeError("plate wells have inconsistent volumes!")

        max_vol, dead_vol = volumes_set.pop()
        working_vol = max_vol - dead_vol

        if working_vol < 0:
            raise AttributeError(f"dead volume {dead_vol} is larger than "
                                 f"working volume {working_vol}")

        return {
            'max_vol': max_vol,
            'dead_vol': dead_vol,
            'working_vol': working_vol
        }

    @well_volumes.setter
    def well_volumes(self, volumes: Dict[str, pint.Quantity]) -> None:
        if all(
                x not in {'max_vol', 'dead_vol', 'working_vol'}
                for x in volumes.keys()
        ):
            raise AttributeError("must input dict with the following keys: "
                                 "'max_vol', 'dead_vol', 'working_vol'")

        for well in self.dict.values():
            well.max_vol = volumes['max_vol']
            well.dead_vol = volumes['dead_vol']

        if self.well_volumes != volumes:
            raise AttributeError(
                f"something happened during volume reassingment. \n"
                f"new volumes = {volumes} \n"
                f"self.well_volumes = {self.well_volumes} \n"
            )

    @property
    def num_wells(self) -> int:
        return prod(self.shape)

    @property
    def empty_wells(self) -> Dict[str, Well]:
        return {
            name: well
            for name, well in self.dict.items()
            if well.volume == 0 * u.uL
        }

    @property
    def full_wells(self) -> Dict[str, Well]:
        empty = self.empty_wells.keys()
        return {
            name: well
            for name, well in self.dict.items()
            if name not in empty
        }

    @staticmethod
    def get_row_names(num_rows: int) -> List[str]:

        def name_generator():
            alphabet = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
            num_letters = 1

            while True:
                yield from (
                    ''.join(chars)
                    for chars in itertools.product(alphabet, repeat=num_letters)
                )
                num_letters += 1

        return list(itertools.islice(name_generator(), num_rows))

    @staticmethod
    def get_col_names(num_cols: int) -> List[str]:
        return list(str(i + 1) for i in range(num_cols))


PlateLocationType = TypeVar("PlateLocationType", Tuple[int, int], str)


@dataclass
class ExperimentalCondition:
    name: str
    content: Mixture
    replicates: int = 3
    volume: pint.Quantity = 1 * u.uL


@dataclass
class Experiment:
    name: str
    source_plate: Plate
    dest_plate: Plate
    conditions: List[ExperimentalCondition]
    name_map: Dict[str, str]
    data: pd.DataFrame = None
    date: date = date.today()


# Functions
def make_96_well(name: Optional[str] = None) -> Plate:
    if name is None:
        name = "96well"

    shape = (8, 12)
    max_vol = 200 * u.uL
    dead_vol = 40 * u.uL

    return Plate(name, shape, max_vol, dead_vol)


def make_384_well(name: Optional[str] = None) -> Plate:
    if name is None:
        name = "384well"

    shape = (16, 24)
    max_vol = 65 * u.uL
    dead_vol = 20 * u.uL

    return Plate(name, shape, max_vol, dead_vol)


def make_384_ldv_well(name: Optional[str] = None) -> Plate:
    if name is None:
        name = "384LDV"

    shape = (16, 24)
    max_vol = 14 * u.uL
    dead_vol = 6 * u.uL

    return Plate(name, shape, max_vol, dead_vol)


def make_1536_ldv_well(name: Optional[str] = None) -> Plate:
    if name is None:
        name = "1536LDV"

    shape = (32, 48)
    max_vol = 5.5 * u.uL
    dead_vol = 1 * u.uL

    return Plate(name, shape, max_vol, dead_vol)
