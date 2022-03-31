import itertools
from dataclasses import dataclass
from functools import cached_property
from itertools import count
from math import prod
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd
import pint

from synbio.reagents import Reagent, add_recipes, calculate_reagent_volumes
from synbio.units import unit_registry as u

__all__ = [
    # Dataclasses
    "Well", "Plate",
    # Functions
    "make_96_well",
    "make_384_well", "make_384_ldv_well",
    "make_1536_ldv_well",
    "make_plate_like"
]


@dataclass
class Well:
    name: str
    content: Reagent
    max_vol: pint.Quantity
    dead_vol: pint.Quantity
    location: str = None

    vol: pint.Quantity = 0 * u.uL

    def __repr__(self):
        if self.location is None:
            var_str = f"('{self.name}', {self.volume})"
        else:
            var_str = f" @ {self.location}('{self.name}', {self.volume})"

        return f"{self.__class__.__name__}{var_str}"

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
    def reagents(self) -> Set[Reagent]:
        return {rgt for rgt in self.content.recipe.keys()}

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        # PREVIOUS IMPLEMENTATION
        # vol_factor = self.volume / sum(self.content.recipe.values())
        # return {
        #     rgt: num * vol_factor
        #     for rgt, num in self.content.recipe.items()
        # }
        return calculate_reagent_volumes(self.content, self.volume)

PlateLocationType = Union[Tuple[int, int], int, str]

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
                content=Reagent(),
                max_vol=max_vol,
                dead_vol=dead_vol,
                location=f"{rows[i]}{cols[j]}"
            )
        self.array = array

        # self._full_wells = self.__full_wells()

    def __repr__(self) -> str:
        cls_ = self.__class__.__name__
        name = self.name
        shp = self.shape
        n_full = len(self.full_wells)

        vols = self.well_volumes

        return f'{cls_}(name="{name}", shape={shp}, {n_full} full wells, ' \
               f'max_vol={vols["max_vol"]}), dead_vol={vols["dead_vol"]})'

    def __getitem__(self, key: PlateLocationType):
        if isinstance(key, str):
            val = self.dict[key]
        elif isinstance(key, tuple):
            i, j = key
            val = self.array[i, j]
        elif isinstance(key, int):
            val = self.array[key]

        else:
            raise TypeError(f"key is not a str, int or tuple of int")

        return val

    def __setitem__(self, key: PlateLocationType, value: Well):
        if not isinstance(value, Well):
            raise TypeError(f"value is not a Well")
        if isinstance(key, str):
            well = self.dict[key]
        elif isinstance(key, tuple):
            i, j = key
            well = self.array[i, j]
        else:
            raise TypeError(f"key is not a str or tuple of int")

        well.name = value.name
        well.content = value.content
        well.volume = value.volume

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
        """
        dict attributes
            'max_vol': maximum volume of well,
            'dead_vol': volume lost at bottom of well,
            'working_vol': total useable volume
        """
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
    def contents(self) -> Dict[str, Reagent]:
        return {
            well.content.name: well.content
            for well in self._full_wells.values()
        }

    @property
    def wells_by_content(self) -> Dict[Reagent, List[Well]]:
        contents = set(self.contents.values())
        d = {rgt: [] for rgt in contents}

        for well in self.dict.values():
            if well.content in contents:
                d[well.content].append(well)

        return d

    @property
    def content_volumes(self) -> Dict[Reagent, pint.Quantity]:
        return {
            rgt: sum(
                well.volume for well in well_list
            ) for rgt, well_list in self.wells_by_content.items()
        }

    @property
    def reagents(self) -> Dict[str, Reagent]:
        """
        differs from self.contents, in as much as it separates Mixtures
        into the Reagents that compose them
        """
        return {
            rgt.name: rgt
            for content in self.contents.values()
            for rgt in content.recipe.keys()
        }

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        # rgt_vols_list_by_content = {
        #     content: [well.reagent_volumes for well in well_list]
        #     for content, well_list in self.wells_by_content.items()
        # }
        # rgt_vols_by_content = {
        #     content: add_recipes(rgt_vol_dict_list)
        #     for content, rgt_vol_dict_list in rgt_vols_list_by_content.items()
        # }
        # return add_recipes(rgt_vols_by_content.values())
        recipe_list: List[Dict[Reagent, pint.Quantity]] = [
            well.reagent_volumes
            for well in self.full_wells.values()
        ]
        return add_recipes(recipe_list)

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

    @cached_property
    def _full_wells(self) -> Dict[str, Well]:
        """
        Alternate to self.full_wells. Returns any Wells with nonzero
        contents. Will return wells with well.volume == 0 if well.content !=
        Reagent()
        """
        return {
            name: well
            for name, well in self.dict.items()
            if well.content != Reagent()
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

    def load_plate_map(
            self,
            filepath: str,
            reagents: Dict[str, Reagent] = None
    ) -> None:
        df = pd.read_csv(filepath)
        index = df.values[:, 0]
        name_dataframe = df.set_index(index).drop("Rows/Cols", axis=1)


        if reagents is None:
            reagents = {
                name: Reagent(name)
                for row in name_dataframe.values
                for name in row
                if name != "-"
            }

        counters = {
            name: count(1)
            for name in reagents.keys()
        }

        # Reassign wells
        for i, row in enumerate(name_dataframe.values):
            for j, name in enumerate(row):
                if name == "-":
                    continue
                reagent = reagents[name]

                well = self.array[i, j]
                well.content = reagent
                well.name = f"{name}_{next(counters[name])}"




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

def make_plate_like(plate: Plate, name=None) -> Plate:
    if name is None:
        name = f"like_{plate.name}"
    return Plate(
        name=name,
        shape=plate.shape,
        dead_vol=plate.well_volumes["dead_vol"],
        max_vol = plate.well_volumes["max_vol"]
    )


def diff_plates(plate1: Plate, plate2: Plate) -> Plate:
    if not plate1.shape == plate2.shape:
        raise ValueError("compared plates must have the same shape")
    if not plate1.well_volumes == plate2.well_volumes:
        raise ValueError("compared plates must have the same well volumes")
