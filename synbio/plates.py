from typing import List, Tuple, Dict
from dataclasses import dataclass
import itertools

import numpy as np
import pint

from synbio.reagents import Reagent


__all__ = [
    # Dataclasses
    "Well", "Plate",
    # Functions
    "make_96_well",
    "make_384_well", "make_384_ldv_well",
    "make_1536_well", "make_1536_ldv_well"
]


@dataclass
class Well:
    name: str
    content: Reagent
    max_vol: pint.Quantity
    dead_vol: pint.Quantity

    _vol: pint.Quantity = 0 * pint.UnitRegistry().uL

    def __str__(self):
        return self.name

    @property
    def volume(self) -> pint.Quantity:
        return self._vol

    @volume.setter
    def volume(self, v: pint.Quantity) -> None:
        if v > self.max_vol:
            raise ValueError(
                f"set volume exceeds max well volume: ({self.max_vol})"
            )
        else:
            self._vol = v


class Plate:
    def __init__(
            self, name: str, shape:Tuple[int, int],
            max_vol: pint.Quantity, dead_vol: pint.Quantity
    ):
        self.name = name
        self.shape = shape
        self.max_vol = max_vol
        self.dead_vol = dead_vol

        rows = self.get_row_names(self.shape[0])
        cols = self.get_col_names(self.shape[1])
        self.name_map = {
            f"{r}{c}": (i, j)
            for i, r in enumerate(rows)
            for j, c in enumerate(cols)
        }

        array = np.empty(self.shape, Well)
        for well_name, (i, j) in self.name_map.items():
            array[i][j] = Well(
                name=well_name,
                content=Reagent("empty"),
                max_vol=self.max_vol,
                dead_vol=self.dead_vol
            )
        self.array = array

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


# Functions
def make_96_well() -> Plate:
    return


def make_384_well() -> Plate:
    return


def make_384_ldv_well() -> Plate:
    return


def make_1536_well() -> Plate:
    return


def make_1536_ldv_well() -> Plate:
    return
