from __future__ import annotations

from dataclasses import dataclass
from itertools import count
from typing import Callable, Dict, List, Tuple

import pandas as pd
import pint

from synbio.experiment import Condition, Data, Experiment
from synbio.plates import Plate, make_96_well
from synbio.reagents import Reagent, add_recipes, calculate_reagent_volumes
from synbio.units import unit_registry as u

__all__ = [
    "PlateReaderCondition", "PlateReaderData", "PlateReaderExperiment",
    "load_plate_map"
]


class PlateReaderCondition(Condition):
    def __init__(
            self, content: Reagent = None,
            volume: pint.Quantity = 10 * u.uL,
            replicate: int = 3,
            plate: Plate = None
    ):
        self.content = content
        self.volume = volume
        self.replicate = replicate
        self.plate = plate

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        # PREVIOUS IMPLEMENTATION
        # recipe = self.content.recipe
        #
        # total = sum(recipe.values())
        # vol_factor = self.volume / total * self.replicate
        #
        # return {
        #     reagent: num * vol_factor
        #     for reagent, num in recipe.items()
        # }
        total_vol = self.volume * self.replicate
        return calculate_reagent_volumes(self.content, total_vol)


@dataclass
class PlateReaderData(Data):
    def to_csv(self):
        pass

    def subtract_background(self):
        pass

    def plot_results(self, print_flag: bool = False, plot_rc: dict = None):
        pass


class PlateReaderExperiment(Experiment):
    def __init__(
            self,
            conditions: List[PlateReaderCondition],
            data: PlateReaderData = None,
            *exp_args, **exp_kwargs,
    ):
        super().__init__(*exp_args, **exp_kwargs)
        self.conditions = conditions
        self.data = data

    @classmethod
    def from_plate_map(
            cls, filepath: str,
            cond_name_map: Callable[[str, Plate], PlateReaderCondition],
            plate: Plate = make_96_well()
    ) -> PlateReaderExperiment:
        conditions, _ = load_plate_map(filepath, cond_name_map, plate)
        return cls(conditions=conditions)

    @property
    def reagents(self):
        return {
            reagent.name: reagent
            for cond in self.conditions
            for reagent in cond.content.recipe.keys()
        }

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        # list of dicts describing contents for each condition
        recipe_list: List[Dict[Reagent, pint.Quantity]] = [
            c.reagent_volumes for c in self.conditions
        ]
        return add_recipes(recipe_list)


def load_plate_map(
        filepath: str,
        cond_name_map: Callable[[str, Plate], PlateReaderCondition],
        plate: Plate = make_96_well()
) -> Tuple[List[PlateReaderCondition], Plate]:
    """
    """
    df = pd.read_csv(filepath)
    index = df.values[:, 0]
    name_dataframe = df.set_index(index).drop("Rows/Cols", axis=1)

    conditions = {
        name: cond_name_map(name, plate)
        for row in name_dataframe.values
        for name in row
        if name != "-"
    }

    counters = {
        name: count(1)
        for name in conditions.keys()
    }

    # Reassign wells
    for i, row in enumerate(name_dataframe.values):
        for j, name in enumerate(row):
            if name == "-":
                continue
            replicate = next(counters[name])
            condition = conditions[name]
            condition.replicate = replicate

            well = plate.array[i, j]
            well.content = condition.content
            well.volume = condition.volume
            well.name = f"{name}_{replicate}"

    for cond in conditions.values():
        cond.plate = plate

    return list(conditions.values()), plate
