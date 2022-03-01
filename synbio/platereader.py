from typing import Dict, List
from dataclasses import dataclass, field
from functools import reduce

import pint

from synbio.plates import Plate
from synbio.experiment import *
from synbio.reagents import Reagent
from synbio.units import unit_registry as u


__all__ = [
    "PlateReaderCondition", "PlateReaderData", "PlateReaderExperiment"
]

@dataclass
class PlateReaderCondition(Condition):
    content: Reagent = None
    volume: pint.Quantity = 10 * u.uL
    replicate: int = 3
    plate: Plate = field(default=None, repr=False)

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        recipe = self.content.recipe

        total = sum(recipe.values())
        vol_factor = self.volume / total * self.replicate

        return {
            reagent: num * vol_factor
            for reagent, num in recipe.items()
        }


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


    @property
    def reagents(self):
        return {
            reagent.name: reagent
            for cond in self.conditions
            for reagent in cond.content.recipe.keys()
        }


    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        # list of dicts describing reagents for each condition
        adj_recipes = [c.reagent_volumes for c in self.conditions]

        def merge_dict(dict1, dict2):
            dict3 = {**dict1, **dict2}
            for key, value in dict3.items():
                if key in dict1 and key in dict2:
                    dict3[key] = dict1[key] + dict2[key]
            return dict3

        return reduce(merge_dict, adj_recipes)

