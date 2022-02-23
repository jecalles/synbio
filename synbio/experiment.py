from typing import List, Dict
from datetime import date
from dataclasses import dataclass
from functools import reduce

import pint

from synbio.units import unit_registry as u
from synbio.reagents import Mixture, Reagent

__all__ = [
    # Dataclasses
    "Condition", "ExperimentDesign",
    # functions
    "calc_reagents"
]


@dataclass
class Condition:
    name: str
    content: Mixture
    replicates: int = 3
    volume: pint.Quantity = 1 * u.uL

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        recipe = self.content.recipe

        total = sum(recipe.values())
        vol_factor = self.volume / total * self.replicates

        return {
            reagent: num * vol_factor
            for reagent, num in recipe.items()
        }


@dataclass
class ExperimentDesign:
    name: str
    conditions: List[Condition]
    date: date = date.today()


def calc_reagents(conditions: List[Condition]) -> Dict[Reagent, pint.Quantity]:
    """
    TODO: 1. write docstring 2. move this to plates.Well
    """


    # list of dicts describing reagents for each condition
    adj_recipes = [c.reagent_volumes for c in conditions]

    def merge_dict(dict1, dict2):
        dict3 = {**dict1, **dict2}
        for key, value in dict3.items():
            if key in dict1 and key in dict2:
                dict3[key] = dict1[key] + dict2[key]
        return dict3

    return reduce(merge_dict, adj_recipes)
