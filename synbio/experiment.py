from typing import List, Dict, Set
from datetime import date
from dataclasses import dataclass, field
from functools import reduce

import pint

from synbio.units import unit_registry as u
from synbio.reagents import Reagent

__all__ = [
    # Dataclasses
    "Condition", "Experiment", "Data"
]


@dataclass
class Condition:
    name: str
    content: Reagent
    replicates: int = 3
    volume: pint.Quantity = 10 * u.uL


    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        recipe = self.content.recipe

        total = sum(recipe.values())
        vol_factor = self.volume / total * self.replicates

        return {
            reagent: num * vol_factor
            for reagent, num in recipe.items()
        }


class Data: pass

class Experiment:
    def __init__(
            self, name: str,
            conditions: List[Condition] = None,
            data: Data = None,
            meta: dict = None
    ):
        today = date.today()
        if meta is None:
            meta = {}

        self.name = name
        self.conditions = conditions
        self.data = data
        self.meta = meta


    @property
    def reagents(self):
        return {
            reagent
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


    def add_conditions(
            self, contents: Dict[str, Reagent],
            **cond_kwargs
    ) -> "self":
        """of dubious use, and questionable validity"""
        self.conditions = [
            Condition(
                name=name, content=cont,
                **cond_kwargs

            ) for name, cont in contents.items()
        ]

        return self
