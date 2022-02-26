from typing import List, Dict
from datetime import date
from dataclasses import dataclass
from functools import reduce

import pint

from synbio.units import unit_registry as u
from synbio.reagents import Mixture, Reagent, ReagentRegistry

__all__ = [
    # Dataclasses
    "Condition", "Experiment", "Data"
]


@dataclass
class Condition:
    name: str
    content: Reagent
    reagent_registry: ReagentRegistry
    replicates: int = 3
    volume: pint.Quantity = 10 * u.uL


    def __post_init__(self):
        if not self._check_registry():
            raise ValueError("Content not found in ReagentRegistry")

    @property
    def reagent_volumes(self) -> Dict[Reagent, pint.Quantity]:
        recipe = self.content.recipe

        total = sum(recipe.values())
        vol_factor = self.volume / total * self.replicates

        return {
            reagent: num * vol_factor
            for reagent, num in recipe.items()
        }

    def _check_registry(self) -> bool:
        if isinstance(self.content, Mixture):
            chk = all([
                r in self.reagent_registry.reagents
                for r in self.content.reagents
            ])


        else:
            chk = self.content in self.reagent_registry

        return chk


@dataclass
class Experiment:
    name: str
    reagent_registry: ReagentRegistry
    conditions: List[Condition] = None
    date: date = date.today()

    def __post_init__(self):
        self._check_conditions()


    def add_conditions(
            self, contents: Dict[str, Reagent],
            **cond_kwargs
            #replicates: int = None, volume: pint.Quantity = None
    ) -> "self":
        self.conditions = [
            Condition(
                name=name, content=cont,
                reagent_registry=self.reagent_registry,
                **cond_kwargs

            ) for name, cont in contents.items()
        ]

        self._check_conditions()
        return self

    def _check_conditions(self) -> None:
        if self.conditions is not None:
            all_same_registry = all(
                cond.reagent_registry == self.reagent_registry
                for cond in self.conditions
            )
        else:
            all_same_registry = True

        if not all_same_registry:
            raise ValueError("Conditions using multiple registries")



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

class Data: pass
