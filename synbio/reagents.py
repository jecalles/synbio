from dataclasses import dataclass
from typing import Dict, List

__all__ = [
    # dataclasses
    "Reagent", "Mixture", "ReagentRegistry",
    # functions
    "mixtures_from_recipes",
    # constants
    "pure_registry", "PURE"
]

"""
TODO: make ReagentsRegistry class to track reagents during experiment
"""

@dataclass
class Reagent:
    name: str

    def __post_init__(self):
        self._recipe = None

    def __hash__(self):
        return hash((self.name, self._recipe))

    @property
    def recipe(self) -> Dict['OwnType', float]:
        """
        A "recipe" is a dictionary that maps component reagents to their
        relative ratios in the resulting Mixture. Pure Reagents
        return {self: 1}. All units are nondimensional.
        """
        if self._recipe is None:
            rec = {self: 1}
        else:
            rec = self._recipe

        return rec


class Mixture(Reagent):
    """
    Mixtures represent combinations of reagents (simple or mixtures
    themselves). amounts in each recipe represent "parts" and not absolute
    volumes, e.g.,

    PURE.recipe = {
        'Sol A'     : 4,
        'Sol B'     : 3,
        'DNA'       : 1,
        'RNase Inh' : 1,
        'H20'       : 1,
    }
    """

    def __init__(self, name: str, recipe: Dict[Reagent, float]):
        super().__init__(name)
        self._recipe = recipe

    def __repr__(self) -> str:
        cls_ = self.__class__.__name__
        name = self.name
        rcp = {
            rgnt.name: amt
            for rgnt, amt in self.recipe.items()
        }
        return f'{cls_}(name="{name}", recipe="{rcp}")'


    @property
    def reagents(self) -> List[Reagent]:
        return [r for r in self.recipe.keys()]


class ReagentRegistry(dict):
    @property
    def reagents(self) -> List[Reagent]:
        return [r for r in self.values()]

    def add_by_name(self, names: "str|List[str]") -> "self":
        """
        updates self in place from list of new reagent names. returns self
        """
        if not isinstance(names, list):
            names = [names]

        new_dict = {name: Reagent(name) for name in names}
        self.update(new_dict)

        return self


def mixtures_from_recipes(
        recipe_dict: Dict[str, Dict[Reagent, int]]
) -> Dict[str, Mixture]:
    return {
        name: Mixture(name=name, recipe=rec)
        for name, rec in recipe_dict.items()
    }


# Constants
pure_registry = (ReagentRegistry()
                 .add_by_name(
                    "Sol-A Sol-B DNA RNase-Inh H20".split()
                ))

PURE = Mixture(name="PURE", recipe={
    pure_registry['Sol-A']: 4,
    pure_registry['Sol-B']: 3,
    pure_registry['DNA']: 1,
    pure_registry['RNase-Inh']: 1,
    pure_registry['H20']: 1
})
