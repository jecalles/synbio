from __future__ import annotations

import logging
from collections import UserDict
from typing import Callable, Dict, Iterable, List, Set, Union

import numpy as np
import pint

from synbio.units import QuantityType

__all__ = [
    # dataclasses
    "Reagent", "Recipe", "RecipeType",
    "Mixture",
    # functions
    "add_recipes", "calculate_reagent_volumes", "flatten_reagents", "get_recipe",
    # constants
    "reagent_registry", "PURE_reagents", "PURE",
]


class ReagentRegistry(UserDict):
    def __init__(self, name: str):
        super().__init__()
        self.name = name

    def __repr__(self):
        cls_name = self.__class__.__name__
        instance_name = self.name
        return f"{cls_name}({instance_name})"

    @staticmethod
    def __get_key(item: ReagentType) -> str:
        key = item.name if isinstance(item, Reagent) else item
        return key

    def __getitem__(self, item: ReagentType) -> Reagent:
        key = self.__get_key(item)
        return self.data[key]

    def __contains__(self, item: ReagentType) -> bool:
        if isinstance(item, Reagent):
            val = item in self.reagents
        else:
            val = item in self.data
        return val

    @property
    def reagents(self) -> Set[Reagent]:
        return set(self.values())

    def add_reagent(self, rgt: Reagent) -> None:
        if rgt in self.reagents:
            logging.warning(f"{rgt.name} already in registry")
        else:
            self[rgt.name] = rgt


reagent_registry = ReagentRegistry("synbio default")


class Reagent:
    """
    Mixtures can be made that represent combinations of contents (simple
    reagents or mixtures of reagents themselves). amounts in each recipe
    represent "parts" and not absolute volumes, e.g.,

    PURE.recipe = {
        'Sol A'     : 4,
        'Sol B'     : 3,
        'DNA'       : 1,
        'RNase Inh' : 1,
        'H2O'       : 1,
    }
    """

    def __init__(
            self,
            name: str = "NullReagent",
            registry: ReagentRegistry = reagent_registry
    ):
        self.name = name
        self._registry = registry

        self._registry.add_reagent(self)

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.name}')"

    def __eq__(self, other: Reagent) -> bool:
        same_registry = (self._registry == other._registry)
        same_name = (self.name == other.name)
        return same_registry and same_name

    def __hash__(self):
        hashables = (self.name, id(self._registry))
        return hash(hashables)


ReagentType = Union[Reagent, str]


class Recipe(UserDict):
    def __init__(
            self,
            recipe: RecipeType,
            registry: ReagentRegistry = reagent_registry
    ):
        super().__init__()
        self.registry = registry
        for key, quant in recipe.items():
            rgt = self.__get_reagent(key)
            self[rgt] = quant

    def __get_reagent(self, item: ReagentType) -> Reagent:
        if isinstance(item, Reagent):
            rgt = item
        elif isinstance(item, str):
            rgt = self.registry[item]
        else:
            raise TypeError(f"{item} is not ReagentType")
        return rgt

    def __getitem__(self, item: ReagentType) -> QuantityType:
        rgt = self.__get_reagent(item)
        return self.data[rgt]

    def __setitem__(self, key: ReagentType, value: QuantityType) -> None:
        # type checking
        if not isinstance(key, ReagentType):
            raise TypeError(f"{key} is not a Reagent")

        if not isinstance(value, QuantityType):
            raise TypeError(f"{value} is not a Quantity")

        rgt = self.__get_reagent(key)

        # delete item if value < 0
        mag = value.magnitude if isinstance(value, pint.Quantity) else value
        if mag <= 0:
            # delete if there, ignore if not
            if rgt in self:
                self.data.pop(rgt)
            else:
                pass
        # otherwise, update value
        else:
            self.data[rgt] = value

    @staticmethod
    def _merge_dict(
            dict1: RecipeType, dict2: RecipeType,
            op: Callable[[QuantityType, QuantityType], QuantityType]
    ) -> RecipeType:
        dict1 = Recipe(dict1)
        dict2 = Recipe(dict2)
        dict3 = {**dict1, **dict2}

        for key, value in dict3.items():
            if key in dict1 and key in dict2:
                dict3[key] = op(dict1[key], dict2[key])

        rcp = Recipe(dict3)
        return rcp

    def __add__(self, other: RecipeType):
        add = lambda x, y: x + y
        return self._merge_dict(self, other, add)

    def __radd__(self, other: RecipeType):
        return self + other

    def __sub__(self, other: RecipeType) -> Recipe:
        add = lambda x, y: x - y
        return self._merge_dict(self, other, add)

    def __rsub__(self, other: RecipeType):
        return self - other

    def __mul__(self, scalar: float) -> Recipe:
        return Recipe({
            rgt: quant * scalar
            for rgt, quant in self.items()
        })

    def __rmul__(self, scalar: float) -> Recipe:
        return self * scalar

    def __truediv__(self, scalar: float) -> Recipe:
        return Recipe({
            rgt: quant / scalar
            for rgt, quant in self.items()
        })

    def __eq__(self, other: RecipeType):
        this_flat = self.flat
        that_flat = other.flat
        same_rgts = this_flat.keys() == that_flat.keys()
        quants = [
            (this_flat[rgt], that_flat[rgt])
            for rgt in this_flat.keys()
        ]
        close_quants = np.allclose(*zip(*quants))
        return same_rgts and close_quants

    @property
    def total(self) -> QuantityType:
        return sum(self.data.values())

    @property
    def normalized(self) -> Recipe:
        return self / self.total

    def _flatten(self) -> List[Recipe]:
        """Recursive helper method for .flat"""
        simple_reagents = Recipe({
            rgt: quant
            for rgt, quant in self.items()
            if type(rgt) == Reagent
        })
        complex_reagents = [
            rcp
            for mix, quant in self.items() if type(mix) == Mixture
            for rcp in (mix.recipe.normalized * quant)._flatten()
        ]
        return [simple_reagents] + complex_reagents

    @property
    def flat(self) -> Recipe:
        """should return a new recipe with solely pure reagents"""
        return sum(self._flatten(), {})


RecipeType = Union[Recipe, Dict[ReagentType, QuantityType]]


class Mixture(Reagent):
    def __init__(
            self, name: str,
            recipe: RecipeType,
            registry: ReagentRegistry = reagent_registry
    ):
        if not all(rgt in registry for rgt in recipe.keys()):
            raise ValueError()

        self.recipe = Recipe(recipe)
        super().__init__(name, registry)

    # @property
    # def recipe(self) -> Dict[Reagent, float]:
    #     """
    #     A "recipe" is a dictionary that maps component contents to their
    #     relative ratios in the resulting Mixture. Pure Reagents
    #     return {self: 1}. All units are nondimensional.
    #     """
    #     if self._recipe is None:
    #         rec = {self: 1}
    #     else:
    #         rec = self._recipe
    #
    #     return rec
    #
    # @property
    # def reagents(self) -> List[Reagent]:
    #     return [r for r in self.recipe.keys()]


def add_recipes(recipes: Iterable[RecipeType]):
    return sum(recipes, {})

def calculate_reagent_volumes(
        rgt: Reagent, total_vol: pint.Quantity, flat: bool = True
) -> Recipe:
    """this needs to recursively walk the mixture tree"""
    if isinstance(rgt, Mixture):
        raw_rcp = rgt.recipe.flat if flat else rgt.recipe
        rcp = raw_rcp.normalized
    else:
        rcp = Recipe({rgt: 1})

    return rcp * total_vol


def flatten_reagents(
        reagents: Iterable[Reagent],
) -> Set[Reagent]:
    """Function that flattens a list of Reagents (can include mixtures) to
    extract all unique Reagents in recipe"""

    # handle None input
    reagents = [] if reagents is None else reagents

    # extract Reagents
    simple_reagents = set(rgt for rgt in reagents if type(rgt) == Reagent)

    # flatten recipes for
    subcomponents = set(
        sub for rgt in reagents if isinstance(rgt, Mixture)
        for sub in rgt.recipe.flat.keys()
    )
    return {*simple_reagents, *subcomponents}

def get_recipe(reagent: Reagent, flat: bool = False) -> Recipe:
    if isinstance(reagent, Mixture):
        rcp = reagent.recipe.flat if flat else reagent.recipe
    elif type(reagent) == Reagent:
        rcp = Recipe({reagent: 1})
    else:
        raise TypeError(f"{reagent} is not a Reagent")
    return rcp


# Constants
PURE_reagents = {
    name: Reagent(name)
    for name in "Sol-A Sol-B DNA RNase-Inh H2O".split()
}

PURE = Mixture(name="PURE", recipe={
    'Sol-A': 4,
    'Sol-B': 3,
    'DNA': 1,
    'RNase-Inh': 1,
    'H2O': 1
})

if __name__ == "__main__":
    new_rgts = {
        name: Reagent(name)
        for name in "A B C D".split()
    }
    test_recipe = Recipe({
        "A": 2,
        "B": 3,
        Mixture("C+D", {"C": 2, "D": 1}): 3,
        Mixture("A+'C+D'", {"A": 1, "C+D": 1}): 2,
        "D": 1
    })
    flat = test_recipe.flat
    assert flat == Recipe({
        "A": 2 + 1,
        "B": 3,
        "C": (2 / 3) * 3 + (2 / 3) * (1 / 2) * 2,
        "D": (1 / 3) * 3 + (1 / 3) * (1 / 2) * 2 + 1
    })
