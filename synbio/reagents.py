from typing import Dict, List

from synbio.interfaces import HashableMixin

__all__ = [
    # dataclasses
    "Reagent", "Mixture",
    # functions
    "mixture_from_recipes",
    # constants
    "PURE_reagents", "PURE"
]


class Reagent(HashableMixin):
    def __init__(self, name: str):
        self.name = name
        self._recipe = None

    def _comparables(self) -> List[str]:
        return ["name"]

    # def __hash__(self):
    #     return hash((self.name, self._recipe))

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.name}')"

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

    # def __repr__(self) -> str:
    """
    old repr for Mixture. more detailed, but leads to clutter when printing
    """

    #     cls_ = self.__class__.__name__
    #     name = self.name
    #     rcp = {
    #         rgnt.name: amt
    #         for rgnt, amt in self.recipe.items()
    #     }
    #     return f'{cls_}(name="{name}", recipe="{rcp}")'

    @property
    def reagents(self) -> List[Reagent]:
        return [r for r in self.recipe.keys()]


def mixture_from_recipes(
        recipe_dict: Dict[str, Dict[Reagent, int]]
) -> Dict[str, Mixture]:
    return {
        name: Mixture(name=name, recipe=rec)
        for name, rec in recipe_dict.items()
    }


# Constants
PURE_reagents = {
    name: Reagent(name)
    for name in "Sol-A Sol-B DNA RNase-Inh H20".split()
}

PURE = Mixture(name="PURE", recipe={
    PURE_reagents['Sol-A']: 4,
    PURE_reagents['Sol-B']: 3,
    PURE_reagents['DNA']: 1,
    PURE_reagents['RNase-Inh']: 1,
    PURE_reagents['H20']: 1
})
