from typing import Dict
from dataclasses import dataclass

__all__ = [
    # dataclasses
    "Reagent", "Mixture",
    # constants
    "PURE"
]


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


# Constants
pure_reagents = {
    name: Reagent(name)
    for name in "Sol-A Sol-B DNA RNase-Inh H20".split()
}
PURE = Mixture(name="PURE", recipe={
    pure_reagents['Sol-A']: 4,
    pure_reagents['Sol-B']: 3,
    pure_reagents['DNA']: 1,
    pure_reagents['RNase-Inh']:1,
    pure_reagents['H20']: 1
})