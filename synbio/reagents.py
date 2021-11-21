from typing import Dict, List
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

    def __hash__(self):
        return hash(self.name)


@dataclass
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
    recipe: Dict[Reagent, float]

    @property
    def components(self) -> List[Reagent]:
        return list(self.recipe.keys())

# Constants
pure_reagents = {
    name: Reagent(name)
    for name in "Sol_A Sol_B DNA RNase_Inh H20".split()
}
PURE = Mixture(name="PURE", recipe={
    pure_reagents['Sol_A']: 4,
    pure_reagents['Sol_B']: 3,
    pure_reagents['DNA']: 1,
    pure_reagents['RNase_Inh']:1,
    pure_reagents['H20']: 1
})