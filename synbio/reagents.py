from typing import Dict, List, Optional
from dataclasses import dataclass

import pint

__all__ = [
    "Reagent", "PureReagent", "Mixture"
]

@dataclass
class Reagent:
    name: str

    def __hash__(self):
        return hash(self.name)


@dataclass
class PureReagent(Reagent):
    molar_mass: Optional[pint.Quantity] = None


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
