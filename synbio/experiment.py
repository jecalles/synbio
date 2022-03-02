from datetime import date
from typing import Any, List

from synbio.interfaces import HashableMixin

__all__ = [
    # Dataclasses
    "Condition", "Experiment", "Data"
]


class Condition(HashableMixin):
    def _comparables(self) -> List[str]:
        return list(self.variables.keys())

    def __repr__(self) -> str:
        class_name = str(self.__class__.__name__)
        var_string = " ".join(f"{k}={v}" for k, v in self.variables.items())
        return f"{class_name}({var_string})"

    @property
    def variables(self) -> dict:
        return {
            k: v for k, v in self.__dict__.items()
            if k[0] != "_"  # at the moment, this is all the filtering we're
            # doing
        }

    @variables.setter
    def variables(self, new_dict: dict):
        self.__dict__.update(new_dict)


class Data: pass


class Experiment:
    def __init__(
            self, name: str = None,
            conditions: List[Condition] = None,
            data: Any = None,
            meta: dict = None
    ):
        today = date.today()

        if name is None:
            name = str(today)

        if meta is None:
            meta = {}
        meta["date"] = today

        self.name = name
        self.conditions = conditions
        self.data = data
        self.meta = meta

    @property
    def variables(self):
        return {
            var for condition in self.conditions
            for var in condition.variables
        }
