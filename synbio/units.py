from typing import Union

import pint

__all__ = ["unit_registry", "QuantityType"]

unit_registry = pint.UnitRegistry()
QuantityType = Union[int, float, pint.Quantity]
