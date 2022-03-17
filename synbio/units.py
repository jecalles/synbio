from typing import TypeVar

import pint

__all__ = ["unit_registry", "QuantityType"]

unit_registry = pint.UnitRegistry()
QuantityType = TypeVar("QuantityType", int, float, pint.Quantity)
