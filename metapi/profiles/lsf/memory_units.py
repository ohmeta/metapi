import re
from enum import Enum
from typing import Union
from collections import namedtuple


class InvalidSuffix(Exception):
    pass


class InvalidPower(Exception):
    pass


class InvalidMemoryString(Exception):
    pass


Scale = namedtuple("Scale", ["power", "metric_suffix"])


SCALE_MAP = {
    "B": Scale(0, "B"),
    "K": Scale(1, "KB"),
    "M": Scale(2, "MB"),
    "G": Scale(3, "GB"),
    "T": Scale(4, "TB"),
    "P": Scale(5, "PB"),
    "E": Scale(6, "EB"),
    "Z": Scale(7, "ZB"),
}


class Unit(Enum):
    BYTES = SCALE_MAP["B"]
    KILO = SCALE_MAP["K"]
    MEGA = SCALE_MAP["M"]
    GIGA = SCALE_MAP["G"]
    TERA = SCALE_MAP["T"]
    PETA = SCALE_MAP["P"]
    EXA = SCALE_MAP["E"]
    ZETTA = SCALE_MAP["Z"]

    @staticmethod
    def from_suffix(suffix: str) -> "Unit":
        first_letter = suffix[0].upper()
        if first_letter not in SCALE_MAP:
            valid_suffixes = ",".join(
                scale.metric_suffix for scale in SCALE_MAP.values()
            )
            raise InvalidSuffix(
                "{suffix}. Valid suffixes are: {valid_suffixes}".format(
                    suffix=suffix, valid_suffixes=valid_suffixes
                )
            )
        return Unit(SCALE_MAP[first_letter])

    @staticmethod
    def from_power(power: int) -> "Unit":
        valid_powers = []
        for scale in SCALE_MAP.values():
            if scale.power == power:
                return Unit(scale)
            valid_powers.append(scale.power)

        raise InvalidPower(
            "{power}. Valid powers are: {valid}".format(
                power=power, valid=",".join(str(p) for p in valid_powers)
            )
        )

    @property
    def power(self) -> int:
        return self.value.power

    @property
    def suffix(self) -> str:
        return self.value.metric_suffix


Number = Union[int, float]


class Memory:
    def __init__(self, value: Number = 1, unit: Unit = Unit.BYTES):
        self.value = value
        self.unit = unit
        self._decimal_scaling_factor = 1000
        self._binary_scaling_factor = 1024

    def __eq__(self, other: "Memory") -> bool:
        return self.bytes() == other.bytes()

    def __repr__(self) -> str:
        val = (
            int(self.value)
            if isinstance(self.value, int) or self.value.is_integer()
            else self.value
        )
        return "{val}{sfx}".format(val=val, sfx=self.suffix)

    @property
    def power(self) -> int:
        return self.unit.power

    @property
    def suffix(self) -> str:
        return self.unit.suffix

    def _scaling_factor(self, decimal: bool = True) -> int:
        return self._decimal_scaling_factor if decimal else self._binary_scaling_factor

    def bytes(self, decimal_multiples: bool = True) -> float:
        scaling_factor = self._scaling_factor(decimal_multiples)
        return float(self.value * (scaling_factor ** self.power))

    def to(self, unit: Unit, decimal_multiples: bool = True) -> "Memory":
        scaling_factor = self._scaling_factor(decimal_multiples) ** unit.power
        size = self.bytes(decimal_multiples)
        size /= scaling_factor

        return Memory(size, unit)

    @staticmethod
    def from_str(s: str) -> "Memory":
        valid_suffixes = "".join(scale.metric_suffix for scale in SCALE_MAP.values())
        regex = re.compile(
            r"^(?P<size>[0-9]*\.?[0-9]+)\s*(?P<sfx>[{}]B?)?$".format(valid_suffixes),
            re.IGNORECASE,
        )
        match = regex.search(s)

        if not match:
            raise InvalidMemoryString("{s} is an invalid memory string.".format(s=s))

        size = float(match.group("size"))
        suffix = match.group("sfx") or "B"
        unit = Unit.from_suffix(suffix)

        return Memory(size, unit)
