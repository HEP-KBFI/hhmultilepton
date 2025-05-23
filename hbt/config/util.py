# coding: utf-8

"""
Config-related object definitions and utils.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from order import UniqueObject, TagMixin
from order.util import typed

from columnflow.types import Callable, Any, Sequence, Hashable, ClassVar


class TriggerLeg(object):
    """
    Container class storing information about trigger legs:

        - *pdg_id*: The id of the object that should have caused the trigger leg to fire.
        - *min_pt*: The minimum transverse momentum in GeV of the triggered object.
        - *trigger_bits*: Integer bit mask or masks describing whether the last filter of a trigger fired.
          See https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py.
          Per mask, any of the bits should match (*OR*). When multiple masks are configured, each of
          them should match (*AND*).

    For accepted types and conversions, see the *typed* setters implemented in this class.
    """

    def __init__(
        self,
        pdg_id: int | None = None,
        min_pt: float | int | None = None,
        trigger_bits: int | Sequence[int] | None = None,
    ):
        super().__init__()

        # instance members
        self._pdg_id = None
        self._min_pt = None
        self._trigger_bits = None

        # set initial values
        self.pdg_id = pdg_id
        self.min_pt = min_pt
        self.trigger_bits = trigger_bits

    def __repr__(self):
        return (
            f"<{self.__class__.__name__} "
            f"'pdg_id={self.pdg_id}, min_pt={self.min_pt}, trigger_bits={self.trigger_bits}' "
            f"at {hex(id(self))}>"
        )

    @typed
    def pdg_id(self, pdg_id: int | None) -> int | None:
        if pdg_id is None:
            return None

        if not isinstance(pdg_id, int):
            raise TypeError(f"invalid pdg_id: {pdg_id}")

        return pdg_id

    @typed
    def min_pt(self, min_pt: int | float | None) -> float | None:
        if min_pt is None:
            return None

        if isinstance(min_pt, int):
            min_pt = float(min_pt)
        if not isinstance(min_pt, float):
            raise TypeError(f"invalid min_pt: {min_pt}")

        return min_pt

    @typed
    def trigger_bits(
        self,
        trigger_bits: int | Sequence[int] | None,
    ) -> list[int] | None:
        if trigger_bits is None:
            return None

        # cast to list
        if isinstance(trigger_bits, tuple):
            trigger_bits = list(trigger_bits)
        elif not isinstance(trigger_bits, list):
            trigger_bits = [trigger_bits]

        # check bit types
        for bit in trigger_bits:
            if not isinstance(bit, int):
                raise TypeError(f"invalid trigger bit: {bit}")

        return trigger_bits


class Trigger(UniqueObject, TagMixin):
    """
    Container class storing information about triggers:

        - *name*: The path name of a trigger that should have fired.
        - *id*: A unique id of the trigger.
        - *run_range*: An inclusive range describing the runs where the trigger is to be applied
          (usually only defined by data). None in the tuple means no lower or upper boundary.
        - *legs*: A dictionary mapping arbitrary keys to :py:class:`TriggerLeg` objects containing
          additional information and constraints of particular trigger legs.
        - *applies_to_dataset*: A function that obtains an ``order.Dataset`` instance to decide
          whether the trigger applies to that dataset. Defaults to *True*.

    For accepted types and conversions, see the *typed* setters implemented in this class.

    In addition, a base class from *order* provides additional functionality via mixins:

        - *tags*: Trigger objects can be assigned *tags* that can be checked later on, e.g. to
          describe the type of the trigger ("single_mu", "cross", ...).
    """

    def __init__(
        self,
        name: str,
        id: int,
        run_range: tuple[int | None, int | None] | None = None,
        legs: dict[Hashable, TriggerLeg] | Sequence[TriggerLeg] | None = None,
        applies_to_dataset: Callable | bool | Any = True,
        tags: Any = None,
    ):
        UniqueObject.__init__(self, name, id)
        TagMixin.__init__(self, tags=tags)

        # force the id to be positive
        if self.id < 0:
            raise ValueError(f"trigger id must be positive, but found {self.id}")

        # instance members
        self._run_range = None
        self._leg = None
        self._applies_to_dataset = None

        # set initial values
        self.run_range = run_range
        self.legs = legs
        self.applies_to_dataset = applies_to_dataset

    def __repr__(self):
        return (
            f"<{self.__class__.__name__} 'name={self.name}, id={self.id}, nlegs={self.n_legs}' "
            f"at {hex(id(self))}>"
        )

    @typed
    def name(self, name: str) -> str:
        if not isinstance(name, str):
            raise TypeError(f"invalid name: {name}")
        if not name.startswith("HLT_"):
            raise ValueError(f"invalid name: {name}")

        return name

    @typed
    def run_range(
        self,
        run_range: Sequence[int] | None,
    ) -> tuple[int] | None:
        if run_range is None:
            return None

        # cast list to tuple
        if isinstance(run_range, list):
            run_range = tuple(run_range)

        # run_range must be a tuple with two integers
        if not isinstance(run_range, tuple):
            raise TypeError(f"invalid run_range: {run_range}")
        if len(run_range) != 2:
            raise ValueError(f"invalid run_range length: {run_range}")
        if not isinstance(run_range[0], int):
            raise ValueError(f"invalid run_range start: {run_range[0]}")
        if not isinstance(run_range[1], int):
            raise ValueError(f"invalid run_range end: {run_range[1]}")

        return run_range

    @typed
    def legs(
        self,
        legs: dict[Hashable, TriggerLeg] | Sequence[TriggerLeg] | TriggerLeg | None,
    ) -> dict[Hashable, TriggerLeg] | None:
        if legs is None:
            return None

        # cast to dict
        if isinstance(legs, TriggerLeg):
            legs = [legs]
        if isinstance(legs, (list, tuple)):
            legs = dict(enumerate(legs))

        # validate
        for key, leg in legs.items():
            if not isinstance(leg, TriggerLeg):
                raise TypeError(f"invalid trigger leg with key {key}: {leg}")

        return legs or None

    @typed
    def applies_to_dataset(self, func: Callable | bool | Any) -> Callable:
        if not callable(func):
            if func is not None:
                raise TypeError(f"invalid applies_to_dataset: {func}")
            decision = True if func is None else bool(func)
            func = lambda dataset_inst: decision

        return func

    @property
    def has_legs(self):
        return bool(self._legs)

    @property
    def n_legs(self):
        return len(self.legs) if self.has_legs else 0

    @property
    def hlt_field(self):
        # remove the first four "HLT_" characters
        return self.name[4:]


@dataclass
class TriggerBits:
    """
    Lightweight container wrapping trigger bits for different versions of NanoAOD.
    """

    v12: int | None = None
    v14: int | None = None

    supported_versions: ClassVar[set[int]] = {12, 14}

    def __post_init__(self) -> None:
        # versions might be strings such as "v12" that act as references
        cre = re.compile(r"^v(\d+)$")
        for v in self.supported_versions:
            attr = f"v{v}"
            # only check strings
            if not isinstance((val := getattr(self, attr)), str):
                continue
            # check format
            if not (m := cre.match(val)):
                raise ValueError(f"invalid reference {attr} -> {val}")
            # check if not circular
            ref_v = int(m.group(1))
            if v == ref_v:
                raise ValueError(f"reference to same version {attr} -> {val}")
            # check type of referred value
            ref_val = getattr(self, f"v{ref_v}")
            if not isinstance(ref_val, int) and v is not None:
                raise ValueError(f"wrong reference value in {attr} -> {val}")
            # set it
            setattr(self, attr, ref_val)

    def get(self, nano_version: int) -> int:
        if nano_version not in self.supported_versions:
            raise ValueError(f"nano_version {nano_version} not supported")
        return getattr(self, f"v{nano_version}")
