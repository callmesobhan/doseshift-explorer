"""The Perturbation contract: a pure, composable ``PlanState -> PlanState`` transform."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

from ..state import PlanState


class Perturbation(ABC):
    """``kind`` separates what really happens to the patient from what the pipeline assumes.

    * ``physical``  - the patient is somewhere other than planned (setup error, rotation,
      deformation). Real clinical risk.
    * ``numerical`` - the patient is where planned but the dose/structure grids are interpreted
      differently (voxel-centre vs corner, resolution, interpolation, mask sampling). Pure
      software/convention risk.
    """

    kind: ClassVar[str]

    @abstractmethod
    def apply(self, state: PlanState) -> PlanState: ...

    @abstractmethod
    def describe(self) -> str: ...

    def __call__(self, state: PlanState) -> PlanState:
        return self.apply(state)


class Compose(Perturbation):
    """Apply perturbations left to right. Numerical sampling steps must precede motion."""

    kind = "mixed"

    def __init__(self, *steps: Perturbation):
        self.steps = steps

    def apply(self, state: PlanState) -> PlanState:
        for step in self.steps:
            state = step.apply(state)
        return state

    def describe(self) -> str:
        return " + ".join(s.describe() for s in self.steps) or "identity"
