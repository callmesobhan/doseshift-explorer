"""Plan-quality metrics, addressed by short strings so they can live in YAML configs.

Grammar (``x`` is a number):

======================  ==========================================================  =====
spec                    meaning                                                     unit
======================  ==========================================================  =====
``Dmean`` ``Dmax``      mean / maximum sampled dose                                 Gy
``Dmin``
``D{x}%``               dose received by at least x % of the volume                 Gy
``D{x}cc``              dose received by at least x cc                              Gy
``V{x}%``               % of volume receiving >= x % of the prescription            %
``V{x}Gy``              % of volume receiving >= x Gy                               %
``V{x}Gy_cc``           volume (cc) receiving >= x Gy                               cc
``gEUD(a={x})``         generalised equivalent uniform dose                         Gy
``HI``                  ICRU 83 homogeneity index (D2% - D98%) / D50%               -
``CI``                  Paddick conformity index TV_PIV^2 / (TV * PIV)              -
``GI``                  gradient index V50%Rx / V100%Rx (dose-cloud volumes)       -
======================  ==========================================================  =====

``CI`` and ``GI`` use the isodose volumes of the *dose grid* (PIV, V50), which do not move with
the patient; only ``TV_PIV`` does. That is why grid assumptions and setup shifts affect them
differently.
"""

from __future__ import annotations

import operator
import re
from dataclasses import dataclass

import numpy as np

from .dvh import DVH

_NUM = r"(?P<x>\d+(?:\.\d+)?)"
_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(rf"^D{_NUM}%$"), "dose_at_pct"),
    (re.compile(rf"^D{_NUM}cc$"), "dose_at_cc"),
    (re.compile(rf"^V{_NUM}%$"), "vol_pct_rx"),
    (re.compile(rf"^V{_NUM}Gy_cc$"), "vol_cc_gy"),
    (re.compile(rf"^V{_NUM}Gy$"), "vol_pct_gy"),
    (re.compile(r"^gEUD\(a=(?P<x>-?\d+(?:\.\d+)?)\)$"), "geud"),
]
_SIMPLE = {"Dmean", "Dmax", "Dmin", "HI", "CI", "GI"}


@dataclass(frozen=True)
class MetricContext:
    """Everything a metric may need for one structure in one perturbed state."""

    dvh: DVH
    doses: np.ndarray
    volumes: np.ndarray
    rx_gy: float
    piv_mm3: float  # volume of the dose grid receiving >= Rx
    v50_mm3: float  # volume of the dose grid receiving >= 50 % Rx


def metric_unit(spec: str) -> str:
    if spec in ("Dmean", "Dmax", "Dmin") or spec.startswith(("gEUD",)):
        return "Gy"
    if spec in ("HI", "CI", "GI"):
        return "-"
    if spec.startswith("D"):
        return "Gy"
    if spec.endswith("_cc"):
        return "cc"
    return "%"


def compute_metric(spec: str, ctx: MetricContext) -> float:
    d = ctx.dvh
    if spec in _SIMPLE:
        return _simple(spec, ctx)
    for pattern, kind in _PATTERNS:
        m = pattern.match(spec)
        if not m:
            continue
        x = float(m.group("x"))
        if kind == "dose_at_pct":
            return d.dose_at_volume_fraction(x / 100.0)
        if kind == "dose_at_cc":
            return d.dose_at_volume_cc(x)
        if kind == "vol_pct_rx":
            return 100.0 * float(d.volume_fraction_at_dose(x / 100.0 * ctx.rx_gy))
        if kind == "vol_pct_gy":
            return 100.0 * float(d.volume_fraction_at_dose(x))
        if kind == "vol_cc_gy":
            return d.volume_cc_at_dose(x)
        if kind == "geud":
            return d.generalized_eud(x)
    raise ValueError(f"unknown metric spec {spec!r}")


def _simple(spec: str, ctx: MetricContext) -> float:
    d = ctx.dvh
    if spec == "Dmean":
        return d.dmean
    if spec == "Dmax":
        return d.dmax
    if spec == "Dmin":
        return d.dmin
    if spec == "HI":
        return (d.dose_at_volume_fraction(0.02) - d.dose_at_volume_fraction(0.98)) / max(
            d.dose_at_volume_fraction(0.50), 1e-9
        )
    tv = d.total_volume_mm3
    tv_piv = float(ctx.volumes[ctx.doses >= ctx.rx_gy].sum())
    if spec == "CI":
        return tv_piv**2 / (tv * ctx.piv_mm3) if ctx.piv_mm3 > 0 else float("nan")
    return ctx.v50_mm3 / ctx.piv_mm3 if ctx.piv_mm3 > 0 else float("nan")  # GI


DEFAULT_METRICS: dict[str, list[str]] = {
    "ptv": ["D98%", "D95%", "D50%", "D2%", "Dmean", "V95%", "HI", "CI", "GI"],
    "ctv": ["D98%", "D95%", "Dmean", "V95%", "V100%"],
    "oar": ["Dmean", "D2%", "V20%", "V50%"],
}

_OPS = {">=": operator.ge, "<=": operator.le, ">": operator.gt, "<": operator.lt}


@dataclass(frozen=True)
class Constraint:
    """A pass/fail criterion, e.g. ``CTV D98% >= 95 % of Rx``."""

    structure: str
    metric: str
    op: str
    value: float
    pct_of_rx: bool = False  # interpret ``value`` as % of prescription (Gy-valued metrics)

    def threshold(self, rx_gy: float) -> float:
        return self.value / 100.0 * rx_gy if self.pct_of_rx else self.value

    def passes(self, values: np.ndarray, rx_gy: float) -> np.ndarray:
        return _OPS[self.op](np.asarray(values, dtype=float), self.threshold(rx_gy))

    @property
    def key(self) -> str:
        return f"{self.structure}:{self.metric}"

    def __str__(self) -> str:
        unit = "% Rx" if self.pct_of_rx else metric_unit(self.metric)
        return f"{self.structure} {self.metric} {self.op} {self.value:g} {unit}"
