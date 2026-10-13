"""Exact (un-binned) cumulative DVH from weighted dose samples."""

from __future__ import annotations

import numpy as np


class DVH:
    """Cumulative DVH of ``doses`` (Gy) weighted by ``volumes`` (mm^3).

    All metrics are evaluated on the sorted samples rather than on a histogram, so they carry no
    bin-width error and are comparable across perturbed states.
    """

    def __init__(self, doses: np.ndarray, volumes: np.ndarray):
        doses = np.asarray(doses, dtype=float)
        volumes = np.asarray(volumes, dtype=float)
        if doses.shape != volumes.shape or doses.ndim != 1 or doses.size == 0:
            raise ValueError("doses and volumes must be equal-length, non-empty 1D arrays")
        if not (np.isfinite(doses).all() and np.isfinite(volumes).all()):
            raise ValueError(
                "doses and volumes must be finite (NaN/inf dose would corrupt every metric)"
            )
        order = np.argsort(doses, kind="stable")
        self._d = doses[order]
        w = volumes[order]
        self._suffix = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])  # volume with index >= i
        self.total_volume_mm3 = float(self._suffix[0])
        self._mean = float((self._d * w).sum() / self.total_volume_mm3)
        self._w = w

    @property
    def total_volume_cc(self) -> float:
        return self.total_volume_mm3 / 1000.0

    @property
    def dmin(self) -> float:
        return float(self._d[0])

    @property
    def dmax(self) -> float:
        return float(self._d[-1])

    @property
    def dmean(self) -> float:
        return self._mean

    def volume_fraction_at_dose(self, dose_gy) -> np.ndarray:
        """Fraction of volume receiving >= dose (vectorised)."""
        i = np.searchsorted(self._d, np.asarray(dose_gy, dtype=float), side="left")
        return self._suffix[i] / self.total_volume_mm3

    def dose_at_volume_fraction(self, fraction: float) -> float:
        """Highest dose D such that at least ``fraction`` of the volume receives >= D."""
        fraction = float(np.clip(fraction, 0.0, 1.0))
        target = fraction * self.total_volume_mm3 * (1.0 - 1e-12)
        i = np.searchsorted(-self._suffix, -target, side="right") - 1
        return float(self._d[min(i, self._d.size - 1)])

    def dose_at_volume_cc(self, cc: float) -> float:
        """Dose received by at least ``cc`` of absolute volume; NaN if the structure is smaller."""
        if cc > self.total_volume_cc:
            return float("nan")
        return self.dose_at_volume_fraction(cc / self.total_volume_cc)

    def volume_cc_at_dose(self, dose_gy: float) -> float:
        return float(self.volume_fraction_at_dose(dose_gy)) * self.total_volume_cc

    def generalized_eud(self, a: float) -> float:
        """gEUD = (sum v_i d_i^a)^(1/a); a>1 serial-like OARs, a<0 tumours (cold-spot weighted)."""
        v = self._w / self.total_volume_mm3
        d = np.maximum(self._d, 1e-9)
        if a == 0:  # the limit a -> 0 is the geometric mean
            return float(np.exp((v * np.log(d)).sum()))
        return float((v * d**a).sum() ** (1.0 / a))

    def curve(self, dose_axis: np.ndarray) -> np.ndarray:
        """Cumulative volume fraction on a shared dose axis (for plotting / MC envelopes)."""
        return self.volume_fraction_at_dose(dose_axis)
