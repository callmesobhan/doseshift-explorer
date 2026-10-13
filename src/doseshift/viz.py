"""Matplotlib figures. Every function returns a Figure; none touches the backend.

Colour follows the dataviz method: categorical hues are assigned to *entities* (a structure keeps
its colour whatever else is plotted) in a fixed, validated slot order and never cycled; dose
magnitude uses a one-hue light-to-dark ramp; text always wears text tokens, never series colour.
Light/dark steps and the slot order were checked with the dataviz skill's palette validator
(adjacent-pair CVD dE >= 8, normal-vision dE >= 15).
"""

from __future__ import annotations

import textwrap
import warnings

import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib import patheffects
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .evaluate import Evaluation
from .geometry import Grid, sample_volume
from .metrics import Constraint, metric_unit
from .montecarlo import MCResult
from .state import PlanState

THEMES = {
    "light": dict(
        surface="#fcfcfb",
        ink="#0b0b0b",
        ink2="#52514e",
        muted="#8a8984",
        grid="#e6e5e1",
        series=[
            "#2a78d6",
            "#eb6834",
            "#1baf7a",
            "#eda100",
            "#e87ba4",
            "#008300",
            "#4a3aa7",
            "#e34948",
        ],
    ),
    "dark": dict(
        surface="#1a1a19",
        ink="#ffffff",
        ink2="#c3c2b7",
        muted="#8e8d85",
        grid="#2f2f2d",
        series=[
            "#3987e5",
            "#d95926",
            "#199e70",
            "#c98500",
            "#d55181",
            "#008300",
            "#9085e9",
            "#e66767",
        ],
    ),
}
SEQUENTIAL_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DOSE_CMAP = mcolors.LinearSegmentedColormap.from_list("dose_blue", SEQUENTIAL_BLUE)
DPI = 160


def entity_colors(names: list[str], theme: str = "light") -> dict[str, str]:
    """A fixed colour slot per entity, in the order given (pass the plan's full structure list so a
    structure keeps its colour in every figure). Beyond the 8 validated slots the rest share a
    neutral grey and are told apart by the legend."""
    series = THEMES[theme]["series"]
    if len(names) > len(series):
        warnings.warn(
            f"{len(names)} structures exceed the {len(series)} validated colours; "
            f"{', '.join(names[len(series) :])} share a neutral grey",
            stacklevel=2,
        )
    return {
        n: series[i] if i < len(series) else THEMES[theme]["muted"] for i, n in enumerate(names)
    }


def _fig(width: float, height: float, theme: str, **kw):
    t = THEMES[theme]
    fig = Figure(figsize=(width, height), dpi=DPI, facecolor=t["surface"], **kw)
    return fig, t


def _style(ax, t, grid_axis="y") -> None:
    ax.set_facecolor(t["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
        ax.spines[side].set_linewidth(1)
    ax.tick_params(colors=t["ink2"], labelsize=8, length=3, width=1, color=t["grid"])
    ax.xaxis.label.set_color(t["ink2"])
    ax.yaxis.label.set_color(t["ink2"])
    ax.xaxis.label.set_size(9)
    ax.yaxis.label.set_size(9)
    ax.set_axisbelow(True)
    if grid_axis:
        ax.grid(axis=grid_axis, color=t["grid"], linewidth=1, linestyle="-")


def _title(fig, t, title, subtitle=None) -> None:
    """Title and subtitle placed in inches from the top, so short figures do not collide."""
    h = fig.get_figheight()
    if subtitle:  # never run off the figure: shorten to what fits
        subtitle = textwrap.shorten(
            subtitle, width=max(20, int(fig.get_figwidth() * 15)), placeholder=" ..."
        )
    fig.text(
        0.01,
        1 - 0.07 / h,
        title,
        color=t["ink"],
        fontsize=11,
        fontweight="bold",
        va="top",
        ha="left",
    )
    if subtitle:
        fig.text(0.01, 1 - 0.36 / h, subtitle, color=t["ink2"], fontsize=8.5, va="top", ha="left")


def _legend(ax, handles, t, **kw):
    leg = ax.legend(handles=handles, frameon=False, fontsize=8, **kw)
    for text in leg.get_texts():
        text.set_color(t["ink"])
    return leg


def _label(ax, x, y, s, t, **kw):
    return ax.text(x, y, s, color=t["ink2"], fontsize=8, **kw)


_LINESTYLES = ["-", (0, (5, 2.5)), (0, (1.2, 2.2)), (0, (5, 2, 1.2, 2))]


def plot_dvh_overlay(
    evals: dict[str, Evaluation],
    rx_gy: float,
    structures: list[str] | None = None,
    title: str = "DVH under perturbation",
    subtitle: str | None = None,
    theme: str = "light",
    rx_label: str = "Rx",
    color_names: list[str] | None = None,
) -> Figure:
    """Cumulative DVHs; colour = structure, line style = scenario (first is the baseline).

    ``color_names`` fixes the colour slots (pass the plan's full structure list to keep colours
    identical across figures that show different subsets)."""
    labels = list(evals)
    first = evals[labels[0]]
    structures = structures or list(first.dvhs)
    colors = entity_colors(color_names or list(first.dvhs), theme)
    fig, t = _fig(7.2, 4.4, theme)
    ax = fig.add_axes([0.085, 0.13, 0.89, 0.72])
    _style(ax, t)
    top = max(1.25 * rx_gy, 1.05 * max(e.dvhs[s].dmax for e in evals.values() for s in structures))
    axis = np.linspace(0.0, top, 600)
    for li, label in enumerate(labels):
        for s in structures:
            ax.plot(
                axis,
                100 * evals[label].dvhs[s].curve(axis),
                color=colors[s],
                linewidth=2,
                linestyle=_LINESTYLES[li % len(_LINESTYLES)],
                solid_capstyle="round",
            )
    ax.axvline(rx_gy, color=t["muted"], linewidth=1)
    _label(ax, rx_gy + 0.6, 104.5, f"{rx_label} {rx_gy:g} Gy", t, va="bottom")
    ax.set_xlim(0, axis[-1])
    ax.set_ylim(0, 104)
    ax.set_yticks(range(0, 101, 20))
    ax.set_xlabel("Dose (Gy)")
    ax.set_ylabel("Volume (%)")
    handles = [Line2D([], [], color=colors[s], linewidth=2, label=s) for s in structures]
    if len(labels) > 1:
        handles += [
            Line2D(
                [],
                [],
                color=t["ink2"],
                linewidth=2,
                label=lab,
                linestyle=_LINESTYLES[i % len(_LINESTYLES)],
            )
            for i, lab in enumerate(labels)
        ]
    _legend(ax, handles, t, loc="best")
    _title(fig, t, title, subtitle)
    return fig


def plot_mc_bands(
    result: MCResult,
    structures: list[str] | None = None,
    q: tuple[float, float] = (5, 95),
    theme: str = "light",
    rx_label: str = "Rx",
    color_names: list[str] | None = None,
) -> Figure:
    """Nominal DVH (line) inside the pointwise q-percentile envelope of the course DVHs (wash)."""
    structures = structures or list(result.dvh_curves)
    colors = entity_colors(color_names or list(result.dvh_curves), theme)
    fig, t = _fig(7.2, 4.4, theme)
    ax = fig.add_axes([0.085, 0.13, 0.89, 0.72])
    _style(ax, t)
    x = result.dvh_axis
    for s in structures:
        lo, hi = np.percentile(result.dvh_curves[s], q, axis=0) * 100
        ax.fill_between(x, lo, hi, color=colors[s], alpha=0.14, linewidth=0)
        ax.plot(
            x,
            100 * result.nominal.dvhs[s].curve(x),
            color=colors[s],
            linewidth=2,
            solid_capstyle="round",
        )
    rx = result.prescription_gy
    ax.axvline(rx, color=t["muted"], linewidth=1)
    _label(ax, rx + 0.6, 104.5, f"{rx_label} {rx:g} Gy", t, va="bottom")
    ax.set_xlim(0, x[-1])
    ax.set_ylim(0, 104)
    ax.set_yticks(range(0, 101, 20))
    ax.set_xlabel("Dose (Gy)")
    ax.set_ylabel("Volume (%)")
    handles = [Line2D([], [], color=colors[s], linewidth=2, label=s) for s in structures]
    handles.append(
        Patch(
            facecolor=t["muted"],
            alpha=0.3,
            label=f"pointwise {q[0]:g}-{q[1]:g}th percentile of course DVHs",
        )
    )
    _legend(ax, handles, t, loc="best")
    m = result.model
    _title(
        fig,
        t,
        "Delivered DVH under setup error",
        f"{result.n_courses} courses, {m.n_fractions} fractions ({result.method}); "
        f"Sigma={m.sys_sigma_mm[0]:g}, sigma={m.rand_sigma_mm[0]:g} mm; line = planned",
    )
    return fig


def plot_metric_distribution(
    result: MCResult,
    key: str,
    constraint: Constraint | None = None,
    theme: str = "light",
    color_names: list[str] | None = None,
) -> Figure:
    """Histogram of one metric across courses, with the planned value and an optional threshold."""
    vals = result.metrics[key].to_numpy()
    fig, t = _fig(6.4, 3.8, theme)
    ax = fig.add_axes([0.1, 0.15, 0.87, 0.68])
    _style(ax, t)
    structure, metric = key.split(":", 1)
    colors = entity_colors(color_names or list(result.dvh_curves), theme)
    unit = metric_unit(metric)
    ax.hist(vals, bins=30, color=colors[structure], edgecolor=t["surface"], linewidth=1)
    top = ax.get_ylim()[1]
    nominal = result.nominal.flat()[key]
    ax.axvline(nominal, color=t["ink"], linewidth=1.5)
    _label(ax, nominal, top * 1.01, "planned ", t, va="bottom", ha="right")
    sub = f"{result.n_courses} simulated courses"
    if constraint is not None:
        thr = constraint.threshold(result.prescription_gy)
        ax.axvline(thr, color=t["muted"], linewidth=1.5)
        _label(ax, thr, top * 1.01, f"threshold {thr:.1f} {unit} ", t, va="bottom", ha="right")
        p = result.probability(constraint)
        sub = f"P({constraint}) = {p:.0%} over {result.n_courses} courses"
    ax.set_xlabel(f"{key} ({unit})")
    ax.set_ylabel("Courses")
    _title(fig, t, f"Distribution of {key}", sub)
    return fig


def _signed(sweeps: dict[str, pd.DataFrame], axis: str, key: str) -> tuple[np.ndarray, np.ndarray]:
    neg, pos = sweeps[f"-{axis}"], sweeps[f"+{axis}"]
    x = np.concatenate([neg["shift_mm"].to_numpy()[:0:-1], pos["shift_mm"].to_numpy()])
    y = np.concatenate([neg[key].to_numpy()[:0:-1], pos[key].to_numpy()])
    return x, y


def plot_sweeps(
    sweeps: dict[str, pd.DataFrame],
    keys: list[str],
    thresholds: dict[str, float] | None = None,
    theme: str = "light",
) -> Figure:
    """Small multiples: one panel per metric, one line per patient-shift axis."""
    thresholds = thresholds or {}
    axes_names = {"x": "x (+ = patient left)", "y": "y (+ = posterior)", "z": "z (+ = superior)"}
    fig, t = _fig(3.0 * len(keys) + 0.6, 3.9, theme)
    # slots 5-7: kept apart from the slots structures use in the other figures
    colors = dict(zip("xyz", THEMES[theme]["series"][4:7], strict=True))
    panels = np.atleast_1d(fig.subplots(1, len(keys)))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.78, bottom=0.17, wspace=0.32)
    for ax, key in zip(panels, keys, strict=True):
        _style(ax, t)
        for a in "xyz":
            x, y = _signed(sweeps, a, key)
            ax.plot(x, y, color=colors[a], linewidth=2, solid_capstyle="round")
        if key in thresholds:
            ax.axhline(thresholds[key], color=t["muted"], linewidth=1)
            _label(ax, ax.get_xlim()[0], thresholds[key], " limit", t, va="bottom")
        ax.set_title(key, color=t["ink"], fontsize=9, loc="left")
        ax.set_xlabel("Patient shift (mm)")
        ax.set_ylabel("Gy" if key.split(":")[1].startswith(("D", "gEUD")) else "")
    handles = [Line2D([], [], color=colors[a], linewidth=2, label=axes_names[a]) for a in "xyz"]
    fig.legend(
        handles=handles,
        loc="upper right",
        ncol=3,
        frameon=False,
        fontsize=8,
        labelcolor=t["ink"],
        bbox_to_anchor=(0.99, 0.985),
    )
    _title(
        fig,
        t,
        "Sensitivity to rigid patient shift",
        "Dose is fixed; the patient moves. Positive x = patient left, y = posterior, z = superior",
    )
    return fig


def plot_equivalence(
    eq: pd.DataFrame,
    metrics: list[str],
    max_mm: float = 12.0,
    theme: str = "light",
    title: str = "How big a rigid shift is each grid assumption worth?",
    subtitle: str = "Smallest axis-aligned patient shift with the same effect (lower bound)",
    na_label: str = "n/a: no rigid shift reproduces this",
    color_names: list[str] | None = None,
) -> Figure:
    """Grid-assumption audit: each convention choice expressed as an equivalent rigid shift.

    Bars are lower bounds (smallest axis-aligned shift of the right sign). For dose-grid offsets
    a black tick marks the *actual* displacement of the dose, so the gap between bar and tick is
    visible rather than hidden.
    """
    scen = list(dict.fromkeys(eq["scenario"]))
    pivot = eq.pivot(index="scenario", columns="metric", values="equivalent_shift_mm").loc[scen]
    order = pivot[metrics].replace(np.inf, max_mm * 1.5).max(axis=1).sort_values().index.tolist()
    height = 1.5 + 0.27 * len(order)
    fig, t = _fig(4.2 + 3.6 * len(metrics), height, theme)
    structures_here = list(dict.fromkeys(k.split(":", 1)[0] for k in metrics))
    slot = entity_colors(color_names or structures_here, theme)
    panels = np.atleast_1d(fig.subplots(1, len(metrics), sharey=True))
    fig.subplots_adjust(
        left=0.36 / len(metrics) ** 0.5 * (1 if len(metrics) == 1 else 0.8),
        right=0.98,
        top=1 - 1.0 / height,
        bottom=0.55 / height,
        wspace=0.06,
    )
    ys = np.arange(len(order))
    for ax, key in zip(panels, metrics, strict=True):
        _style(ax, t, grid_axis="x")
        color = slot[key.split(":", 1)[0]]
        sub = eq[eq["metric"] == key].set_index("scenario")
        for y, name in zip(ys, order, strict=True):
            row = sub.loc[name]
            mm = row["equivalent_shift_mm"]
            if np.isnan(mm):
                _label(ax, 0.1, y, na_label, t, va="center", ha="left")
                continue
            if mm == 0:
                _label(ax, 0.1, y, "negligible", t, va="center", ha="left")
                continue
            length = max_mm if np.isinf(mm) else mm
            ax.barh(y, length, height=0.62, color=color, linewidth=0)
            disp = row.get("displacement_mm", np.nan)
            if np.isfinite(disp):
                ax.plot([disp, disp], [y - 0.42, y + 0.42], color=t["ink"], linewidth=1.8)
            text = (
                f">{max_mm:g} mm"
                if np.isinf(mm)
                else "<0.1 mm"
                if mm < 0.05
                else f"{mm:.1f} mm ({row['direction']})"
            )
            reach = max(length, disp) if np.isfinite(disp) else length
            _label(ax, reach + 0.25, y, text, t, va="center", ha="left")
        ax.set_yticks(ys)
        ax.set_yticklabels(order, fontsize=8, color=t["ink"])
        ax.set_ylim(len(order) - 0.4, -0.6)
        finite = [
            v
            for v in list(sub["equivalent_shift_mm"]) + list(sub.get("displacement_mm", []))
            if np.isfinite(v)
        ]
        any_inf = bool(np.isinf(sub["equivalent_shift_mm"]).any())
        ax.set_xlim(0, max_mm * 1.5 if any_inf or not finite else max(2.0, 1.8 * max(finite)))
        ax.set_xlabel("Equivalent rigid shift (mm)")
        ax.set_title(key, color=t["ink"], fontsize=9, loc="left")
        ax.tick_params(axis="y", length=0)
    if eq.get("displacement_mm") is not None and eq["displacement_mm"].notna().any():
        panels[-1].legend(
            handles=[
                Line2D([], [], color=t["ink"], linewidth=1.8, label="actual dose displacement")
            ],
            loc="lower right",
            frameon=False,
            fontsize=8,
            labelcolor=t["ink"],
        )
    _title(fig, t, title, subtitle)
    return fig


def _plane(grid: Grid, array, plane: str, center, half_u, half_v, step, order=1, cval=0.0):
    """Sample ``array`` on an axial ('ax': u=x, v=y) or sagittal ('sag': u=y, v=z) plane."""
    u = np.arange(-half_u, half_u + step, step)
    v = np.arange(-half_v, half_v + step, step)
    uu, vv = np.meshgrid(u, v)
    if plane == "ax":
        pts = np.stack([center[0] + uu, center[1] + vv, np.full_like(uu, center[2])], -1)
    else:
        pts = np.stack([np.full_like(uu, center[0]), center[1] + uu, center[2] + vv], -1)
    out = sample_volume(grid, array, pts.reshape(-1, 3), order, cval).reshape(uu.shape)
    return u, v, out


def _image_window(img: np.ndarray) -> tuple[float, float]:
    """CT (air at about -1000 HU) gets a soft-tissue/bone window; anything else (MR) percentiles."""
    if np.nanmin(img) < -500:
        return -250.0, 500.0
    lo, hi = np.percentile(img, [1, 99.5])
    return float(lo), float(hi)


def plot_slices(
    state: PlanState,
    ct: tuple[Grid, np.ndarray] | None = None,
    structures: list[str] | None = None,
    theme: str = "light",
    title: str = "Plan geometry",
) -> Figure:
    """Axial and sagittal slices through the isocentre: CT, dose wash, structure outlines."""
    structures = structures or state.names
    colors = entity_colors(state.names, theme)
    iso = state.isocenter()
    pts = np.concatenate([state.clouds[n].points for n in structures])
    reach = np.maximum(np.abs(pts.min(axis=0) - iso), np.abs(pts.max(axis=0) - iso)) + 25.0
    half = {"ax": (reach[0], reach[1]), "sag": (reach[1], reach[2])}
    step = max(0.25, float(reach.max()) / 220.0)
    dmax = float(state.dose.data.max())
    window = _image_window(ct[1]) if ct is not None else None
    fig, t = _fig(9.0, 4.6, theme)
    halo = [patheffects.withStroke(linewidth=3.6, foreground=t["surface"])]
    for n, plane in enumerate(("ax", "sag")):
        ax = fig.add_axes([0.04 + 0.5 * n, 0.04, 0.44, 0.80])
        ax.set_facecolor(t["surface"])
        ax.set_aspect("equal")
        ax.set_anchor("N")
        u, v, dose = _plane(state.dose.grid, state.dose.data, plane, iso, *half[plane], step)
        ext = [u[0], u[-1], v[-1], v[0]] if plane == "ax" else [u[0], u[-1], v[0], v[-1]]
        if ct is not None:
            _, _, img = _plane(
                ct[0], ct[1], plane, iso, *half[plane], step, cval=float(ct[1].min())
            )
            ax.imshow(
                img,
                extent=ext,
                origin="upper" if plane == "ax" else "lower",
                cmap="gray",
                vmin=window[0],
                vmax=window[1],
                interpolation="bilinear",
            )
        shown = np.ma.masked_less(dose, 0.1 * dmax)
        im = ax.imshow(
            shown,
            extent=ext,
            origin="upper" if plane == "ax" else "lower",
            cmap=DOSE_CMAP,
            vmin=0,
            vmax=dmax,
            alpha=0.6,
            interpolation="bilinear",
        )
        for s in structures:
            m = state.masks[s]
            _, _, mm = _plane(m.grid, m.mask.astype(float), plane, iso, *half[plane], step, order=0)
            cs = ax.contour(
                u, v, mm, levels=[0.5], colors=[colors[s]], linewidths=1.8, linestyles="solid"
            )
            cs.set_path_effects(halo)  # 2px surface ring keeps outlines legible over the wash
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_title(
            "Axial (left on the right of the image)"
            if plane == "ax"
            else "Sagittal (anterior on the left)",
            color=t["ink2"],
            fontsize=9,
            loc="left",
        )
    cax = fig.add_axes([0.30, 0.94, 0.20, 0.018])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Dose (Gy)", color=t["ink2"], fontsize=8, labelpad=-30)
    cb.ax.tick_params(colors=t["ink2"], labelsize=7, length=2)
    cb.outline.set_visible(False)
    handles = [Line2D([], [], color=colors[s], linewidth=2, label=s) for s in structures]
    fig.legend(
        handles=handles,
        loc="upper right",
        ncol=len(structures),
        frameon=False,
        fontsize=8,
        labelcolor=t["ink"],
        bbox_to_anchor=(0.99, 0.995),
    )
    _title(fig, t, title)
    return fig


def save(fig: Figure, path) -> None:
    fig.savefig(path, dpi=DPI, facecolor=fig.get_facecolor())


def plot_margin_scan(
    df: pd.DataFrame,
    recipe_mm: float | None = None,
    target: float = 0.9,
    theme: str = "light",
    label: str = "P(constraint met)",
    color: str | None = None,
) -> Figure:
    """Probability of meeting a constraint versus CTV->PTV margin (``color``: the structure's)."""
    fig, t = _fig(6.2, 3.8, theme)
    ax = fig.add_axes([0.1, 0.16, 0.86, 0.66])
    _style(ax, t)
    color = color or THEMES[theme]["series"][0]
    ax.axhline(target, color=t["muted"], linewidth=1)
    _label(ax, df["margin_mm"].min(), target + 0.015, f" target {target:.0%}", t, va="bottom")
    if "se" in df:
        ax.fill_between(
            df["margin_mm"],
            np.clip(df["probability"] - 2 * df["se"], 0, 1),
            np.clip(df["probability"] + 2 * df["se"], 0, 1),
            color=color,
            alpha=0.14,
            linewidth=0,
        )
    ax.plot(
        df["margin_mm"],
        df["probability"],
        color=color,
        linewidth=2,
        marker="o",
        markersize=5,
        markeredgecolor=t["surface"],
        markeredgewidth=2,
    )
    if recipe_mm is not None:
        ax.axvline(recipe_mm, color=t["ink2"], linewidth=1)
        _label(ax, recipe_mm + 0.1, 0.04, f"van Herk recipe {recipe_mm:.1f} mm", t, va="bottom")
    ax.set_ylim(-0.04, 1.05)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("CTV-PTV margin (mm)")
    ax.set_ylabel(label)
    _title(
        fig,
        t,
        "Margin needed for robust coverage",
        "Monte Carlo over simulated courses (band = 2 standard errors); synthetic phantom only",
    )
    return fig


def plot_resolution_scan(
    df: pd.DataFrame, theme: str = "light", color_names: list[str] | None = None
) -> Figure:
    """Cost of a half-voxel convention error against the displacement it amounts to (one plan).

    Each point is the same plan with its dose coarsened by 1x, 2x, 3x ... so only voxel size varies.
    A straight line through the origin would mean the cost is linear in voxel size.
    """
    from .experiments import loglog_slope

    keys = list(dict.fromkeys(df["metric"]))
    structures = list(dict.fromkeys(k.split(":")[0] for k in keys))
    colors = entity_colors(color_names or structures, theme)
    styles = ["-", (0, (5, 2.5)), (0, (1.2, 2.2))]
    fig, t = _fig(6.6, 4.2, theme)
    ax = fig.add_axes([0.1, 0.15, 0.86, 0.66])
    _style(ax, t)
    seen: dict[str, int] = {}
    handles = []
    for key in keys:
        d = df[df["metric"] == key].sort_values("displacement_mm")
        s_name = key.split(":")[0]
        style = styles[seen.get(s_name, 0) % len(styles)]
        seen[s_name] = seen.get(s_name, 0) + 1
        ax.plot(
            d["displacement_mm"],
            d["cost"],
            color=colors[s_name],
            linewidth=2,
            linestyle=style,
            marker="o",
            markersize=5,
            markeredgecolor=t["surface"],
            markeredgewidth=2,
        )
        handles.append(
            Line2D(
                [],
                [],
                color=colors[s_name],
                linewidth=2,
                linestyle=style,
                label=f"{key} (log-log slope {loglog_slope(df, key):.1f})",
            )
        )
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Dose displacement of a half-voxel offset, |A v| (mm)")
    ax.set_ylabel("Mean |change| of the metric (Gy)")
    _legend(ax, handles, t, loc="upper left")
    _title(
        fig,
        t,
        "Does the cost of a half-voxel error scale with voxel size?",
        "One plan, dose coarsened 1x to 4x; slope 1 = linear in voxel size",
    )
    return fig
