"""DoseShift Explorer: interactive DVH robustness on the synthetic phantom.

    pip install -e ".[app]"
    streamlit run app/streamlit_app.py

Everything shown comes from the same library calls the CLI and tests use.
"""

import numpy as np
import pandas as pd
import streamlit as st

from doseshift import viz
from doseshift.evaluate import evaluate
from doseshift.experiments import directional_sweeps, equivalent_shift_mm
from doseshift.perturb import (
    Compose,
    DoseCoarsen,
    DoseGridOffset,
    InterpolationOrder,
    RigidShift,
    SmoothDeformation,
    StructureSampling,
)
from doseshift.synthetic import make_case

st.set_page_config(page_title="DoseShift Explorer", layout="wide")


@st.cache_resource
def load():
    case = make_case(seed=0)
    state = case.to_state()
    return case, state, directional_sweeps(state, max_mm=12.0, step_mm=0.5)


case, state, sweeps = load()
rx = state.prescription_gy

st.title("DoseShift Explorer")
st.caption(
    "Radiotherapy dose robustness under geometric uncertainty and grid-alignment assumptions"
)

with st.sidebar:
    st.header("Physical: patient moves")
    t = [st.slider(f"Shift {a} (mm)", -10.0, 10.0, 0.0, 0.5) for a in "xyz"]
    r = [st.slider(f"Rotate about {a} (deg)", -5.0, 5.0, 0.0, 0.5) for a in "xyz"]
    deform = st.slider("Smooth deformation RMS (mm)", 0.0, 5.0, 0.0, 0.5)
    st.header("Numerical: grid assumptions")
    off = [st.slider(f"Dose offset {a} (voxels)", -1.0, 1.0, 0.0, 0.25) for a in "xyz"]
    coarsen = st.selectbox(
        "Dose resolution",
        [1, 2, 3],
        format_func=lambda f: "x1 (native)" if f == 1 else f"coarsened x{f}",
    )
    interp = st.selectbox(
        "Interpolation", [1, 0, 3], format_func={0: "nearest", 1: "trilinear", 3: "cubic"}.get
    )
    sampling = st.selectbox("Structure sampling", ["native", "dose_center", "dose_fractional"])

steps = []
if sampling != "native":
    steps.append(StructureSampling(sampling))
if coarsen > 1:
    steps.append(DoseCoarsen(coarsen))
if any(off):
    steps.append(DoseGridOffset(tuple(off)))
if interp != 1:
    steps.append(InterpolationOrder(interp))
if any(t) or any(r):
    steps.append(RigidShift(tuple(t), tuple(r)))
if deform:
    steps.append(SmoothDeformation(deform, seed=0))
perturbed = Compose(*steps)(state)

planned_ev, new_ev = evaluate(state), evaluate(perturbed)
left, right = st.columns([3, 2])
with left:
    st.pyplot(
        viz.plot_dvh_overlay(
            {"planned": planned_ev, "perturbed": new_ev},
            rx,
            subtitle=perturbed.history and " + ".join(perturbed.history) or "no perturbation",
        )
    )
with right:
    a, b = planned_ev.flat(), new_ev.flat()
    df = pd.DataFrame({"planned": a, "perturbed": b})
    df["change"] = df["perturbed"] - df["planned"]
    st.dataframe(df.round(2), height=420, use_container_width=True)

st.subheader("What is this worth as a setup error?")
key = st.selectbox("Metric", ["PTV:D95%", "CTV:D98%", "Rectum:D2%", "Rectum:Dmean", "Bladder:D2%"])
delta = b[key] - a[key]
mm, direction = equivalent_shift_mm(sweeps, key, delta, min_delta=0.1)
if np.isinf(mm):
    st.write(f"{key} changed by {delta:+.2f}: more than any rigid shift up to 12 mm would cause.")
elif mm == 0:
    st.write(f"{key} changed by {delta:+.2f}: negligible.")
else:
    st.write(
        f"{key} changed by **{delta:+.2f}**, about what a **{mm:.1f} mm** rigid shift "
        f"({direction}) causes (magnitude-matched)."
    )
with st.expander("Plan geometry"):
    st.pyplot(viz.plot_slices(state, ct=(case.ct_grid, case.ct_hu)))
