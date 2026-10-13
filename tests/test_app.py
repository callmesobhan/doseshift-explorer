from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = Path(__file__).parent.parent / "app" / "streamlit_app.py"


def test_app_runs_and_reacts_to_a_shift():
    at = AppTest.from_file(str(APP), default_timeout=120).run()
    assert not at.exception
    sliders = {s.label: s for s in at.sidebar.slider}
    sliders["Shift y (mm)"].set_value(3.0).run()
    assert not at.exception
    assert any("mm" in m.value and "rigid shift" in m.value for m in at.markdown)
