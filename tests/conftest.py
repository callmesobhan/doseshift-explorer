import pytest

from doseshift.synthetic import make_case


@pytest.fixture(scope="session")
def case():
    return make_case(seed=0)


@pytest.fixture(scope="session")
def state(case):
    return case.to_state()
