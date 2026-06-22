
import math
from requirement_scratch import calc_mass
from partition import trans_parti, parti, rate
from rate_calc import top_type

class DummyReader:
    def __init__(self, lines):
        self.lines = [str(x).split() for x in lines]
        self.idx = 0

    def next_data_line(self):
        value = self.lines[self.idx]
        self.idx += 1
        return value


def test_calc_mass_diatomic():
    reader = DummyReader([2, 1.0, 35.0])
    num_atom, masses, mass_tot, mass_tot_kg = calc_mass(1, reader)

    assert num_atom == 2
    assert masses == [1.0, 35.0]
    assert mass_tot == 36.0
    assert mass_tot_kg > 0


def test_trans_partition_positive():
    q = trans_parti(1.0e-26, 1.0e-3, 300.0)
    assert q > 0


def test_total_partition():
    q = parti(10.0, 20.0, 30.0)
    assert q == 6000.0


def test_rate_positive():
    k = rate(1.0, 10.0, 300.0)
    assert k > 0


def test_top_type_single_atom():
    assert top_type(0, 0.0, 0.0, 0.0) == "single atom"
