from modelling.utils import get_locations, get_time_windows, get_service_times, get_demands, get_parameters, compute_distance
import pytest


@pytest.fixture
def instance():
    return [
        ["id", "kind", "x", "y", "demand", "earliest", "latest", "service_time"],
        ["1", "depot", "0.0", "0.0", "0.0", "0.0", "10.0", "0.0"],
        ["2", "station", "0.0", "0.0", "0.0", "0.0", "10.0", "0.0"],
        ["3", "customer", "1.0", "1.0", "1.0", "3.0", "8.0", "1.0"],
        ["4", "customer", "2.0", "2.0", "2.0", "2.0", "4.0", "2.0"],
        ["5", "customer", "3.0", "3.0", "3.0", "6.0", "9.0", "3.0"],
        ["fuel_capacity", "/100.0/"],
        ["average_velocity", "/50.0/"],
    ]


def test_compute_distance():
    loc1 = (0.0, 0.0)
    loc2 = (0.0, 0.0)
    assert compute_distance(loc1, loc2) == 0.0


def test_get_locations_depot(instance):
    kind = "depot"
    expected = [(0.0, 0.0)]
    assert get_locations(instance, kind) == expected


def test_get_locations_customer(instance):
    kind = "customer"
    expected = [(1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    assert get_locations(instance, kind) == expected


def test_get_locations_station(instance):
    kind = "station"
    expected = [(0.0, 0.0)]
    assert get_locations(instance, kind) == expected


def test_get_time_windows(instance):
    kind = "customer"
    expected = [(3.0, 8.0), (2.0, 4.0), (6.0, 9.0)]
    assert get_time_windows(instance, kind) == expected


def test_get_service_times(instance):
    kind = "customer"
    expected = [1.0, 2.0, 3.0]
    assert get_service_times(instance, kind) == expected


def test_get_demands(instance):
    kind = "customer"
    expected = [1.0, 2.0, 3.0]
    assert get_demands(instance, kind) == expected


def test_get_parameters(instance):
    expected = {
        "fuel_capacity": 100.0,
        "average_velocity": 50.0,
    }
    assert get_parameters(instance) == expected