from modelling.model import Model
from modelling.utils import get_locations, get_demands, get_parameters
from pprint import pprint


def time_window_generator(num_windows, window_length=8):
    result = []
    for i in range(num_windows):
        if i == 0:
            start = 0
        else:
            start = result[-1][1] + 1
        end = start + window_length
        result.append([start, end])
    return result


def main():
    instance = [
        # ["id", "kind", "x", "y", "demand", "earliest", "latest"],
        ["1", "station", "40.0", "50.0", "0.0", "0.0", "10.0"],
        ["2", "station", "70.0", "20.0", "0.0", "0.0", "10.0"],
        ["3", "customer", "25.0", "41.0", "1.0", "3.0", "8.0"],
        ["4", "customer", "62.0", "43.0", "1.5", "2.0", "4.0"],
        ["5", "customer", "34.0", "33.0", "0.5", "6.0", "9.0"],
        ["6", "customer", "10.0", "35.0", "1.0", "6.0", "9.0"],
        ["charge_capacity", "/50.0/"],
        ["average_velocity", "/10.0/"],
        ["consumption_rate", "/0.1/"],
        ["fixed_costs_station", "/1000.0/"],
        ["fixed_costs_vessel", "/1000.0/"],
        ["max_cumulative_travel_time", "/8.0/"],
    ]
    customers = get_locations(instance, "customer")
    stations = get_locations(instance, "station")
    vessels = [1, 2, 3]
    turbine_demands = get_demands(instance, "customer")
    station_demands = get_demands(instance, "station")
    parameters = get_parameters(instance)

    model = Model(
        name="test_model",
        turbines=customers,
        stations=stations,
        vessels=vessels,
        turbine_demands=turbine_demands,
        station_demands=station_demands,
        parameters=parameters,
    )
    print("DISTANCES")
    pprint(model.distances)
    print("\n")
    print("TRAVEL TIMES")
    pprint(model.travel_times)
    print("\n")
    print("DEMANDS")
    pprint(model.demands)

    model.run()


if __name__ == "__main__":
    main()
