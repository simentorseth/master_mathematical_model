from modelling.model import Model
from modelling.utils import get_locations, get_demands, get_parameters


def main():
    instance = [
        # ["id", "kind", "x", "y", "demand", "earliest", "latest"],
        ["1", "station", "0.0", "0.0", "0.0", "0.0", "10.0"],
        ["2", "station", "0.5", "0.5", "0.0", "0.0", "10.0"],
        ["3", "customer", "1.0", "1.0", "5.0", "3.0", "8.0"],
        ["4", "customer", "2.0", "2.0", "2.0", "2.0", "4.0"],
        ["5", "customer", "3.0", "3.0", "3.0", "6.0", "9.0"],
        ["charge_capacity", "/100.0/"],
        ["average_velocity", "/50.0/"],
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

    model.run()


if __name__ == "__main__":
    main()
