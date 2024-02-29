import itertools
import gurobipy as gp
from gurobipy import GRB
from modelling.utils import compute_distance


class Model:
    def __init__(
        self,
        *,
        name: str,
        turbines: list[tuple[float, float]],
        stations: list[tuple[float, float]],
        vessels: list[int],
        turbine_demands: list[float],
        station_demands: list[float],
        parameters: dict[str, float],
        num_time_windows: int = 7,
    ):
        self.model = gp.Model(name)
        # Sets
        self.turbines = {k: v for k, v in enumerate(turbines)}
        # self.stations = {len(self.turbines) + k: v for k, v in enumerate(stations)}  # Set of station locations, (x, y)
        self.origins = {len(self.turbines) + k: v for k, v in enumerate(stations)}
        self.destinations = {
            len(self.turbines) + k + len(stations): v for k, v in enumerate(stations)
        }
        self.nodes = self.turbines | self.origins | self.destinations
        # self.vertices = set(self.turbines.keys()).union(set(self.stations.keys()))
        self.vertices_w_origin = set(self.turbines.keys()).union(
            set(self.origins.keys())
        )
        self.vertices_w_destination = set(self.turbines.keys()).union(
            set(self.destinations.keys())
        )
        self.vertices_w_origin_and_destination = (
            set(self.turbines.keys())
            .union(set(self.origins.keys()))
            .union(set(self.destinations.keys()))
        )
        self.vessels = vessels  # Set of vessels
        self.arcs = [
            (i, j)
            for i, j in itertools.product(
                self.vertices_w_origin_and_destination, repeat=2
            )
            if i != j
        ]
        self.num_time_windows = num_time_windows  # Number of time windows
        self.time_windows = [num for num in range(self.num_time_windows)]
        # Parameters
        self.fixed_costs_station = {
            j: parameters["fixed_costs_station"] for j in self.destinations.keys()
        }
        self.fixed_costs_vessel = {
            k: parameters["fixed_costs_vessel"] for k in self.vessels
        }
        self.turbine_demands = {
            key: turbine_demands[i] for i, key in enumerate(self.turbines.keys())
        }
        # self.station_demands = {
        #     key: station_demands[i] for i, key in enumerate(self.origins.keys())
        # }
        self.origin_demands = {
            key: station_demands[i] for i, key in enumerate(self.origins.keys())
        }
        self.destination_demands = {
            key: station_demands[i] for i, key in enumerate(self.destinations.keys())
        }
        self.demands = (
            self.turbine_demands | self.origin_demands | self.destination_demands
        )
        self.distances = {
            (i, j): compute_distance(self.nodes[i], self.nodes[j]) for i, j in self.arcs
        }
        self.travel_times = {
            (i, j): self.distances[i, j] / parameters["average_velocity"]
            for i, j in self.arcs
        }
        self.consumption_rate = parameters["consumption_rate"]
        self.charge_capacity = parameters["charge_capacity"]
        self.time_window_parameters = {i: [0.0, 10.0] for i in self.time_windows}
        self.max_cumulative_travel_time = parameters["max_cumulative_travel_time"]
        self.big_M = 1000  # Big M
        # Variables
        self.x = None
        self.y = None
        self.z = None
        # Auxiliary variables
        self.b = None
        self.t = None
        self.e_a = None
        self.e_d = None

    def create_variables(self):
        self.x = self.model.addVars(
            [(i, j, k) for (i, j), k in itertools.product(self.arcs, self.vessels)],
            vtype=GRB.BINARY,
            name="x",
        )
        self.y = self.model.addVars(
            self.destinations.keys(),
            vtype=GRB.BINARY,
            name="y",
        )
        self.z = self.model.addVars(self.vessels, vtype=GRB.BINARY, name="z")
        self.b = self.model.addVars(
            [(k, p) for k in self.vessels for p in self.time_windows],
            vtype=GRB.BINARY,
            name="b",
        )
        self.t = self.model.addVars(
            [
                (i, k)
                for i, k in itertools.product(
                    self.vertices_w_origin_and_destination, self.vessels
                )
            ],
            vtype=GRB.CONTINUOUS,
            name="t",
        )
        self.e_a = self.model.addVars(
            [
                (i, k)
                for i, k in itertools.product(
                    self.vertices_w_origin_and_destination, self.vessels
                )
            ],
            vtype=GRB.CONTINUOUS,
            name="e_a",
        )
        self.e_d = self.model.addVars(
            [
                (i, k)
                for i, k in itertools.product(
                    self.vertices_w_origin_and_destination, self.vessels
                )
            ],
            vtype=GRB.CONTINUOUS,
            name="e_d",
        )

    def create_constraints(self):
        # single assignment
        self.model.addConstrs(
            (
                gp.quicksum(
                    self.x[i, j, k]
                    for k in self.vessels
                    for j in self.vertices_w_destination
                    if i != j
                )
                == 1
                for i in self.turbines.keys()
            ),
            name="turbine_assignment",
        )
        # flow conservation 1
        self.model.addConstrs(
            (
                gp.quicksum(self.x[i, j, k] for i in self.turbines.keys() if i != j)
                <= 1
                for j in self.destinations
                for k in self.vessels
            ),
            name="start_and_end_at_station",
        )
        # flow conservation 2
        self.model.addConstrs(
            (
                gp.quicksum(
                    self.x[i, j, k] for j in self.vertices_w_destination if i != j
                )
                - gp.quicksum(self.x[j, i, k] for j in self.vertices_w_origin if i != j)
                == 0
                for k in self.vessels
                for i in self.turbines.keys()
            ),
            name="flow_conservation",
        )
        # station assignment
        self.model.addConstrs(
            (
                gp.quicksum(self.x[i, j, k] for i in self.turbines.keys() if i != j)
                <= self.y[j] * len(self.arcs)
                for j in self.destinations.keys()
                for k in self.vessels
            ),
            name="station_assignment",
        )
        # vessel assignment
        self.model.addConstrs(
            (
                gp.quicksum(self.x[i, j, k] for i, j in self.arcs if i != j)
                <= self.z[k] * len(self.arcs)
                for k in self.vessels
            ),
            name="vessel assignment",
        )
        # time window assignment
        self.model.addConstrs(
            (
                gp.quicksum(self.b[k, p] for p in self.time_windows) == 1
                for k in self.vessels
            ),
            name="time_window_assignment",
        )
        # time window constraints
        self.model.addConstrs(
            (
                gp.quicksum(
                    self.b[k, p] * self.time_window_parameters[p][0]  # lower bound
                    for p in self.time_windows
                )
                <= self.t[i, k]
                for i in self.vertices_w_origin_and_destination
                for k in self.vessels
            ),
            name="time_window_earliest",
        )
        self.model.addConstrs(
            (
                gp.quicksum(
                    self.b[k, p] * self.time_window_parameters[p][-1]  # upper bound
                    for p in self.time_windows
                )
                >= self.t[i, k]
                for i in self.vertices_w_origin_and_destination
                for k in self.vessels
            ),
            name="time_window_latest",
        )
        # maximum cumulative travel time
        # self.model.addConstrs((gp.quicksum((self.demands[i] + self.travel_times[i, j]) * self.x[i, j, k] for i, j in self.arcs) <= gp.quicksum(self.b[k, p] * self.max_cumulative_travel_time for p in self.time_windows[i]) for i in self.stations.keys() for k in self.vessels), name="cumulative_travel_time")
        # arrival time
        self.model.addConstrs(
            (
                self.t[j, k]
                >= self.t[i, k]
                + (self.demands[i] + self.travel_times[i, j]) * self.x[i, j, k]
                - self.big_M * (1 - self.x[j, i, k])
                for i, j in self.arcs
                for k in self.vessels
            ),
            name="arrival_time",
        )
        # arrival time
        self.model.addConstrs(
            (
                self.t[i, k] >= 0
                for i in self.vertices_w_origin_and_destination
                for k in self.vessels
            ),
            name="non_negative_arrival_time",
        )

        # remaining charge at arrival
        self.model.addConstrs(
            (
                self.e_a[j, k]
                <= self.e_d[i, k]
                - self.x[i, j, k] * self.distances[i, j] * self.consumption_rate
                + self.big_M * (1 - self.x[i, j, k])
                for i, j in self.arcs
                for k in self.vessels
            ),
            name="remaining_charge_arrival",
        )
        # no energy loss at turbines
        self.model.addConstrs(
            (
                self.e_d[i, k] == self.e_a[i, k]
                for i in self.turbines.keys()
                for k in self.vessels
            ),
            name="no_loss_at_node",
        )
        # charge within capacity
        self.model.addConstrs(
            (
                self.e_a[i, k] >= 0
                for i in self.vertices_w_origin_and_destination
                for k in self.vessels
            ),
            name="charge_within_capacity",
        )
        # full charge at station
        self.model.addConstrs(
            (
                self.e_d[j, k] == self.charge_capacity
                for j in self.origins.keys()
                for k in self.vessels
            ),
            name="full_charge_at_station",
        )

    def create_objective(self):
        self.model.setObjective(
            gp.quicksum(
                self.fixed_costs_station[j] * self.y[j]
                for j in self.destinations.keys()
            )
            + gp.quicksum(self.fixed_costs_vessel[k] * self.z[k] for k in self.vessels),
            GRB.MINIMIZE,
        )

    # def subtour_elimination(self, model, where):
    #     if where == GRB.Callback.MIPSOL:
    #         # make a list of edges selected in the solution
    #         xs = model.cbGetSolution(self.model._x)
    #         selected = gp.tuplelist((i, j, k) for i, j, k in model._x.keys() if xs[i, j, k] > 0.5)
    #         # find the shortest cycle in the selected edge list
    #         tour = self.subtour(selected)
    #         if len(tour) < len(self.turbines.keys()) - 2:
    #             # add subtour elimination constr. for every pair of cities in subtour
    #             [(i, j, k) for (i, j), k in itertools.product(self.arcs, self.vessels)]
    #             model.cbLazy(gp.quicksum(model._x[i, j, k] for (i, j), k in itertools.product(tour, self.vessels)) <= len(tour)-1)

    # def subtour(self, edges):
    #     unvisited = self.turbines.keys()[:]
    #     cycle = self.turbines.keys()[:] # Dummy - guaranteed to be replaced
    #     while unvisited:  # true if list is non-empty
    #         thiscycle = []
    #         neighbors = unvisited
    #         while neighbors:
    #             current = neighbors[0]
    #             thiscycle.append(current)
    #             unvisited.remove(current)
    #             neighbors = [j for i, j in edges.select(current, '*')
    #                         if j in unvisited]
    #         if len(thiscycle) <= len(cycle):
    #             cycle = thiscycle # New shortest subtour
    #     return cycle

    def solve(self):
        # self.model.Params.LazyConstraints = 1
        self.model.optimize()

    def print_solution(self):
        for vessel, time_window in self.b.keys():
            if self.b[vessel, time_window].x > 0.5:
                print(f"Vessel {vessel} used time window {time_window}")
                for origin, destination in self.arcs:
                    if self.x[origin, destination, vessel].x > 0.5:
                        arrival_origin = self.t[origin, vessel].x
                        arrival_destination = self.t[destination, vessel].x
                        print(
                            f"Vessel {vessel} traveled from {origin} at time {arrival_origin} and arrived at {destination} at time {arrival_destination}"
                        )

    def run(self):
        self.create_variables()
        self.create_constraints()
        self.create_objective()

        assert self.x is not None
        assert self.y is not None
        assert self.z is not None
        assert self.t is not None
        assert self.b is not None
        assert self.e_a is not None
        assert self.e_d is not None

        self.solve()
        self.print_solution()
        # Free resources safely
        if self.model.Status == 3:
            print("Model is proven infeasible by Gurobi")
        self.model.dispose()
        gp.disposeDefaultEnv()
