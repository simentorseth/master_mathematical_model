import itertools
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from modelling.utils import compute_distance
from termcolor import colored, cprint


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
        # self.time_window_parameters = {i: [0.0, 10.0] for i in self.time_windows}
        self.time_window_parameters = self.create_time_windows()
        self.max_cumulative_travel_time = parameters["max_cumulative_travel_time"]
        self.time_window_max = max(
            [value[1] for value in self.time_window_parameters.values()]
        )
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
                gp.quicksum(
                    self.x[i, j, k]
                    for i in self.turbines.keys()
                    for k in self.vessels
                    if i != j
                )
                <= self.y[j] * len(self.turbines)
                for j in self.destinations.keys()
            ),
            name="station_assignment",
        )
        # self.model.addConstrs(
        #     (
        #         self.e_d[j, k] <= self.charge_capacity * self.y[j]
        #         for j in self.destinations.keys()
        #         for k in self.vessels
        #     )
        # )
        # vessel assignment
        self.model.addConstrs(
            (
                gp.quicksum(
                    self.x[i, j, k]
                    for i in self.vertices_w_origin
                    for j in self.vertices_w_destination
                    if i != j
                )
                <= self.z[k] * (len(self.turbines) + 1)
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
                    self.b[k, p] * self.time_window_parameters[p][1]  # upper bound
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
                - self.time_window_max * (1 - self.x[i, j, k])
                for i, j in self.arcs
                for k in self.vessels
                if i != j
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
                - (self.x[i, j, k] * self.distances[i, j] * self.consumption_rate)
                + self.charge_capacity * (1 - self.x[i, j, k])
                for i in self.vertices_w_origin
                for j in self.vertices_w_destination
                for k in self.vessels
                if i != j
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

        self.model.addConstrs(
            (
                self.e_d[j, k] >= 0
                for j in self.vertices_w_origin_and_destination
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

    def create_time_windows(self, window_length=8):
        result = []
        for i in range(self.num_time_windows):
            if i == 0:
                start = 0
            else:
                start = result[-1][1] + 1
            end = start + window_length
            result.append([start, end])
        return {i: v for i, v in enumerate(result)}

    def solve(self):
        # self.model.Params.LazyConstraints = 1
        self.model.optimize()

    def print_solution(self):
        for k in self.vessels:
            for i, j in self.arcs:
                if self.x[i, j, k].x > 0.5:
                    print(
                        colored(
                            f"x[{i}, {j}, {k}] = {self.x[i, j, k].x}",
                            "light_green",
                            attrs=["bold"],
                        )
                    )
                # else:
                #     cprint(
                #         f"x[{i}, {j}, {k}] = {self.x[i, j, k].x}", "red", attrs=["bold"]
                #     )
        # for i, j, k in self.x.keys():
        #     if self.x[i, j, k].x > 0.5:
        #         print(
        #             colored(
        #                 f"x[{i}, {j}, {k}] = {self.x[i, j, k].x}",
        #                 "light_green",
        #                 attrs=["bold"],
        #             )
        #         )
        #     else:
        #         cprint(f"x[{i}, {j}, {k}] = {self.x[i, j, k].x}", "red", attrs=["bold"])
        # print("\n")
        # y[j]
        for j in self.y.keys():
            if self.y[j].x > 0.5:
                print(colored(f"y[{j}] = {self.y[j].x}", "light_green", attrs=["bold"]))
            else:
                cprint(f"y[{j}] = {self.y[j].x}", "red", attrs=["bold"])
        print("\n")
        # z[k]
        for k in self.z.keys():
            if self.z[k].x > 0.5:
                print(colored(f"z[{k}] = {self.z[k].x}", "light_green", attrs=["bold"]))
            else:
                cprint(f"z[{k}] = {self.z[k].x}", "red", attrs=["bold"])
        print("\n")
        # b[k, p]
        for k, p in self.b.keys():
            if self.b[k, p].x > 0.5:
                print(
                    colored(
                        f"b[{k}, {p}] = {self.b[k, p].x}", "light_green", attrs=["bold"]
                    )
                )
            else:
                cprint(f"b[{k}, {p}] = {self.b[k, p].x}", "red", attrs=["bold"])
        print("\n")
        # t[i, k]
        for k in self.vessels:
            if self.z[k].x > 0.5:
                for i in self.vertices_w_origin_and_destination:
                    print(
                        colored(
                            f"t[{i}, {k}] = {self.t[i, k].x}",
                            "light_green",
                            attrs=["bold"],
                        )
                    )
            else:
                print("\n")
        print("\n")
        # for i, k in self.t.keys():
        #     if self.t[i, k].x > 0.0:
        #         print(
        #             colored(
        #                 f"t[{i}, {k}] = {self.t[i, k].x}",
        #                 "light_green",
        #                 attrs=["bold"],
        #             )
        #         )
        #     else:
        #         cprint(f"t[{i}, {k}] = {self.t[i, k].x}", "red", attrs=["bold"])
        print("\n")
        # e_a[i, k]
        for i, k in self.e_a.keys():
            print(f"e_a[{i}, {k}] = {self.e_a[i, k].x}")
        print("\n")
        # e_d[i, k]
        for i, k in self.e_d.keys():
            print(f"e_d[{i}, {k}] = {self.e_d[i, k].x}")

        for i, window in self.time_window_parameters.items():
            print(f"T_{i}: {window}")

    def visualise_solution(self):
        G = nx.DiGraph()
        G.add_edges_from((i, j) for i, j, k in self.x.keys() if self.x[i, j, k].x > 0.5)
        # nodes
        options = {"edgecolors": "black", "node_size": 300}
        nx.draw_networkx_nodes(
            G,
            self.turbines,
            nodelist=self.turbines.keys(),
            node_color="tab:blue",
            label="turbines",
            **options,
        )
        nx.draw_networkx_nodes(
            G,
            self.destinations,
            nodelist=[
                j
                for j in self.y.keys()
                if self.y[j].x > 0.5 and j in self.destinations.keys()
            ],
            node_color="tab:green",
            label="built substation",
            **options,
        )
        nx.draw_networkx_nodes(
            G,
            self.destinations,
            nodelist=[j for j in self.y.keys() if self.y[j].x < 0.5],
            node_color="tab:gray",
            label="unbuilt substation",
            **options,
        )
        # edges
        nx.draw_networkx_edges(G, self.nodes, edgelist=G.edges(), edge_color="tab:gray")
        # labels
        nx.draw_networkx_labels(G, self.nodes, font_size=12, font_color="white")
        # plt.legend()
        plt.show()

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
        if self.model.Status != 3:  # infeasible
            self.print_solution()
            self.visualise_solution()
        else:
            print("Model is proven infeasible by Gurobi")
        # Free resources safely
        self.model.dispose()
        gp.disposeDefaultEnv()
