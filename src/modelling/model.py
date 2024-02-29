import itertools
import gurobipy as gp
from gurobipy import GRB
from modelling.utils import compute_distance


class Model:
    def __init__(self, *, 
                 name: str, 
                 turbines: list[tuple[float, float]],
                 stations: list[tuple[float, float]],
                 vessels: list[int],
                 turbine_demands: list[float],
                 station_demands: list[float],
                 parameters: dict[str, float],
                 planning_horizon: int = 7,
                 ):
        self.model = gp.Model(name)
        self.planning_horizon = planning_horizon
        # Sets
        self.turbines = {k: v for k, v in enumerate(turbines)}  # Set of turbine locations, (x, y)
        self.stations = {len(self.turbines) + k: v for k, v in enumerate(stations)}  # Set of station locations, (x, y)
        self.num_time_windows = {i: self.planning_horizon for i in self.turbines.keys()}
        self.nodes = self.turbines | self.stations  # Set of nodes (turbines and stations), (x, y)
        self.vessels = vessels  # Set of vessels
        self.vertices = set(self.turbines.keys()).union(set(self.stations.keys()))  # Set of vertices (turbines and stations)
        self.arcs = [(i, j) for i, j in itertools.product(self.vertices, repeat=2) if i != j]
        self.time_windows = {i: [num for num in range(self.num_time_windows[i])] for i in self.turbines.keys()}
        # Parameters
        self.fixed_costs_station = {j: parameters["fixed_costs_station"] for j in self.stations.keys()}  # Fixed costs of stations
        self.fixed_costs_vessel = {k: parameters["fixed_costs_vessel"] for k in self.vessels}  # Fixed costs of vessels
        self.turbine_demands = {key: turbine_demands[i] for i, key in enumerate(self.turbines.keys())} # Demand of turbines
        self.station_demands = {key: station_demands[i] for i, key in enumerate(self.stations.keys())} # Demand of stations
        self.demands = self.turbine_demands | self.station_demands  # Demand of all nodes
        self.distances = {(i, j): compute_distance(self.nodes[i], self.nodes[j]) for i, j in self.arcs}  # Distance between nodes
        self.travel_times = {(i, j): self.distances[i, j] / parameters["average_velocity"] for i, j in self.arcs}  # Travel time between nodes
        self.consumption_rate = parameters["consumption_rate"]  # Consumption rate of vessels
        self.charge_capacity = parameters["charge_capacity"]  # Charge capacity of vessels
        self.time_window_parameters = {i: [[0.0, 10.0] for _ in range(self.num_time_windows[i])] for i in self.turbines.keys()}  # Time window parameters
        self.max_cumulative_travel_time = parameters["max_cumulative_travel_time"]  # Maximum cumulative travel time
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
        self.x = self.model.addVars([(i, j, k) for (i, j), k in itertools.product(self.arcs, self.vessels)], vtype=GRB.BINARY, name="x")
        self.y = self.model.addVars(self.stations.keys(), vtype=GRB.BINARY, name="y")
        self.z = self.model.addVars(self.vessels, vtype=GRB.BINARY, name="z")
        self.b = self.model.addVars([(k, p) for k in self.vessels for p in range(self.planning_horizon)], vtype=GRB.BINARY, name="b")
        self.t = self.model.addVars([(i, k) for i, k in itertools.product(self.vertices, self.vessels)], vtype=GRB.CONTINUOUS, name="t")
        self.e_a = self.model.addVars([(i, k) for i, k in itertools.product(self.nodes, self.vessels)], vtype=GRB.CONTINUOUS, name="e_a")
        self.e_d = self.model.addVars([(i, k) for i, k in itertools.product(self.nodes, self.vessels)], vtype=GRB.CONTINUOUS, name="e_d")
        
    def create_constraints(self):
        self.model.addConstrs((gp.quicksum(self.x[i, j, k] for j in self.vertices if i != j) == 1 for i in self.turbines.keys() for k in self.vessels), name="turbine_assignment")
        self.model.addConstrs((gp.quicksum(self.x[i ,j ,k] for i in self.vertices if i != j) - gp.quicksum(self.x[j, i, k] for i in self.vertices if i != j) == 0 for j in self.vertices for k in self.vessels), name="flow_conservation")
        self.model.addConstrs((gp.quicksum(self.x[i, j ,k] for i in self.vertices if i != j) <= self.y[j] * self.big_M for j in self.stations.keys() for k in self.vessels), name="station_assignment")
        self.model.addConstrs((gp.quicksum(self.x[i, j, k] for i in self.vertices for j in self.vertices if i != j) <= self.z[k] * self.big_M for k in self.vessels), name="vessel assignment")
        
        self.model.addConstrs((gp.quicksum(self.x[i, j, k] for i in self.vertices if i != j) <= 1 for j in self.turbines.keys() for k in self.vessels),name="start_and_end_at_station")

        self.model.addConstrs((gp.quicksum(self.b[k, p] for p in range(self.planning_horizon)) == 1 for k in self.vessels), name="time_window_assignment")
        self.model.addConstrs((gp.quicksum(self.b[k, p] * self.time_window_parameters[i][p][0] for p in self.time_windows[i]) <= self.t[i, k] for i in self.turbines.keys() for k in self.vessels), name="time_window_earliest")
        self.model.addConstrs((gp.quicksum(self.b[k, p] * self.time_window_parameters[i][p][-1] for p in self.time_windows[i]) >=self.t[i, k] for i in self.turbines.keys() for k in self.vessels), name="time_window_latest")
        
        #self.model.addConstrs((gp.quicksum((self.demands[i] + self.travel_times[i, j]) * self.x[i, j, k] for i, j in self.arcs) <= gp.quicksum(self.b[k, p] * self.max_cumulative_travel_time for p in self.time_windows[i]) for i in self.stations.keys() for k in self.vessels), name="cumulative_travel_time")

        self.model.addConstrs((self.t[j, k] >= self.t[i, k] + self.demands[i] * self.travel_times[i, j] * self.x[i, j, k] - self.big_M * (1 - self.x[j, i, k]) for i, j in self.arcs for k in self.vessels), name="arrival_time")

        self.model.addConstrs((self.e_a[j, k] <= self.e_d[i, k] - self.x[i, j, k] * self.distances[i, j] * self.consumption_rate + self.big_M * (1 - self.x[i, j, k]) for i, j in self.arcs for k in self.vessels), name="remaining_charge_arrival")
        self.model.addConstrs((self.e_d[i, k] == self.e_a[i, k] for i in self.turbines.keys() for k in self.vessels), name="no_loss_at_node")
        self.model.addConstrs((self.e_a[i, k] >= 0 for i in self.vertices for k in self.vessels), name="charge_within_capacity")
        self.model.addConstrs((self.e_d[j, k] == self.charge_capacity for j in self.stations.keys() for k in self.vessels), name="full_charge_at_station")

    def create_objective(self):
        self.model.setObjective(gp.quicksum(self.fixed_costs_station[j] * self.y[j] for j in self.stations.keys()) + gp.quicksum(self.fixed_costs_vessel[k] * self.z[k] for k in self.vessels), GRB.MINIMIZE)

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
        #self.model.Params.LazyConstraints = 1
        self.model.optimize()


    def print_solution(self):
        routes = []
        for origin, destination, vessel in self.x.keys():
            if self.x[origin, destination, vessel].x > 0.5:
                routes.append((origin, destination, vessel))

        for origin, destination, vessel in sorted(routes, key=lambda x: x[2]):
                print(f"Vessel {vessel} travels from {origin} to {destination}")
        
        for station in self.y.keys():
            if self.y[station].x > 0.5:
                print(f"Station {station} is open")

        for vessel in self.z.keys():
            if self.z[vessel].x > 0.5:
                print(f"Vessel {vessel} is used")

        for vessel, time in self.b.keys():
            if self.b[vessel, time].x > 0.5:
                for i in self.turbines.keys():
                    print(f"Vessel {vessel} is active at time {self.time_window_parameters[i][time]} and visited turbine {i}")

        for node, vessel in self.t.keys():
            print(f"Vessel {vessel} arrives at node {node} at time {self.t[node, vessel].x}")

        for node, vessel in self.e_a.keys():
            print(f"Vessel {vessel} has {self.e_a[node, vessel].x} charge at node {node}")

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

        self.model.Params.Presolve = 0
        self.solve()
        self.print_solution()
        # Free resources safely
        if self.model.Status == 3:
            print("Model is proven infeasible by Gurobi")
        self.model.dispose()
        gp.disposeDefaultEnv()