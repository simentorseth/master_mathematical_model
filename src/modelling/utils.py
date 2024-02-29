from math import sqrt
import re


def get_instance(filename: str) -> list[list[str]]:
    with open(filename, "r") as file:
        lines = file.readlines()
        lines = [line.strip().split() for line in lines]
    instance = lines[1:] # skip header
    return instance
       


def compute_distance(loc1: tuple[float, float], loc2: tuple[float, float]) -> float:
    dx = loc1[0] - loc2[0]
    dy = loc1[1] - loc2[1]
    return sqrt(dx * dx + dy * dy)


def get_locations(instance: list[list[str]], kind: str) -> list[tuple[float, float]]:
    locations = []
    for row in instance:
        if row:
            if row[1] == kind:
                x, y = float(row[2]), float(row[3])
                locations.append((x, y))
    return locations


def get_time_windows(instance: list[list[str]], kind: str) -> list[tuple[float, float]]:
    time_windows = []
    for row in instance:
        if row:
            if row[1] == kind:
                earliest, latest = float(row[5]), float(row[6])
                time_windows.append((earliest, latest))
    return time_windows


def get_service_times(instance: list[list[str]], kind: str) -> list[float]:
    service_times = []
    for row in instance[1:]:  # skip header
        if row:
            if row[1] == kind:
                service_time = float(row[-1])
                service_times.append(service_time)
    return service_times


def get_demands(instance: list[list[str]], kind: str) -> list[float]:
    demands = []
    for row in instance:
        if row:
            if row[1] == kind:
                demand = float(row[4])
                demands.append(demand)
    return demands


def get_parameters(instance: list[list[str]]) -> dict[str, float]:
    parameters = {}
    for row in instance[-6:]:
        key = row[0]
        if (match := re.search(r'/(?P<number>\d+\.\d+|\d+)/', row[-1])):
            value = match.group("number")
            parameters[key] = float(value)
    return parameters
