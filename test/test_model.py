from modelling.model import Model
import pytest

@pytest.fixture
def model():
    return Model("test_model", [], [], [])


def test_create_variables(model):
    model.create_variables()
    # Add assertions to check if the variables are created correctly

def test_create_constraints(model):
    model.create_constraints()
    # Add assertions to check if the constraints are created correctly

def test_create_objective(model):
    model.create_objective()
    # Add assertions to check if the objective is created correctly

def test_solve(model):
    model.solve()
    # Add assertions to check if the model is solved correctly

def test_print_solution(model):
    model.print_solution()
    # Add assertions to check if the solution is printed correctly

def test_run(model):
    model.run()
    # Add assertions to check if the entire model execution is correct
