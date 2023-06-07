import unittest
from optstoicpy.script.solver import (
    load_pulp_solver,
    ORDERED_SOLVERS,
    GLPK_CMD_OPTIONS,
    GUROBI_CMD_OPTIONS)


class TestSolver(unittest.TestCase):

    def test_at_least_one_pulp_solver_loading(self):
        pulp_solver = load_pulp_solver(solver_names=ORDERED_SOLVERS)
        self.assertNotEqual(pulp_solver, None)
        self.assertIn(pulp_solver.name, ORDERED_SOLVERS)

    def test_load_scip_cmd(self):
        solver = load_pulp_solver(solver_names=['SCIP_CMD'])

        if not solver:
            self.skipTest("SCIP_CMD is not available!")
        else:
            self.assertIn("-s", solver.options)

    def test_load_gurobi_cmd(self):
        solver = load_pulp_solver(solver_names=['GUROBI_CMD'])

        if not solver:
            self.skipTest("GUROBI_CMD is not available!")
        else:
            self.assertListEqual(GUROBI_CMD_OPTIONS, solver.options)

    def test_load_glpk_cmd(self):
        solver = load_pulp_solver(solver_names=['GLPK_CMD'])

        if not solver:
            self.skipTest("GLPK_CMD is not available!")
        else:
            self.assertListEqual(GLPK_CMD_OPTIONS, solver.options)
