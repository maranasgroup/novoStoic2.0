import sys
import unittest
from optstoicpy.core.database import (
    load_db_v3,
    Database)
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import (
    load_pulp_solver,
    ORDERED_SOLVERS)
import optstoicpy.script.optstoic_glycolysis as optsg
import optstoicpy.script.optstoic as opts


class TestOptStoic(unittest.TestCase):
    def setUp(self):
        self.logger = create_logger(name='Test generalized optstoic')
        self.pulp_solver = load_pulp_solver(
            solver_names=ORDERED_SOLVERS)
        if self.pulp_solver.name not in ['GUROBI', 'GUROBI_CMD', 'CPLEX_CMD']:
            self.skipTest("Skip because it will take a long time to solve.")

    def load_database(self):
        self.logger.info("Test loading Database")

        user_defined_export_rxns_Sji = {
            'EX_glc': {'C00031': -1.0},
            'EX_nad': {'C00003': -1.0},
            'EX_adp': {'C00008': -1.0},
            'EX_phosphate': {'C00009': -1.0},
            'EX_pyruvate': {'C00022': -1.0},
            'EX_nadh': {'C00004': -1.0},
            'EX_atp': {'C00002': -1.0},
            'EX_h2o': {'C00001': -1.0},
            'EX_hplus': {'C00080': -1.0},
            'EX_nadp': {'C00006': -1.0},
            'EX_nadph': {'C00005': -1.0}
        }

        DB = load_db_v3(
            reduce_model_size=True,
            user_defined_export_rxns_Sji=user_defined_export_rxns_Sji)
        DB.validate()
        # assert_equal(len(DB.metabolites), 5969)
        # assert_equal(len(DB.reactions), 7175)
        return DB

    def test_optstoic_glycolysis(self):
        """Test optstoic for glycolysis pathway."""
        model = optsg.OptStoicGlycolysis(
            objective='MinFlux',
            nATP=1,
            zlb=10,
            add_loopless_constraints=False, # Setting this to False, would speed up.
            max_iteration=1,
            pulp_solver=self.pulp_solver,
            result_filepath=None,
            M=1000)

        lp_prob, pathways = model.solve(
            outputfile='test_optstoic.txt')

        # if sys.platform == 'cygwin':
        #     lp_prob, pathways = model.solve_gurobi_cl(
        #         outputfile='test_optstoic_cyg.txt', cleanup=False)
        #     test.max_iteration = test.max_iteration + 2
        #     lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic_cyg.txt', exclude_existing_solution=True, cleanup=False)
        #else:
        lp_prob, pathways = model.solve(outputfile='test_optstoic.txt')
        #test.max_iteration = test.max_iteration + 1
        #lp_prob, pathways = test.solve(outputfile='test_optstoic.txt', exclude_existing_solution=True)

        self.assertEqual(pathways[1].note['modelstat'], 'Optimal')

    def test_general_optstoic(self):
        """Test optstoic analysis with standard setup

        How to add custom flux constraints:
        E.g.,
        v('EX_nadph') + v('EX_nadh') = 2;
        v('EX_nadp') + v('EX_nad') = -2;
        v('EX_nadh') + v('EX_nad') = 0;
        v('EX_nadph') + v('EX_nadp') = 0;
        became
        lp_prob += v['EX_nadph'] + v['EX_nadh'] == 2, 'nadphcons1'
        lp_prob += v['EX_nadp'] + v['EX_nad'] == -2, 'nadphcons2'
        lp_prob += v['EX_nadh'] + v['EX_nad'] == 0, 'nadphcons3'
        lp_prob += v['EX_nadph'] + v['EX_nadp'] == 0, 'nadphcons4'
        """
        self.DB = self.load_database()

        custom_flux_constraints = [
            {'constraint_name': 'nadphcons1',
             'reactions': ['EX_nadph', 'EX_nadh'],
             'UB': 2,
             'LB': 2},
            {'constraint_name': 'nadphcons2',
             'reactions': ['EX_nadp', 'EX_nad'],
             'UB': -2,
             'LB': -2},
            {'constraint_name': 'nadphcons3',
             'reactions': ['EX_nadh', 'EX_nad'],
             'UB': 0,
             'LB': 0},
            {'constraint_name': 'nadphcons4',
             'reactions': ['EX_nadph', 'EX_nadp'],
             'UB': 0,
             'LB': 0}]

        specific_bounds = {
            'EX_glc': {
                'LB': -1,
                'UB': -1},
            'EX_pyruvate': {
                'LB': 2,
                'UB': 2},
            'EX_nad': {
                'LB': -2,
                'UB': 0},
            'EX_nadh': {
                'LB': 0,
                'UB': 2},
            'EX_nadp': {
                'LB': -2,
                'UB': 0},
            'EX_nadph': {
                'LB': 0,
                'UB': 2},
            'EX_adp': {
                'LB': -1,
                'UB': -1},
            'EX_phosphate': {
                'LB': -1,
                'UB': -1},
            'EX_atp': {
                'LB': 1,
                'UB': 1},
            'EX_h2o': {
                'LB': 1,
                'UB': 1},
            'EX_hplus': {
                'LB': -10,
                'UB': 10}}  # pulp/gurobi has issue with "h+"

        model = opts.OptStoic(
            database=self.DB,
            objective='MinFlux',
            zlb=None,
            specific_bounds=specific_bounds,
            custom_flux_constraints=custom_flux_constraints,
            add_loopless_constraints=False,
            max_iteration=1,
            pulp_solver=self.pulp_solver,
            result_filepath='./result/',
            M=1000,
            logger=self.logger)

        # if sys.platform == 'cygwin':
        #     lp_prob, pathways = model.solve_gurobi_cl(
        #         outputfile='test_optstoic_general_cyg.txt', cleanup=False)
        #     test.max_iteration = test.max_iteration + 2
        #     lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic_general_cyg.txt', exclude_existing_solution=True, cleanup=False)
        # else:
        lp_prob, pathways = model.solve(
            outputfile='test_optstoic_general.txt')

        self.assertEqual(pathways[1].note['modelstat'], 'Optimal')
