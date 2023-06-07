import unittest
from optstoicpy.core.database import (
    load_custom_reactions_to_be_excluded,
    load_base_reaction_db
)
from optstoicpy.script.database_preprocessing import blocked_reactions_analysis
from optstoicpy.script.solver import (
    load_pulp_solver,
    ORDERED_SOLVERS)


class TestDatabasePreprocessing(unittest.TestCase):
    def test_blocked_reactions_analysis(self):
        """Test blocked reactions analysis.
        A Pulp solver must be loaded for this test to complete.
        Otherwise this test will be skipped.
        """
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

        db = load_base_reaction_db(
            user_defined_export_rxns_Sji=user_defined_export_rxns_Sji
        )

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

        specific_bounds = {'EX_glc': {'LB': -1, 'UB': -1},
                           'EX_pyruvate': {'LB': 2, 'UB': 2},
                           'EX_nad': {'LB': -2, 'UB': 0},
                           'EX_nadh': {'LB': 0, 'UB': 2},
                           'EX_nadp': {'LB': -2, 'UB': 0},
                           'EX_nadph': {'LB': 0, 'UB': 2},
                           'EX_adp': {'LB': -5, 'UB': -1},  # range 1-5 ATP
                           'EX_phosphate': {'LB': -5, 'UB': -1},
                           'EX_atp': {'LB': 1, 'UB': 5},  # range 1-5 ATP
                           'EX_h2o': {'LB': 1, 'UB': 5},
                           'EX_hplus': {'LB': -10, 'UB': 10}}  # pulp/gurobi has issue with "h+"

        pulp_solver = load_pulp_solver(
            solver_names=ORDERED_SOLVERS,
            logger=None)

        if not pulp_solver:
            self.skipTest(
                "Solver is not found. Please install at least one solver.")

        exclude_reactions = load_custom_reactions_to_be_excluded()

        blocked_reactions_list, FVA_res = blocked_reactions_analysis(
            database=db,
            pulp_solver=pulp_solver,
            specific_bounds=specific_bounds,
            custom_flux_constraints=custom_flux_constraints,
            excluded_reactions=exclude_reactions,
            target_reactions_list=['R01266', 'R07882', 'R00658', 'R01059'])

        self.assertEqual(set(blocked_reactions_list),
                         set(['R01266', 'R07882']))
