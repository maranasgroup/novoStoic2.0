# /usr/bin/python
"""
Loopless OptStoic program to identify glycolytic pathway
(glucose to pyruvate) for n ATP production.
It read input files that are used for GAMS.
Currently, it has been tested with SCIP, GLPK, Gurobi and CPLEX solvers.

Tip:
Change the number of Thread for gurobi if needed.

"""
from __future__ import absolute_import
import os
import time
import sys
import copy
import random
import string  # to generate random hex code
import pulp
# import cPickle as pickle
from . import gams_parser
import json
#import pdb
from optstoicpy.core import database
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import load_pulp_solver
from optstoicpy.script.optstoic import OptStoic
from .gurobi_command_line_solver import *


class OptStoicGlycolysis(OptStoic):

    EXPORT_RXNS_SJI = {
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

    CUSTOM_REDOX_CONSTRAINTS = [
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

    GLYCOLYSIS_BOUNDS = {'EX_glc': {'LB': -1, 'UB': -1},
                         'EX_pyruvate': {'LB': 2, 'UB': 2},
                         'EX_nad': {'LB': -2, 'UB': 0},
                         'EX_nadh': {'LB': 0, 'UB': 2},
                         'EX_nadp': {'LB': -2, 'UB': 0},
                         'EX_nadph': {'LB': 0, 'UB': 2},
                         # 'EX_adp': {'LB': -1, 'UB': -1},
                         # 'EX_phosphate': {'LB': -1, 'UB': -1},
                         # 'EX_atp': {'LB': 1, 'UB': 1},
                         # 'EX_h2o': {'LB': 1, 'UB': 1},
                         'EX_hplus': {'LB': -10, 'UB': 10}}

    def __init__(self,
                 objective='MinFlux',
                 zlb=10,
                 nATP=1,
                 add_loopless_constraints=True,
                 max_iteration=1,
                 pulp_solver=None,
                 result_filepath=None,
                 M=1000,
                 logger=None):
        """An example of the optStoic model for identifying glycolytic pathways
            generating n ATP.

        Args:
            objective (str, optional): Description
            zlb (int, optional): Description
            nATP (int, optional): nATP (int, optional): The number of ATP
            add_loopless_constraints (bool, optional): Description
            max_iteration (int, optional): Description
            pulp_solver (None, optional): Description
            result_filepath (None, optional): Description
            M (int, optional): Description
            logger (None, optional): Description
        """

        self.DBV3 = database.load_db_v3(
            reduce_model_size=True,
            user_defined_export_rxns_Sji=self.EXPORT_RXNS_SJI
        )

        super(OptStoicGlycolysis, self).__init__(
            database=self.DBV3,
            objective='MinFlux',
            zlb=zlb,
            specific_bounds=self.GLYCOLYSIS_BOUNDS,
            custom_flux_constraints=self.CUSTOM_REDOX_CONSTRAINTS,
            add_loopless_constraints=add_loopless_constraints,
            max_iteration=max_iteration,
            pulp_solver=pulp_solver,
            result_filepath=result_filepath,
            M=1000,
            logger=logger)

        self.nATP = nATP

    @property
    def nATP(self):
        return self._nATP

    @nATP.setter
    def nATP(self, value):
        self._nATP = value

        self.specific_bounds['EX_atp'] = {'LB': value, 'UB': value}
        self.specific_bounds['EX_h2o'] = {'LB': value, 'UB': value}
        self.specific_bounds['EX_adp'] = {'LB': -value, 'UB': -value}
        self.specific_bounds['EX_phosphate'] = {'LB': -value, 'UB': -value}

        # When nATP is not integer, change variables v, vf and vb to continuous
        # variables
        if float(value).is_integer():
            self._varCat = 'Integer'
        else:
            self._varCat = 'Continuous'

    def __repr__(self):
        return "<OptStoicGlycolysis(nATP='%s', objective='%s')>" % (
            self.nATP, self.objective)


def test_optstoic_glycolysis():
    logger = create_logger(name='optstoicpy.script.optstoic_glycolysis.main')
    #logger.debug('Testing optstoic output filepath: %s', res_dir)
    logger.info("Test optstoic_glycolysis")

    pulp_solver = load_pulp_solver(
        solver_names=[
            'SCIP_CMD',
            'GUROBI',
            'GUROBI_CMD',
            'CPLEX_CMD',
            'GLPK_CMD'],
        logger=logger)

    test = OptStoicGlycolysis(
        objective='MinFlux',
        nATP=1,
        zlb=10,  # setting this may slow down the optimization, but integer cut constraints will work
        max_iteration=1,
        pulp_solver=pulp_solver,
        result_filepath='./result/',
        M=1000,
        logger=logger)

    if sys.platform == 'cygwin':
        lp_prob, pathways = test.solve_gurobi_cl(
            outputfile='test_optstoic_cyg.txt', cleanup=False)
        #test.max_iteration = test.max_iteration + 2
        #lp_prob, pathways = test.solve_gurobi_cl(outputfile='test_optstoic_cyg.txt', exclude_existing_solution=True, cleanup=False)
    else:
        lp_prob, pathways = test.solve(outputfile='test_optstoic.txt')
        #test.max_iteration = test.max_iteration + 1
        #lp_prob, pathways = test.solve(outputfile='test_optstoic.txt', exclude_existing_solution=True)

    return lp_prob, pathways


if __name__ == '__main__':
    lp_prob, pathways = test_optstoic_glycolysis()
