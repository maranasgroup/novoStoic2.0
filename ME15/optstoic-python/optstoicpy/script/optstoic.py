# /usr/bin/python
"""
Loopless OptStoic program to identify any pathway.
It read input files that are used for GAMS.
Currently, it has been tested with SCIP, GLPK (without loopless constraints), Gurobi and CPLEX solvers.

To-do list:
2. generalize the input for designing other pathways
3. Fix MinRxn
4. Add feature to run MinRxn for a range of zlb

Known issues:
1. MinRxn takes a very long time to solve (gap = 60% after 20 mins)
2. Pulp is having issue with gurobi700

Tip:
Change the number of Thread for gurobi if needed.

"""
from __future__ import absolute_import
from builtins import range
from builtins import object
import os
import time
import sys
import copy
import random
import string  # to generate random hex code
import json
import pulp
# import cPickle as pickle
#import pdb
from optstoicpy.core.database import load_db_v3
from optstoicpy.core.pathway import Pathway
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import load_pulp_solver
from .gurobi_command_line_solver import *

# Global variables/solver options
EPS = 1e-5
GUROBI_OPTIONS = 'Threads=2 TimeLimit=1800 MIPGapAbs=1e-6 MIPGap=1e-6 CliqueCuts=2'


class OptStoic(object):
    """
    An OptStoic problem Class to identify pathways
    using either minFlux or minRxn algorithm.
    """

    def __init__(self,
                 database,
                 objective='MinFlux',
                 zlb=None,
                 specific_bounds=None,
                 custom_flux_constraints=None,
                 add_loopless_constraints=True,
                 max_iteration=2,
                 pulp_solver=None,
                 result_filepath=None,
                 M=1000,
                 logger=None):
        """
        Args:
            database (TYPE): An optStoic Database object (equivalent to GSM model)
            objective (str, optional): The mode for optStoic pathway prospecting.
                Options available are: ['MinFlux', 'MinRxn']
            zlb (int, optional): The lowerbound on objective value z. If not provided,
                the integer cut constraints may not work properly when the objective value
                increases.
            specific_bounds (dict, optional): LB and UB for exchange reactions which defined the
                overall pathway equations. E.g. {'Ex_glc': {'LB': -1, 'UB':-1}}
            custom_flux_constraints (dict, optional): The custom constraints that need to be
                added to the model formulation.
            add_loopless_constraints (bool, optional): If True, use loopless constraints.
                If False, run optStoic without loopless constraint (faster, but the pathway
                may contain loops).
            max_iteration (int, optional): The default maximum number of iteration
            pulp_solver (None, optional): A pulp.solvers object (load any of the user-defined
                solver)
            result_filepath (str, optional): Filepath for result
            M (int, optional): The maximum flux bound (default 1000)
            logger (:obj:`logging.Logger`, optional): A logging.Logger object

        Raises:
            Exception: Description
        """
        if logger is None:
            self.logger = create_logger(name='optstoic.OptStoic')
        else:
            self.logger = logger

        self.objective = objective
        self.zlb = zlb

        if specific_bounds is None:
            raise Exception("specific_bounds must be specified!")
        self.specific_bounds = specific_bounds
        self.add_loopless_constraints = add_loopless_constraints
        self.custom_flux_constraints = custom_flux_constraints
        self.M = M

        self._varCat = 'Integer'
        # self._varCat = 'Continuous'

        self.max_iteration = max_iteration

        if result_filepath is None:
            result_filepath = './result'
        self.result_filepath = result_filepath

        if not os.path.exists(self.result_filepath):
            self.logger.warning(
                "A folder %s is created!" %
                self.result_filepath)
            os.makedirs(self.result_filepath)

        self.database = database
        self.pathways = {}
        self.iteration = 1
        self.lp_prob = None
        self.pulp_solver = pulp_solver
        self.lp_prob_fname = "OptStoic_{0}".format(
            self.generate_random_string(6))

    @staticmethod
    def generate_random_string(N):
        """LP file is appended with random string when using command line mode
            to prevent overwrite/read issues.
        """
        return ''.join(
            random.choice(
                string.ascii_uppercase +
                string.digits) for _ in range(N))

    def change_objective(self, new_objective):
        if new_objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is "
                             "not correctly defined. "
                             "Please use either 'MinFlux' or 'MinRxn'.")
        self.objective = new_objective

    def change_zlb(self, zlb):
        self.zlb = zlb

    def create_minflux_problem(self):
        """
        Create minflux/minRxn LP problem (a pulp LpProblem object)
        minimize z = sum(j\j_exchange, vf(j) + vb(j))
        subject to:
            z = zlb
            v(j) = vf(j) - vb(j)
            sum(j, S(i,j) * v(j)) = 0, for all i
            vf(j) >= yf(j) * eps, for all j
            vf(j) <= yf(j) * M, for all j
            vb(j) >= yb(j) * eps, for all j
            vb(j) <= yb(j) * M, for all j
            yf(j) + yb(j) <= 1, for all j

            for reaction in jloop (not in jblock)
            sum(j, Nint(l,j) * G(j)) = 0, for all l
            G(j) >= -M * a(j) + (1 - a(j))
            G(j) <= -a(j) + M * (1 - a(j))
            v(j) >= -M * (1 - a(j))
            v(j) <= M * a(j)

        """
        self.logger.info("Formulating problem...")
        # Scalar
        M = self.M

        # Initialize variables
        v = pulp.LpVariable.dicts("v", self.database.reactions,
                                  lowBound=-M, upBound=M, cat=self._varCat)
        vf = pulp.LpVariable.dicts("vf", self.database.reactions,
                                   lowBound=0, upBound=M, cat=self._varCat)
        vb = pulp.LpVariable.dicts("vb", self.database.reactions,
                                   lowBound=0, upBound=M, cat=self._varCat)
        yf = pulp.LpVariable.dicts("yf", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        yb = pulp.LpVariable.dicts("yb", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')

        if self.add_loopless_constraints:
            a = pulp.LpVariable.dicts("a", self.database.reactions,
                                      lowBound=0, upBound=1, cat='Binary')
            G = pulp.LpVariable.dicts("G", self.database.reactions,
                                      lowBound=-M, upBound=M, cat='Continuous')
        else:
            a = None
            G = None

        # Update lower and upper bound based on reaction directionality

        for j in self.database.reactions:

            if j in self.database.all_excluded_reactions:
                v[j].lowBound = 0
                v[j].upBound = 0
                vf[j].lowBound = 0
                vf[j].upBound = 0
                yf[j].lowBound = 0
                yf[j].upBound = 0
                yb[j].lowBound = 0
                yb[j].upBound = 0

                continue

            if self.database.rxntype[j] == 0:
                # Forward irreversible
                v[j].lowBound = 0
                v[j].upBound = M
                yb[j].upBound = 0
                vb[j].upBound = 0

            elif self.database.rxntype[j] == 1:
                # Reversible
                v[j].lowBound = -M
                v[j].upBound = M

            elif self.database.rxntype[j] == 2:
                # Reverse irreversible
                v[j].lowBound = -M
                v[j].upBound = 0
                vf[j].upBound = 0
                yf[j].upBound = 0

            elif self.database.rxntype[j] == 4:
                v[j].lowBound = 0
                v[j].upBound = 0

        # Fix stoichiometry of source/sink metabolites
        for rxn, bounds in self.specific_bounds.items():
            v[rxn].lowBound = bounds['LB']
            v[rxn].upBound = bounds['UB']

        LB = {}
        UB = {}

        for j in self.database.reactions:
            LB[j] = v[j].lowBound
            UB[j] = v[j].upBound

        lp_prob = pulp.LpProblem("OptStoic", pulp.LpMinimize)

        # Min-Rxn objective
        if self.objective == 'MinRxn':
            condition = pulp.lpSum([yf[j] + yb[j]
                                    for j in self.database.reactions
                                    if self.database.rxntype[j] != 4])
            lp_prob += condition, "MinRxn"

        # Min-Flux objective
        elif self.objective == 'MinFlux':
            condition = pulp.lpSum([vf[j] + vb[j]
                                    for j in self.database.reactions
                                    if self.database.rxntype[j] != 4])
            lp_prob += condition, "MinFlux"

            if self.zlb is not None:
                # fix lower bound
                lp_prob += condition == self.zlb, 'zLowerBound'

        # Constraints
        # Mass_balance
        for i in self.database.metabolites:
            # If metabolites not involve in any reactions
            if i not in self.database.S:
                continue
            label = "mass_balance_%s" % i
            dot_S_v = pulp.lpSum([self.database.S[i][j] * v[j]
                                  for j in list(self.database.S[i].keys())])
            condition = dot_S_v == 0
            lp_prob += condition, label

        # if self.objective == 'MinRxn':
        # for j in self.database.reactions:
        #     lp_prob += v[j] >= y[j]*LB[j], "cons1_%s"%j
        #     lp_prob += v[j] >= y[j]*LB[j], "cons1_%s"%j
        #     lp_prob += v[j] <= y[j]*UB[j], "cons2_%s"%j

        if self.objective == 'MinFlux':
            for j in self.database.reactions:
                lp_prob += (v[j] == vf[j] - vb[j]), "flux_%s" % j

                # These constraints ensure that when yf=0 and yb=0 ,
                # no flux goes through the reaction
                lp_prob += vf[j] >= yf[j] * 0.5, "cons1_%s" % j
                lp_prob += vf[j] <= yf[j] * M, "cons2_%s" % j
                lp_prob += vb[j] >= yb[j] * 0.5, "cons3_%s" % j
                lp_prob += vb[j] <= yb[j] * M, "cons4_%s" % j
                # Ensure that either yf or yb can be 1, not both
                lp_prob += yf[j] + yb[j] <= 1, 'cons5_%s' % j

        if self.add_loopless_constraints:
            self.logger.info("Loopless constraints are turned on.")

            loop_rxn = list(set(self.database.internal_rxns) -
                            set(self.database.blocked_rxns))

            # Loopless contraints
            for l in self.database.loops:
                label = "loopless_cons_%s" % l
                dot_N_G = pulp.lpSum([self.database.Ninternal[l][j] * G[j]
                                      for j in list(self.database.Ninternal[l].keys())])
                condition = dot_N_G == 0
                lp_prob += condition, label

            for j in loop_rxn:
                lp_prob += G[j] >= -M * a[j] + (1 - a[j]), "llcons1_%s" % j
                lp_prob += G[j] <= -a[j] + M * (1 - a[j]), "llcons2_%s" % j
                lp_prob += v[j] >= -M * (1 - a[j]), "llcons3_%s" % j
                lp_prob += v[j] <= M * a[j], "llcons4_%s" % j

        # Fix nad(p)h production and consumption
        if self.custom_flux_constraints is not None:
            self.logger.info("Adding custom constraints...")

            for group in self.custom_flux_constraints:
                lp_prob += pulp.lpSum(v[rxn] for rxn in group['reactions']
                                      ) <= group['UB'], "%s_UB" % group['constraint_name']
                lp_prob += pulp.lpSum(v[rxn] for rxn in group['reactions']
                                      ) >= group['LB'], "%s_LB" % group['constraint_name']

        return lp_prob, v, vf, vb, yf, yb, a, G

    def solve(
            self,
            exclude_existing_solution=False,
            outputfile="OptStoic_pulp_result.txt",
            max_iteration=None):
        """
        Solve OptStoic problem using pulp.solvers interface

        Args:
            exclude_existing_solution (bool, optional): If True, create and add integer cut
                constraints for pathways that are found using the same OptStoic instance,
                but solved in previous function call.
            outputfile (str, optional): name of outpufile
            max_iteration (None, optional): Externally specified maximum number of pathway to be
                found using OpStoic. If not specified, it will set to the internal max iterations.

        Returns:
            TYPE: Description

        Raises:
            ValueError: Description
        """
        if self.objective not in ['MinFlux', 'MinRxn']:
            raise ValueError(
                "The objective for OptStoic is not correctly defined. Please use either 'MinFlux' or 'MinRxn'.")

        if max_iteration is None:
            max_iteration = self.max_iteration

        self.logger.info(
            "Finding multiple pathways using Optstoic %s...",
            self.objective)
        lp_prob, v, vf, vb, yf, yb, a, G = self.create_minflux_problem()

        # Create integer cut for existing pathways
        if exclude_existing_solution and bool(self.pathways):
            self.iteration = max(self.pathways.keys()) + 1
            if self.iteration > max_iteration:
                raise ValueError('Max iteration is less than current '
                                 'iteration. Increase max_iteration '
                                 'before solving!')

            for ind, pathway in self.pathways.items():
                rxnlist = list(set(pathway.reaction_ids_no_exchange))
                condition = pulp.lpSum(
                    [(1 - yf[j] - yb[j]) for j in rxnlist]) >= 1
                lp_prob += condition, "IntegerCut_%d" % ind

        self.logger.info("Solving problem...")
        # if self.iteration == 1:
        #     result_output = open(os.path.join(self.result_filepath, outputfile), "w+")
        # else:
        #     result_output = open(os.path.join(self.result_filepath, outputfile), "a+")

        while True and self.iteration <= max_iteration:
            self.logger.info("Iteration %s", self.iteration)
            # lp_prob.writeLP("OptStoic.lp", mip=1)  # optional
            e1 = time.time()
            lp_prob.solve(solver=self.pulp_solver)
            e2 = time.time()
            self.logger.info(
                "This iteration solved in %.3f seconds.",
                (e2 - e1))

            # The solution is printed if it was deemed "optimal
            if pulp.LpStatus[lp_prob.status] == "Optimal":
                self.logger.info("Writing result to output file...")
                # result_output.write("\nIteration no.: %d\n" %self.iteration)
                # result_output.write("\nModelstat: %s\n" %pulp.LpStatus[lp_prob.status])
                res = {}
                res['reaction_id'] = []
                res['flux'] = []
                res['iteration'] = self.iteration
                res['time'] = (e2 - e1)
                res['modelstat'] = "Optimal"

                for j in self.database.reactions:
                    if v[j].varValue is not None:
                        if v[j].varValue > EPS or v[j].varValue < -EPS:
                            res['reaction_id'].append(j)
                            res['flux'].append(v[j].varValue)
                #             result_output.write("%s %.8f\n" %(v[j].name, v[j].varValue))

                # result_output.write("%s = %.8f\n" % (self.objective, pulp.value(lp_prob.objective)))
                # result_output.write("----------------------------------\n\n")

                integer_cut_reactions = list(
                    set(res['reaction_id']) - set(self.database.user_defined_export_rxns))

                self.pathways[self.iteration] = Pathway(
                    id=self.iteration,
                    name='Pathway_{:03d}'.format(self.iteration),
                    reaction_ids=res['reaction_id'],
                    fluxes=res['flux'],
                    sourceSubstrateID='C00031',
                    endSubstrateID='C00022',
                    note=res
                )

                self.write_pathways_to_json(json_filename="temp_pathways.json")

                # Integer cut constraint is added so that
                # the same solution cannot be returned again
                condition = pulp.lpSum([(1 - yf[j] - yb[j])
                                        for j in integer_cut_reactions]) >= 1
                lp_prob += condition, "IntegerCut_%d" % self.iteration
                self.iteration += 1

            # If a new optimal solution cannot be found, end the program
            else:
                break

        # result_output.close()

        self.lp_prob = lp_prob

        return self.lp_prob, self.pathways

    def write_pathways_to_json(self, json_filename="temp_pathways.json"):

        temp = {}
        for k in list(self.pathways.keys()):
            temp[k] = self.pathways[k].to_dict()

        json.dump(
            temp,
            open(
                os.path.join(
                    self.result_filepath,
                    json_filename),
                'w+'),
            sort_keys=True,
            indent=4)

    def add_existing_pathways(self, user_defined_pathways):
        """
        Add list of existing solutions (Pathways) to be
        excluded from being identified.

        Args:
            user_defined_pathways (TYPE): pathways output from solve_gurobi_cl()
                                     or solve()

        Raises:
            ValueError: Description
        """
        if (isinstance(user_defined_pathways, dict) and (
                isinstance(list(user_defined_pathways.values())[0], Pathway))):
            self.pathways = copy.deepcopy(user_defined_pathways)
        else:
            raise ValueError(
                "user_defined_pathways must be a dictionary of Pathway instances")

    def reset_pathways(self):
        """
        Reset self.pathways to empty dictionary
        """
        self.pathways = {}

    def solve_gurobi_cl(self,
                        exclude_existing_solution=False,
                        outputfile="OptStoic_pulp_result_gcl.txt",
                        max_iteration=None,
                        cleanup=True,
                        gurobi_options=GUROBI_OPTIONS):
        """
        Solve OptStoic problem using Gurobi command line (gurobi_cl)
        when pulp.solvers.GUROBI_CMD failed.
        Require the module "gurobi_command_line_solver.py".

        Args:
            exclude_existing_solution (bool, optional): If true and if self.pathway is not None,
                exclude the pathways from being identified.
            outputfile (str, optional): name of outpufile
            max_iteration (None, optional): Externally specified maximum number of pathway
                to be found using OpStoic. If not specified,
                it will set to the internal max iterations.
            cleanup (bool, optional): If True, delete the temporary .lp and .sol file. Set as
                False for debugging.
            gurobi_options (TYPE, optional): Description

        Returns:
            TYPE: Description

        Raises:
            ValueError: Description
        """
        if self.objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is not correctly "
                             "defined. Please use either 'MinFlux' or "
                             "'MinRxn'.")

        if max_iteration is None:
            max_iteration = self.max_iteration

        t1 = time.time()

        self.logger.info("Finding multiple pathways using"
                         " Optstoic %s and Gurobi CL...", self.objective)
        lp_prob, v, vf, vb, yf, yb, a, G = self.create_minflux_problem()

        # Create integer cut for existing pathways
        if exclude_existing_solution and bool(self.pathways):
            self.iteration = max(self.pathways.keys()) + 1
            if self.iteration > max_iteration:
                raise ValueError('Max iteration is less than current '
                                 'iteration. Increase max_iteration '
                                 'before solving!')

            for ind, pathway in self.pathways.items():
                rxnlist = list(set(pathway.reaction_ids_no_exchange))
                condition = pulp.lpSum(
                    [(1 - yf[j] - yb[j]) for j in rxnlist]) >= 1
                lp_prob += condition, "IntegerCut_%d" % ind

        # Solve problem
        self.logger.info("Solving problem...")

        # if self.iteration == 1:
        #     result_output = open(os.path.join(
        #         self.result_filepath, outputfile), "w+")
        # else:
        #     result_output = open(os.path.join(
        #         self.result_filepath, outputfile), "a+")

        while True and self.iteration <= max_iteration:
            self.logger.info("Iteration %s", self.iteration)
            lp_prob.writeLP(self.lp_prob_fname + ".lp", mip=1)
            e1 = time.time()
            lp_status, solver_message = solve_with_gurobi_cl_debug(
                self.lp_prob_fname, options=gurobi_options)
            e2 = time.time()
            self.logger.info(
                "This iteration solved in %.3f seconds.",
                (e2 - e1))

            # The solution is printed if it was deemed "optimal
            if lp_status in ["Optimal", "Time_limit"]:
                objective_function, varValue = parse_gurobi_sol(
                    self.lp_prob_fname)

                res = {}
                res['reaction_id'] = []
                res['flux'] = []
                res['iteration'] = self.iteration
                res['time'] = (e2 - e1)
                res['modelstat'] = lp_status
                res['solvestat'] = solver_message

                # result_output.write("\nIteration no.: %d\n" %self.iteration)
                # result_output.write("\nModelstat: %s\n" %lp_status)

                for j in self.database.reactions:
                    if 'v_' + j in varValue:
                        v = varValue['v_' + j]
                        if v > EPS or v < -EPS:
                            res['reaction_id'].append(j)
                            res['flux'].append(v)
                            #result_output.write("%s %.8f\n" %(j, v))

                # result_output.write("%s = %.8f\n" %(self.objective, objective_function))
                # result_output.write("----------------------------------\n\n")

                integer_cut_reactions = list(
                    set(res['reaction_id']) - set(self.database.user_defined_export_rxns))

                self.pathways[self.iteration] = Pathway(
                    id=self.iteration,
                    name='Pathway_{:03d}'.format(self.iteration),
                    reaction_ids=res['reaction_id'],
                    fluxes=res['flux'],
                    sourceSubstrateID='C00031',
                    endSubstrateID='C00022',
                    note=res
                )
                # Keep a copy of pathways in case program terminate midway
                self.write_pathways_to_json(json_filename="temp_pathways.json")

                # Integer cut constraint is added so that
                # the same solution cannot be returned again
                condition = pulp.lpSum([(1 - yf[j] - yb[j])
                                        for j in integer_cut_reactions]) >= 1
                lp_prob += condition, "IntegerCut_%d" % self.iteration
                self.iteration += 1

            # If a new optimal solution cannot be found, end the program
            else:
                break

        # result_output.close()
        # Clean up directory
        if cleanup:
            self.logger.debug("Cleaning up directory...")
            os.remove("./" + self.lp_prob_fname + ".lp")
            os.remove("./" + self.lp_prob_fname + ".sol")
            os.remove("./gurobi.log")

        self.lp_prob = lp_prob

        return self.lp_prob, self.pathways

    def __repr__(self):
        return "<OptStoic(objective='%s')>" % (self.objective)


def test_optstoic():
    """An alternative to the nosetest due to issue with PULP/SCIP_CMD
    """
    logger = create_logger(name='optstoicpy.script.optstoic.main')

    logger.info("Test generalized optstoic")

    # Set the following reactions as allowable export reactions
    db3 = load_db_v3(
        reduce_model_size=True,
        user_defined_export_rxns_Sji={
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
    )

    #logger.debug('Testing optstoic output filepath: %s', res_dir)

    pulp_solver = load_pulp_solver(
        solver_names=[
            'SCIP_CMD',
            'GUROBI',
            'GUROBI_CMD',
            'CPLEX_CMD',
            'GLPK_CMD'],
        logger=logger)

    # How to add custom flux constraints:
    # E.g.,
    # v('EX_nadph') + v('EX_nadh') = 2;
    # v('EX_nadp') + v('EX_nad') = -2;
    # v('EX_nadh') + v('EX_nad') = 0;
    # v('EX_nadph') + v('EX_nadp') = 0;
    # became
    # lp_prob += v['EX_nadph'] + v['EX_nadh'] == 2, 'nadphcons1'
    # lp_prob += v['EX_nadp'] + v['EX_nad'] == -2, 'nadphcons2'
    # lp_prob += v['EX_nadh'] + v['EX_nad'] == 0, 'nadphcons3'
    # lp_prob += v['EX_nadph'] + v['EX_nadp'] == 0, 'nadphcons4'

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
                       'EX_adp': {'LB': -1, 'UB': -1},
                       'EX_phosphate': {'LB': -1, 'UB': -1},
                       'EX_atp': {'LB': 1, 'UB': 1},
                       'EX_h2o': {'LB': 1, 'UB': 1},
                       'EX_hplus': {'LB': -10, 'UB': 10}}  # pulp/gurobi has issue with "h+"

    test = OptStoic(database=db3,
                    objective='MinFlux',
                    zlb=None,
                    specific_bounds=specific_bounds,
                    custom_flux_constraints=custom_flux_constraints,
                    add_loopless_constraints=True,
                    max_iteration=2,
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
    lp_prob, pathways = test()
