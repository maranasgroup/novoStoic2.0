import os
import copy
import json
import pulp
from nose.tools import (
    assert_equal)
import numpy as np
import scipy.io
import pandas as pd
from sympy import (
    Matrix,
    S,
    nsimplify
)
from sympy.matrices import SparseMatrix
from optstoicpy.core.database import (
    load_custom_reactions_to_be_excluded,
    load_base_reaction_db
)
from optstoicpy.script.utils import create_logger
from optstoicpy.script.solver import load_pulp_solver


def blocked_reactions_analysis(
        database,
        pulp_solver,
        specific_bounds,
        custom_flux_constraints,
        excluded_reactions=None,
        target_reactions_list=None,
        logger=None):
    """
    Perform flux variability analysis on the database,
    based on the overall reaction equation of optstoic.
    If a reaction cannot carry flux (i.e., -eps <= v(j) <= eps, where eps = 1e-8),
    then the reaction is considered as a blocked reaction.
    The blocked reactions are then eliminated from the database S matrix.
    Next, the internal loops (excluding cofactors) are identified.
    Then, optStoic analysis can be performed for pathway prospecting.

    max/min v(j)

    subject to:
            sum(j, S(i,j) * v(j)) = 0, for all i
            custom_flux_constraints

        Note: The glycolysis study was done using the GAMS version of this code.
        This is written in attempt to port find_blocked_reactions.gms from GAMS to Python,
        as a part of effort to generalize optstoic analysis.

    Args:
        database (:obj:`BaseReactionDatabase`): The default reaction database
            without blocked reactions/loops.
        pulp_solver (TYPE): The solver for PuLP.
        specific_bounds (dict): LB and UB for exchange reactions which defined the
            overall design equations. E.g. {'Ex_glc': {'LB': -1, 'UB':-1}}
        custom_flux_constraints (TYPE): The custom constraints that need to be
            added to the model formulation.
        excluded_reactions (None, optional): The list of reactions that are manually
            selected to be excluded from optstoic solution.
        target_reactions_list (None, optional): If provided, the blocked reaction analysis is performed
            only on a subset of the reaction provided. If None, the blocked reaction analysis
            will be performed on all reactions in the database. The excluded_reactions set
            can be subtracted(e.g., set(database.reactions) - excluded_reactions), since
            they are blocked reactions.
        logger (:obj:`logging.logger`, optional): The logging instance

    Returns:
        TYPE: Description

    Raises:
        ValueError: Description

    Deleted Parameters:
        user_defined_export_rxns_Sji (dict): The list of export reactions that
            need to be added to the model for metabolite exchange (i.e., any metabolite
            that participate in the design equation)
    """
    if logger is None:
        logger = create_logger(
            name="optstoicpy.script.database_preprocessing.blocked_reactions_analysis")

    logger.warning(
        "This process may take a long time to run. It is recommended to be run in a batch script.")

    M = 1000
    EPS = 1e-8

    # Initialize variables
    v = pulp.LpVariable.dicts("v", database.reactions,
                              lowBound=-M, upBound=M, cat='Continuous')

    for j in database.reactions:
        if database.rxntype[j] == 0:
            # Forward irreversible
            v[j].lowBound = 0
            v[j].upBound = M

        elif database.rxntype[j] == 1:
            # Reversible
            v[j].lowBound = -M
            v[j].upBound = M

        elif database.rxntype[j] == 2:
            # Reverse irreversible
            v[j].lowBound = -M
            v[j].upBound = 0

        elif database.rxntype[j] == 4:
            v[j].lowBound = 0
            v[j].upBound = 0

        else:
            raise ValueError("Reaction type for reaction %s is unknown." % j)

    if excluded_reactions is not None:
        for j in excluded_reactions:
            v[j].lowBound = 0
            v[j].upBound = 0

    # Fix stoichiometry of source/sink metabolites
    for j, bounds in specific_bounds.items():
        v[j].lowBound = bounds['LB']
        v[j].upBound = bounds['UB']

    FVA_res = {}
    blocked_reactions = []
    lp_prob = None

    if target_reactions_list is None:
        target_reactions_list = database.reactions
    num_rxn = len(target_reactions_list)

    for ind, j1 in enumerate(target_reactions_list):
        logger.debug("%s/%s" % (ind, num_rxn))
        FVA_res[j1] = {}

        for obj in ['min', 'max']:

            # Variables (make a copy)
            vt = copy.deepcopy(v)
            del lp_prob

            # Objective function
            if obj == 'min':
                lp_prob = pulp.LpProblem("FVA%s" % obj, pulp.LpMinimize)
                lp_prob += vt[j1], "FVA_min"
            elif obj == 'max':
                lp_prob = pulp.LpProblem("FVA%s" % obj, pulp.LpMaximize)
                lp_prob += vt[j1], "FVA_max"

            # Constraints
            # Mass_balance
            for i in database.metabolites:
                # If metabolites not involve in any reactions
                if i not in database.S:
                    continue
                label = "mass_balance_%s" % i
                dot_S_v = pulp.lpSum([database.S[i][j] * vt[j]
                                      for j in list(database.S[i].keys())])
                condition = dot_S_v == 0
                lp_prob += condition, label

            if custom_flux_constraints is not None:
                logger.info("Adding custom constraints...")

                for group in custom_flux_constraints:
                    lp_prob += pulp.lpSum(vt[rxn] for rxn in group['reactions']
                                          ) <= group['UB'], "%s_UB" % group['constraint_name']
                    lp_prob += pulp.lpSum(vt[rxn] for rxn in group['reactions']
                                          ) >= group['LB'], "%s_LB" % group['constraint_name']

            lp_prob.solve(solver=pulp_solver)

            FVA_res[j1][obj] = pulp.value(lp_prob.objective)

        if (FVA_res[j1]['max'] < EPS) and (FVA_res[j1]['min'] > -EPS):
            blocked_reactions.append(j1)

        json.dump(FVA_res,
                  open("temp_FVA_result.json", 'w+'),
                  sort_keys=True,
                  indent=4)

    return blocked_reactions, FVA_res


def remove_cofactors_from_Sij(Sij_df, cofactors):
    """
    Remove row of cofactors i from Sij matrix.
    Remove reaction j that involved only cofactors from Sij matrix.

    Args:
        Sij_df (TYPE): Description
        cofactors (TYPE): Description

    Returns:
        TYPE: Description
    """
    if len(cofactors) == 0:
        return Sij_df

    # Get a list of cofactors in the model
    cofactors = list(set(cofactors) & set(Sij_df.index.tolist()))

    # Remove row of cofactors
    nSij_df = Sij_df.drop(cofactors)

    allRxns = nSij_df.columns.tolist()

    # Get all columns (j) with all zero entries
    rxns_involving_cofactors_only = nSij_df.columns[(
        nSij_df == 0).all()].tolist()

    remainRxns = list(set(allRxns) - set(rxns_involving_cofactors_only))

    # Drop all columns with zero entries
    nSij_df2 = nSij_df[sorted(remainRxns)]

    return nSij_df2


def internal_loop_analysis(S_df, logger=None):
    """
    Identifies the "rational" basis for the null space of S_df matrix,
    and convert them to internal loops.
    This is an alternative to the Matlab function
    null(S_df, 'r') (which is significantly faster).
    Warning: This is extremely time consuming. Please use the MATLAB
        version.

    M = Matrix([[16, 2, 3,13],
    [5,11,10, 8],
    [9, 7, 6,12],
    [4,14,15, 1]])

    print(nsimplify(M, rational=True).nullspace())
    """
    if logger is None:
        logger = create_logger(
            name="optstoicpy.script.database_preprocessing.internal_loop_analysis")

    raise NotImplementedError("Use the Matlab version as this is too slow.")

    reactions = S_df_no_cofactor.columns.tolist()
    metabolites = S_df_no_cofactor.index.tolist()
    Sint = S_df_no_cofactor.as_matrix()

    # Sint_mat = Matrix(Sint) # too slow
    Smat = SparseMatrix(Sint.astype(int))

    # Get the rational basis of the null space of Sint
    Nint_mat = nsimplify(Smat, rational=True).nullspace()

    # Convert back to numpy array
    Nint = np.array(Nint_mat).astype(np.float64)

    eps = 1e-9
    Nint[Nint < eps] = 0
    # Remove single reaction loop (reaction involving only cofactors)


def write_matfile(Sint_df, outputfilepath='Sint_no_cofactor_20160831.mat'):
    """Convert S matrix to Matlab sparse matrix (.mat file) for null space analysis.
    Coordinate (start from 1, not 0)

    Args:
        Sint_df (TYPE): The Pandas.DataFrame of the internal S matrix
            without cofactor
        outputfilepath (TYPE): Description

    Returns:
        TYPE: Description
    """
    # convert dataframe to matrix
    Smat = Sint_df.as_matrix()

    # get all indices for non-zero elements in Smat (row, col)
    Smat_nzr, Smat_nzc = np.nonzero(Smat)

    # get all non-zero elements from Smat
    Smat_nze = Smat[Smat_nzr, Smat_nzc]

    # Adjust for matlab coordinate
    Smat_nzr = Smat_nzr + 1
    Smat_nzc = Smat_nzc + 1

    # This final line gives the size of the S matrix in matlab
    nr, nc = Smat.shape

    # Create a 2D array
    sparseMat = np.vstack((Smat_nzr, Smat_nzc, Smat_nze)).T
    sparseMat = np.vstack((sparseMat, np.array([[nr, nc, 0]])))

    # Create a numpy object array from dataframe index
    reactionList = Sint_df.columns.ravel()

    # Write only one matlab .mat file
    scipy.io.savemat(outputfilepath,
                     mdict={'Sint_sparse': sparseMat,
                            'reactionList': np.array(reactionList)}
                     )

    return sparseMat, reactionList


def test_internal_loop_analysis():

    logger = create_logger(
        name="optstoicpy.script.database_preprocessing.test_internal_loop_analysis")

    # Load the base database without exchange reactions
    db = load_base_reaction_db(
        user_defined_export_rxns_Sji=None
    )

    # Load cofactors
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.normpath(
        os.path.join(CURRENT_DIR, '../data/'))
    cofactors_df = pd.read_csv(
        os.path.join(
            DATA_DIR,
            'cofactors_to_exclude.csv'))
    cofactors = cofactors_df['KEGG_ID'].tolist()

    # Remove blocked reaction
    blocked_reactions = json.load(
        open(
            os.path.join(
                DATA_DIR,
                'optstoic_db_v3',
                'optstoic_v3_blocked_reactions_0to5ATP.json'),
            'r+'))
    for rxn in blocked_reactions:
        db.remove_reaction(rxn, refresh_database=False)
    db.refresh_database()

    # Remove cofactors
    S_df = copy.deepcopy(db.S_df)
    S_df_no_cofactor = remove_cofactors_from_Sij(S_df, cofactors)

    assert_equal(S_df_no_cofactor.shape, (1844, 3256))

    # Method 1: MATLAB
    # sparseMat, reactionList = write_matfile(S_df_no_cofactor)
    # run find_null_space.m to obtain all the loops

    # Method 2: This function is not yet implemented in Python as it is too slow
    # internal_loop_analysis(S_df_no_cofactor)

    return S_df_no_cofactor
