import os
import pulp
from optstoicpy.script.utils import create_logger

ORDERED_SOLVERS = ['GUROBI', 'GUROBI_CMD', 'CPLEX_CMD', 'SCIP_CMD', 'GLPK_CMD']

# TODO: Solver parameters need to be updated
SCIP_CMD_PARAMETERS = [
    "limits/gap = 1e-6",
    "limits/absgap = 1e-6",
    "lp/threads = 6",
    "limits/time = 600"]

GUROBI_CMD_OPTIONS = [
    ('Threads', 2),
    ('TimeLimit', 1800),
    ('MIPGapAbs', 1e-6),
    ('MIPGap', 1e-6),
    ('CliqueCuts', 2)]

CPLEX_CMD_OPTIONS = [
    'mip tolerances mipgap 1e-6',
    'mip tolerances absmipgap 1e-6']

GLPK_CMD_OPTIONS = ['--clique', '--pcost', '--gomory', '--mipgap', '1e-6']

SOLVER_KWARGS = {
    'SCIP_CMD': dict(
        solver='SCIP_CMD',
        keepFiles=False,
        mip=True,
        msg=True),
    'GUROBI': dict(
        solver='GUROBI',
        mip=True,
        msg=True,
        timeLimit=1800,
        MIPGapAbs=1e-6),
    'GUROBI_CMD': dict(
        solver='GUROBI_CMD',
        path=None,
        keepFiles=False,
        mip=1,
        msg=1,
        options=GUROBI_CMD_OPTIONS),
    'CPLEX_CMD': dict(
        solver='CPLEX_CMD',
        path=None,
        keepFiles=False,
        mip=1,
        msg=1,
        options=CPLEX_CMD_OPTIONS,
        timelimit=1800),
    'GLPK_CMD': dict(
        solver='GLPK_CMD',
        keepFiles=False,
        msg=1,
        mip=1,
        options=GLPK_CMD_OPTIONS)
}


def load_pulp_solver(
        solver_names=ORDERED_SOLVERS,
        logger=None):
    """Load a pulp solver based on what is available.

    Args:
        solver_name (`list` of `str`, optional): A list of solver names in the order of
            loading preferences.
        logger (None, optional): A logging.Logger object

    Returns:
        pulp.apis.core.LpSolver: A pulp solver instance.
    """
    if logger is None:
        logger = create_logger('optstoic.load_pulp_solver')

    if isinstance(solver_names, str):
        solver_names = [solver_names]

    elif not isinstance(solver_names, list):
        raise Exception("Argument solver_names must be a list!")

    # Load solvers in the order of preferences
    for solver_name in solver_names:
        kwargs = SOLVER_KWARGS.get(solver_name, None)
        pulp_solver = pulp.get_solver(**kwargs)

        if pulp_solver.available():
            logger.warning("Pulp solver set to %s." % solver_name)

            if hasattr(pulp_solver, 'tmpDir'):
                pulp_solver.tmpDir = './'

            if solver_name == 'SCIP_CMD':
                scip_parameter_filepath = create_scip_parameter_file(
                    parameters=SCIP_CMD_PARAMETERS, filepath=pulp_solver.tmpDir)
                pulp_solver.options = [
                    "-s", "{}".format(scip_parameter_filepath)]

            if solver_name == 'GLPK_CMD':
                logger.warning(
                    "GLPK takes a significantly longer time to solve "
                    "OptStoic. Please be patient.")

            return pulp_solver

    logger.warning("No solver is available!")
    return None


def create_scip_parameter_file(parameters=SCIP_CMD_PARAMETERS, filepath="./"):
    """Create a setting file called scip_parameters.set for SCIP_CMD.

    Args:
        parameters (list, optional): A list of parameters. Defaults to SCIP_CMD_PARAMETERS.
        filepath (str, optional): The path to store the scip setting files. Defaults to "./".

    Returns:
        str: The path to the full scip parameters.
    """
    fullfilepath = os.path.join(filepath, "scip_parameters.set")
    with open(fullfilepath, "w+") as parameter_file:
        parameter_file.write("\n".join(parameters))
    return fullfilepath
