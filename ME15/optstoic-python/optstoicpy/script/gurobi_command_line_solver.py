# /usr/bin/python
"""
Using Gurobi command line (gurobi_cl) to solve lp problem generated with Pulp.
This is written to solve the issue of unable to install Gurobi in Cygwin.
"""
from __future__ import print_function
import os
import pdb
import sys
from subprocess import Popen, PIPE


def solve_with_gurobi_cl_debug(
        lp_filename,
        options='Threads=2 TimeLimit=1200 MIPGapAbs=1e-6'):
    args = 'gurobi_cl ' + options + \
        ' ResultFile=%s.sol %s.lp' % (lp_filename, lp_filename)
    command = args.split()
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    output = ""
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        output += nextline
        sys.stdout.write(nextline)
        sys.stdout.flush()

    # Iteration limit reached
    if 'Optimal solution found' in output:
        return 'Optimal', None
    elif 'Time limit reached' in output:
        status_message = [
            l for l in output.splitlines() if l.startswith('Best')]
        return 'Time_limit', status_message[0]
    else:
        return 'Not_optimal', None


def solve_with_gurobi_cl(
        lp_filename,
        options='Threads=2 TimeLimit=1200 MIPGapAbs=1e-6',
        verbose=True):
    args = 'gurobi_cl ' + options + \
        ' ResultFile=%s.sol %s.lp' % (lp_filename, lp_filename)
    command = args.split()
    process = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    #exitCode = process.returncode

    if verbose:
        print(stdout)

    if stderr:
        print(stderr)

    # Iteration limit reached
    if 'Optimal solution found' in stdout:
        return 'Optimal', None
    elif 'Time limit reached' in stdout:
        status_message = [
            l for l in stdout.splitlines() if l.startswith('Best')]
        return 'Time_limit', status_message[0]
    else:
        return 'Not_optimal', None


def parse_gurobi_sol(sol_filename):
    try:
        f = open(sol_filename + '.sol', 'rU')
    except IOError:
        print("%s.sol is not in the current directory." % sol_filename)
        return None

    data = f.read().splitlines()
    # Objective function
    objective_function = float(data[0].split()[-1])
    varValueList = [line.split() for line in data[1:] if line is not None]
    varValue = dict((k, float(v)) for (k, v) in varValueList)
    return objective_function, varValue


if __name__ == "__main__":
    solve_with_gurobi_cl('OptStoic')
    objective_function, varValue = parse_gurobi_sol('OptStoic')
