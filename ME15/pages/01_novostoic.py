from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions as Reactions
#from compound_cacher import CompoundCacher
#from compound import Compound
#from chemaxon import *
#import chemaxon
import streamlit as st
import pandas as pd
import numpy as np
import re
from PIL import Image
import webbrowser
import json
import pickle
import sys
import joblib
import sys
import pulp
import os
import pdb

#sys.path.append('./dGPredictor/CC/')

def count_substructures(radius, molecule):
    """Helper function for get the information of molecular signature of a
    metabolite. The relaxed signature requires the number of each substructure
    to construct a matrix for each molecule.
    Parameters
    ----------
    radius : int
        the radius is bond-distance that defines how many neighbor atoms should
        be considered in a reaction center.
    molecule : Molecule
        a molecule object create by RDkit (e.g. Chem.MolFromInchi(inchi_code)
        or Chem.MolToSmiles(smiles_code))
    Returns
    -------
    dict
        dictionary of molecular signature for a molecule,
        {smiles: molecular_signature}
    """
    m = molecule
    smi_count = dict()
    atomList = [atom for atom in m.GetAtoms()]

    for i in range(len(atomList)):
        env = Chem.FindAtomEnvironmentOfRadiusN(m, radius, i)
        atoms = set()
        for bidx in env:
            atoms.add(m.GetBondWithIdx(bidx).GetBeginAtomIdx())
            atoms.add(m.GetBondWithIdx(bidx).GetEndAtomIdx())

        # only one atom is in this environment, such as O in H2O
        if len(atoms) == 0:
            atoms = {i}

        smi = Chem.MolFragmentToSmiles(m, atomsToUse=list(atoms),
                                       bondsToUse=env, canonical=True)

        if smi in smi_count:
            smi_count[smi] = smi_count[smi] + 1
        else:
            smi_count[smi] = 1
    return smi_count


def parse_reaction_formula_side(s):
    """
        Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        Ignores stoichiometry.

        Returns:
            The set of CIDs.
    """
    if s.strip() == "null":
        return {}

    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if len(tokens) == 0:
            continue
        if len(tokens) == 1:
            amount = 1
            key = member
        else:
            amount = float(tokens[0])
            key = tokens[1]

        compound_bag[key] = compound_bag.get(key, 0) + amount

    return compound_bag


def parse_formula(formula, arrow='<=>', rid=None):
    """
        Parses a two-sided formula such as: 2 C00001 => C00002 + C00003

        Return:
            The set of substrates, products and the direction of the reaction
    """
    tokens = formula.split(arrow)
    if len(tokens) < 2:
        print(('Reaction does not contain the arrow sign (%s): %s'
               % (arrow, formula)))
    if len(tokens) > 2:
        print(('Reaction contains more than one arrow sign (%s): %s'
               % (arrow, formula)))

    left = tokens[0].strip()
    right = tokens[1].strip()

    sparse_reaction = {}
    for cid, count in parse_reaction_formula_side(left).items():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count

    for cid, count in parse_reaction_formula_side(right).items():
        sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count

    return sparse_reaction


def draw_rxn_figure(rxn_dict, db_smiles, novel_smiles):
    # db_smiles = load_smiles()

    left = ''
    right = ''

    for met, stoic in rxn_dict.items():
        if met == "C00080" or met == "C00282":
            continue  # hydogen is not considered
        if stoic > 0:
            if met in db_smiles:
                right = right + db_smiles[met] + '.'
            else:
                right = right + novel_smiles[met] + '.'
        else:
            if met in db_smiles:
                left = left + db_smiles[met] + '.'
            else:
                left = left + novel_smiles[met] + '.'
    smarts = left[:-1] + '>>' + right[:-1]
    # print smarts
    smarts = str(smarts)
    rxn = Reactions.ReactionFromSmarts(smarts, useSmiles=True)
    return Draw.ReactionToImage(rxn)  # , subImgSize=(400, 400))


def novoStoic_minFlux_relaxedRule(exchange_mets, novel_mets, project, iterations, n_steps, pulp_solver, use_direction):
    """apply reaction rules generated from a more relaxed manner to search for
    reaction rules that are able to fill the gap between the source and sink
    metabolites.
    - rePrime procedure is more similar to a morgan fingerprints
    - the relaxed rule is generated from substructures without considering the
      bond that connect the atoms at the edge of the substructure to the rest
      of the molecules

    Parameters
    ----------
    exchange_mets : dict
        overall stoichiometry of source and sink metabolites, {met: stoic,...}
        This is a important input for novoStoic to run correctly because the
        method requires that overall moieties are balanced.
    novel_mets : list
        list of novel metabolites that are not in the database (novoStoic/data/
        metanetx_universal_model_kegg_metacyc_rhea_seed_reactome.json)
    filtered_rules : list
        list of rules that are filtered by the user (based on expert knowldedge)
        to reduce the running time of the novoStoic search process
    project : string
        a path to store the tmp information of result from running novoStoic
    iterations : int
        the number of iterations of searching for alternative solutions
    data_dir : type
        Description of parameter `data_dir`.

    Returns
    -------
    None
        all the outputs are saved in the project folder.

    """
    if not os.path.exists(project):
        os.makedirs(project)

    # the maximum flux of a reaction
    M = 2

    data_dir = './dGPredictor/data'

    # read csv files with molecular signatures and reaction rules
    molecular_signature = json.load(open(
        os.path.join(data_dir, 'decompose_vector_ac.json')))
    molsigs = pd.DataFrame.from_dict(molecular_signature).fillna(0)

    rules = pd.read_csv(
        os.path.join(data_dir, "relaxed_rule_noduplic.csv"), index_col=0
    )

    ###### sets ############
    moiety_index = rules.index.tolist()  # moiety sets
    rules_index = rules.columns.values.tolist()
    st.write("Number of unique rules used in this search:", len(rules_index))

    exchange_index = exchange_mets.keys()

    ###### parameters ######
    # T(m,r) contains atom stoichiometry for each rule
    T = rules.to_dict(orient="index")

    # C(m,i) contains moiety cardinality for each metabolite
    C = molsigs.to_dict(orient="index")
    for m in moiety_index:
        C[m]["C00080"] = 0
        C[m]["C00282"] = 0

    # add metabolites that are not present in current database
    for met in novel_mets:
        smiles = novel_mets[met]
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.RemoveHs(mol)
        molsigs_product_dict = count_substructures(1, mol)

        for m in moiety_index:
            if m in molsigs_product_dict.keys():
                C[m][met] = molsigs_product_dict[m]
            else:
                C[m][met] = 0

    ###### variables ######
    v_rule = pulp.LpVariable.dicts(
        "v_rule", rules_index, lowBound=-M, upBound=M, cat="Integer"
    )
    v_rule_obj = pulp.LpVariable.dicts(
        "v_rule_obj", rules_index, lowBound=0, upBound=M, cat="Continuous"
    )

    v_EX = pulp.LpVariable.dicts(
        "v_EX", exchange_index, lowBound=-M, upBound=M, cat="Continuous"
    )
    y_rule = pulp.LpVariable.dicts(
        "y", rules_index, lowBound=0, upBound=1, cat="Binary"
    )

    # create MILP problem
    lp_prob = pulp.LpProblem("novoStoic", pulp.LpMinimize)

    ####### objective function ####
    lp_prob += pulp.lpSum([v_rule_obj[j] for j in rules_index])

    ####### constraints ####
    # constraint 1: moiety change balance
    for m in moiety_index:
        lp_prob += (
            pulp.lpSum([T[m][r] * v_rule[r]
                       for r in rules_index if T[m][r] != 0])
            == pulp.lpSum([C[m][i] * v_EX[i] for i in exchange_index if C[m][i] != 0]),
            "moiety_balance_" + str(moiety_index.index(m)),
        )

    # constraint 2: constraint for exchange reactions
    for i, stoic in exchange_mets.items():
        lp_prob += v_EX[i] == stoic, "exchange" + i

    # constraint 3: control the number of rules

    direction_df = pd.read_csv(
        os.path.join(data_dir, "direction.csv"), index_col=0
    )
    direction_df.index = direction_df['reaction']

    # direction: 0-reversible, 1-backward, 2-forward
    direction = direction_df['direction'].to_dict()

    if use_direction:
        soln_file = os.path.join(project, "solution_use_direction.txt")
        for j in rules_index:
            if direction[j] == 0:
                lp_prob += v_rule[j] >= y_rule[j] * -M, "cons1_%s" % j
                lp_prob += v_rule[j] <= y_rule[j] * M, "cons2_%s" % j
            if direction[j] == 1:
                lp_prob += v_rule[j] >= y_rule[j] * -M, "cons1_%s" % j
                lp_prob += v_rule[j] <= 0, "cons2_%s" % j
            if direction[j] == 2:
                lp_prob += v_rule[j] >= 0, "cons1_%s" % j
                lp_prob += v_rule[j] <= y_rule[j] * M, "cons2_%s" % j
    else:
        soln_file = os.path.join(project, "solution_no_direction.txt")
        for j in rules_index:
            lp_prob += v_rule[j] >= y_rule[j] * -M, "cons1_%s" % j
            lp_prob += v_rule[j] <= y_rule[j] * M, "cons2_%s" % j

    for j in rules_index:
        lp_prob += v_rule_obj[j] >= v_rule[j]
        lp_prob += v_rule_obj[j] >= -v_rule[j]

    # constraint 5: customized constraints
    # the number of steps of the pathway
    lp_prob += pulp.lpSum([v_rule_obj[j] for j in rules_index]) == n_steps

    # solve
    integer_cuts(lp_prob, pulp_solver, iterations, rules_index,
                 y_rule, v_rule, soln_file, direction)


def integer_cuts(lp_prob, pulp_solver, iterations, rules_index, y_rule, v_rule, soln_file, direction):
    """add integer cut constraints to a mixed-integer linear programming problem
    (MILP). The aim of such constraints is to find alternative solutions by
    adding constraints to exclude the already explored solutions.

    Reference: Optimization Methods in Metabolic Networks By Costas D. Maranas,
    Ali R. Zomorrodi, Chapter 4.2.2 Finding alternative optimal integer
    solutions

    Returns
    -------
    type
        Description of returned object.

    """
    for sol_num in range(1, iterations + 1):
        integer_cut_rules = []

        # optinal output: lp file for debug
        lp_prob.writeLP('./test.lp')
        # if pulp_solver = "SCIP":
        # status, values = pulp_solver.solve(lp_prob)
        lp_prob.solve(pulp_solver)
        # pulp_solver.solve(lp_prob)

        st.write("Status:", pulp.LpStatus[lp_prob.status])

        if pulp.LpStatus[lp_prob.status] != 'Optimal':
            break

        st.write('----------- Reaction Rules  --------------')
        write_str_sol_num = '------     iterations: ' + str(sol_num) + '   -----'
        st.write(write_str_sol_num)

        with open(soln_file, 'a') as f:
            f.write('iteration,' + str(sol_num))
            f.write('\n')

        for r in rules_index:
            if (v_rule[r].varValue >= 0.1 or v_rule[r].varValue <= -0.1):

                dG_info = ''
                if (v_rule[r].varValue > 0 and direction[r] == 1) or (v_rule[r].varValue < 0 and direction[r] == 2):
                    dG_info = ' * Thermodynamically infeasible'
                    st.write("##### Found ####: " + str(r) + dG_info)
                integer_cut_rules.append(r)
                #st.write(r, v_rule[r].varValue)

                st.write(r + ',' + str(v_rule[r].varValue) + dG_info + '\n')

                with open(soln_file, 'a') as f:
                    f.write(r + ',' + str(v_rule[r].varValue) + dG_info)
                    f.write('\n')

        length = len(integer_cut_rules) - 1
        lp_prob += (
            pulp.lpSum([y_rule[r] for r in integer_cut_rules]) <= length,
            "integer_cut_" + str(sol_num),
        )


def main():

    st.image('./novoka/figures/novoStoic_header.png')
    st.subheader('Overall Stoichiometric Equation')
    stoic = st.text_input('Enter the pathway stoichiometry here',
                          value='1 C00141 + 1 C00004 <=> 1 C00003 + 1 C14710 + 1 C00011')

    if st.checkbox('Reaction has metabolites not in database'):
        add_info = st.text_area('Additional information (id: SMILES):',
                                '{"14bdo":"OCCCCO"}')
    else:
        add_info = {}

        novel_metab = add_info
    col1, col2 = st.columns(2)

    with col1:
        #st.subheader('Primary Product')
        p_prod = st.text_input('Enter primary product: ', value='C14710')

        with col2:
            p_subs = st.text_input('Enter primary substrate: ', value='C00141')

        st.subheader('Pathway design parameters:')
    max_steps = st.slider('Maximum number of steps',
                          min_value=1.0, max_value=10.0, value=2.0, step=1.0)
    iterations = st.slider('Max pathway search', min_value=0.0,
                           max_value=15.0, value=5.0, step=1.0)

    if st.button("Search"):
        with st.spinner('Searching...'):
            rxn_dict = parse_formula(stoic)
        st.subheader('Reaction lists in pathway design')
        save_sol_folder = './novoka/novoStoic_solutions/' + str(p_prod)
        pulp_solver = pulp.CPLEX_CMD(
            path=None, keepFiles=0, mip=1, msg=1)
        use_direction = False
        novoStoic_minFlux_relaxedRule(
            rxn_dict, novel_metab, save_sol_folder, np.int(iterations), np.int(max_steps), pulp_solver, use_direction)


if __name__ == '__main__':
    main()
