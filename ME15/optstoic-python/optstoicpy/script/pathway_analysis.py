"""Analyze similarity index between different set of pathways.
Compare and combine pathways from different parallel simulation (e.g. python and gams).
"""
from __future__ import print_function
import os
import copy
import pickle as pickle
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from optstoicpy.core.drawpathway import *
from optstoicpy.core.pathway import Pathway, generate_kegg_model
import matplotlib
from builtins import range
from builtins import zip
from future import standard_library
standard_library.install_aliases()
# matplotlib.use('SVG')
matplotlib.use('Agg')


def get_pathway_identity_matrix(pathwayObjList, symmetry=True):
    """
    Matrix[i][j] = 1 if pathway i and pathway j are identical, otherwise 0.

    Arguments:
        pathwayObjList (TYPE): A list of pathway objects
        symmetry (bool, optional): If True, return a full matrix. If False, return only upper triangular matrix.

    Returns:
        TYPE: Description

    """
    nr = len(pathwayObjList)

    # calculate identity score between the two different set of pathways
    similarity_mat = np.identity(nr)
    for i in range(nr):
        for j in range(i + 1, nr):
            similarity_mat[i][j] = pathwayObjList[i].is_same_pathway_with(
                pathwayObjList[j])

    # Fill the lower triangular matrix with the same value
    if symmetry:
        similarity_mat = similarity_mat + similarity_mat.T - \
            np.diag(similarity_mat.diagonal())

    return similarity_mat


def get_unique_pathways_from_list(
        pathwayObjList,
        update_unique_id=True,
        sort_by_total_flux=False,
        debug=False):
    """
    return a set of unique pathways (obj) from a list of Pathway objects

    Arguments:
        pathwayObjList (TYPE): list of pathway objects
        update_unique_id (bool, optional): change the ID to unique id (default True)
        sort_by_total_flux (bool, optional): sort the unique pathways by total flux (default False)
        debug (bool, optional): output a data index to unique pathway id mapping dictionary (default False)

        Note: Currently, "debug=True" does not work with "sort_by_total_flux=True"

    Returns:
        TYPE: Description
    """

    similarity_mat = get_pathway_identity_matrix(
        pathwayObjList, symmetry=False)
    #fig = plot_single_similarity_matrix(similarity_mat, nATP=1)
    # fig.savefig('before_sim.png')
    score = np.sum(similarity_mat, axis=0)

    unique_index = np.where(score == 1)[0]

    all_unique_pathways = list(np.array(pathwayObjList)[unique_index])

    if sort_by_total_flux:
        newlist = sorted(
            all_unique_pathways,
            key=lambda p: p.total_flux_no_exchange,
            reverse=False)
        all_unique_pathways = newlist

    if update_unique_id:
        for i, p in enumerate(all_unique_pathways, start=1):
            p.id = i

    print(len(all_unique_pathways))
    #unique_pathway_similarity_mat = get_pathway_identity_matrix(all_unique_pathways, symmetry=True)
    #fig = plot_single_similarity_matrix(unique_pathway_similarity_mat, nATP=1)
    # fig.savefig('after_sim.png')

    # ---debug function: to check which gams code generate feasible result
    if debug:

        map_unique_index_to_pathid = dict(
            list(zip(unique_index, list(range(1, len(unique_index) + 1)))))
        data_to_id_map = copy.deepcopy(map_unique_index_to_pathid)

        for ind in np.where(score > 1)[0]:
            # get the first occurence of the pathway
            first_occur = np.argmax(similarity_mat[:, ind])
            # get the final id of that
            data_to_id_map[ind] = map_unique_index_to_pathid[first_occur]
        return all_unique_pathways, data_to_id_map
    else:
        return all_unique_pathways


def find_identical_pathways_and_get_unique_pathways(pres, gres):
    """
    Return a similarity matrix between Pathway objects in pres and gres,
    and a set of unique pathways (obj).
    This works only when pres and gres are already unique sets of pathways.
    Arguments:
        pres -- pathway object files generated from python
        gres -- pathway object files generated from gams

    """
    # gres[0].is_same_pathway_with(gres[1])
    lpy = len(pres)
    lgams = len(gres)

    # calculate identity score between the two different set of pathways
    similarity_mat = np.zeros((lpy, lgams))
    for i, p1 in enumerate(pres):
        for j, p2 in enumerate(gres):
            similarity_mat[i][j] = p1.is_same_pathway_with(p2)

    # set of identical pathways (intersection between two sets)
    pmatch, gmatch = np.where(similarity_mat == 1)

    assert len(pmatch) == len(set(pmatch))

    # set of unique pathways
    unique_pathways = [pres[i] for i in range(0, len(pres)) if i not in pmatch]
    all_unique_pathways = copy.deepcopy(gres)

    for i, p in enumerate(unique_pathways, start=1):
        p.id = len(gres) + i
        all_unique_pathways.append(p)

    print(len(all_unique_pathways))

    return similarity_mat, all_unique_pathways


def calculate_jaccard_score_between_pathways(pathway_set):
    """
    Calculate identity score between the same set of pathways
    Arguments:
    pathway_set  --  a list of pathway objects
    """
    nr = len(pathway_set)
    mat = np.identity(nr)

    for i in range(nr):
        for j in range(i + 1, nr):
            mat[i][j] = pathway_set[i].get_pathway_similarity_index(
                pathway_set[j])
            mat[j][i] = mat[i][j]
    return mat


def make_dir_if_not_exist(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
    return 1


def plot_similarity_matrix(similarity_mat, similarity_mat_all_pathways):

    xmax1, ymin1 = np.shape(similarity_mat)

    fig = plt.figure(figsize=(18, 8))
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    im1 = ax1.imshow(
        similarity_mat,
        interpolation='none',
        extent=[
            0,
            xmax1,
            ymin1,
            0],
        vmin=0,
        vmax=1,
        aspect='auto',
        cmap=cm.coolwarm)
    cbar = plt.colorbar(im1)
    cbar.set_label('Jaccard Index')
    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_label_position('top')

    ax1.set_xlabel('Python pathway #')
    ax1.set_ylabel('GAMS pathway #')

    plt.text(
        0.5,
        1.08,
        'Python vs gams',
        horizontalalignment='center',
        fontsize=12,
        transform=ax1.transAxes)
    # ax1.set_title()
    ax1.grid(True)
    #fig.savefig('compare_python_gams_1ATP_pathways', transparent=True)

    xmax2, ymin2 = np.shape(similarity_mat_all_pathways)
    ax2 = plt.subplot2grid((1, 2), (0, 1))
    im2 = ax2.imshow(
        similarity_mat_all_pathways,
        interpolation='none',
        extent=[
            0,
            xmax2,
            ymin2,
            0],
        aspect='auto')
    cbar = plt.colorbar(im2)
    cbar.set_label('Jaccard Index')
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.set_xlabel('Pathway #')
    ax2.set_ylabel('Pathway #')
    #ax2.set_title('All unique 1ATP pathways')
    plt.text(
        0.5,
        1.08,
        'Similarity index for all unique pathways',
        horizontalalignment='center',
        fontsize=12,
        transform=ax2.transAxes)

    return fig


def plot_single_similarity_matrix(similarity_mat_all_pathways, nATP=1):

    fig, ax2 = plt.subplots(figsize=(10, 10))
    xmax2, ymin2 = np.shape(similarity_mat_all_pathways)
    im2 = ax2.imshow(
        similarity_mat_all_pathways,
        interpolation='none',
        extent=[
            0,
            xmax2,
            ymin2,
            0],
        aspect='auto')
    cbar = plt.colorbar(im2)
    cbar.set_label('Jaccard Index')
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.set_xlabel('Pathway #')
    ax2.set_ylabel('Pathway #')
    #ax2.set_title('All unique 1ATP pathways')
    plt.text(0.5, 1.08, 'Similarity index for all unique {0} ATP pathways'.format(
        nATP), horizontalalignment='center', fontsize=12, transform=ax2.transAxes)

    return fig


def draw_all_pathways(pathway_set, outputFilePath, cutoff=0):
    """Draw all the pathway given a list of pathway objects
    Arguments:
    pathway_set  --  a list of pathway objects
    outputFilePath -- output file path
    cutoff -- id cutoff for drawing pathway
    """
    print("Drawing all pathways. Be patient...")
    for p in pathway_set:
        if p.id > cutoff:
            graph_title = "Final_{0}_P{1}".format(p.name, p.id)
            draw_pathway(p,
                         os.path.join(outputFilePath,
                                      '/pathway_{0:03d}'.format(p.id)),
                         imageFormat='png',
                         graphTitle=graph_title,
                         cleanup=True,
                         darkBackgroundMode=False)
    print("Done!")

    return 1


def draw_selected_pathways(pathway_set, outputFilePath, selected_ids=[],
                           file_prefix='selected_pathway_',
                           imageFormat='png', darkBackgroundMode=False):
    """Draw all the pathway given a list of pathway objects
    Arguments:
    pathway_set  --  a list of pathway objects
    outputFilePath -- output file path
    cutoff - id cutoff for drawing pathway
    """
    print("Drawing all selected pathways. Be patient...")
    for p in pathway_set:
        if p.id in selected_ids:
            graph_title = "Final_{0}_P{1}".format(p.name, p.id)
            draw_pathway(
                p,
                os.path.join(
                    outputFilePath,
                    file_prefix +
                    '{0:03d}'.format(
                        p.id)),
                imageFormat=imageFormat,
                graphTitle=graph_title,
                cleanup=True,
                darkBackgroundMode=darkBackgroundMode)
    print("Done!")

    return 1


def extract_pathway_set(pathway_set, selected_ids=[]):
    new_set = []
    for p in pathway_set:
        if p.id in selected_ids:
            new_set.append(p)
    return new_set


def combine_multiple_pathways(pathway_set):

    combined_reactions = {}

    for pathway in pathway_set:

        for rxn in pathway.reactions:
            if rxn.rid not in combined_reactions:
                combined_reactions[rxn.rid] = {}
                combined_reactions[rxn.rid]['count_f'] = 0
                combined_reactions[rxn.rid]['count_r'] = 0
                combined_reactions[rxn.rid]['example_f'] = None
                combined_reactions[rxn.rid]['example_r'] = None

            if rxn.flux > 0:
                combined_reactions[rxn.rid]['count_f'] += 1
                if combined_reactions[rxn.rid]['example_f'] is None:
                    # normalize the flux
                    combined_reactions[rxn.rid]['example_f'] = copy.deepcopy(
                        rxn)
            else:
                combined_reactions[rxn.rid]['count_r'] += 1
                if combined_reactions[rxn.rid]['example_r'] is None:
                    combined_reactions[rxn.rid]['example_r'] = copy.deepcopy(
                        rxn)

    return combined_reactions


def draw_combined_pathway(combined_reactions, fileName):
    reactionObjList = []
    for r, v in combined_reactions.items():
        if v['count_f'] > 0:
            # Create a reaction object with number of reaction occurence in all
            # pathways as flux
            v['example_f'].flux = v['count_f']
            reactionObjList.append(v['example_f'])

        if v['count_r'] > 0:
            v['example_r'].flux = -1 * v['count_r']
            reactionObjList.append(v['example_r'])

    p = Pathway(id='', name='Combined_pathway', reactions=reactionObjList)

    draw_pathway(
        p,
        fileName,
        imageFormat='png',
        graphTitle=fileName,
        scaleLineWidth=True,
        scalingFactor=10.0,
        cleanup=True,
        engine='dot')
    # Layout engines: circo dot fdp neato nop1 nop2 osage patchwork sfdp twopi
    return 1


if __name__ == '__main__':
    pass

    nATP = 1
    outputpath = './{0}ATP_final'.format(nATP)

    #outputpath = './{0}ATP_final_newconstraints'.format(nATP)
    """
    make_dir_if_not_exist(outputpath)

    if float(nATP).is_integer():
        nATP_str = str(nATP)
    else:
        nATP_str = '_'.join(str(nATP).split('.'))

    f1= '{0}ATP/OptStoic_{0}ATP_pathways_obj.pkl'.format(nATP)
    f2 = '{0}ATP_gams/OptStoic_gams_{0}ATP_pathways_obj.pkl'.format(nATP)

    pres = pickle.load(open(f1,'r'))
    gres = pickle.load(open(f2,'r'))

    similarity_mat, all_unique_pathways = find_identical_pathways_and_get_unique_pathways(pres, gres)

    for p in all_unique_pathways:
        p.name = "OptStoic_glycolysis_{0}ATP".format(nATP)

    pickle.dump(all_unique_pathways, open(os.path.join(outputpath, '{0}ATP_all_pathways.pkl'.format(nATP)), 'w'))

    # #Write unique pathways to Kegg Model
    f = open(os.path.join(outputpath, '{0}ATP_{1}_KeggModel_ratio.txt'.format(nATP, len(all_unique_pathways))), 'w+')
    for p in all_unique_pathways:
        generate_kegg_model(p, filehandle=f, add_ratio_constraints=True)
    f.close()

    similarity_mat_all_pathways = calculate_jaccard_score_between_pathways(all_unique_pathways)

    fig = plot_similarity_matrix(similarity_mat, similarity_mat_all_pathways)
    fig.savefig(os.path.join(outputpath,'compare_{0}ATP_pathways.png'.format(nATP_str)), transparent=True)

    fig =  plot_single_similarity_matrix(similarity_mat_all_pathways, nATP)
    fig.savefig(os.path.join(outputpath,'Similarity_plot_{0}ATP_pathways.png'.format(nATP_str)), transparent=True)

    #draw_all_pathways(all_unique_pathways, outputpath, cutoff=600)
    """
