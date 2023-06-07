from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from optstoicpy.core.reaction import Reaction
from optstoicpy.core.pathway import Pathway, generate_kegg_model
from optstoicpy.core.drawpathway import *
import pickle as pickle
import os
import json
import pdb
import logging


def fix_incomplete_json(inputfilename, outputfilename):
    """
    This fix the json file if gams are terminated prior to
    completion of an optstoic iteration.
    """

    f = open(inputfilename, 'rU')
    temp = f.readlines()
    f.close()
    last_complete_entry_index = 0

    if temp[-1] not in ['}\n', '}']:
        #find the last complete iteration
        for ind, line in enumerate(temp[::-1]):
            if line.startswith('   },\n'):
                last_complete_entry_index = ind
                break

        if last_complete_entry_index == 0:
        #if unable to find any of them
            print("There is no complete entry in the json file. Quit program.")
            return None

        newtemp = temp[0: -last_complete_entry_index]
        newtemp[-1] = '   }\n'
        newtemp.append('}\n')

        output = open(outputfilename, 'w+')
        output.writelines(newtemp)
        output.close()
        return newtemp


    else:
        #if the second last line contain trailing comma
        if temp[-2] == '   },\n':
            temp[-2] = '   }\n'
            output = open(outputfilename, 'w+')
            output.writelines(temp)
            output.close()
        else:

            print("The JSON file is valid.")
        return temp

def make_integer_cut(resultDict, outputfname):
    """ This generates gams readable file for integer cut using the resultSet data"""

    fid = open(outputfname, 'w+')
    for ind, res in sorted(resultDict.items()):
        if "pathway" in res:
            if "num_reaction" not in res:
                #if pathway is incomplete
                continue
            for rid in sorted(res['pathway'].keys()):
                fid.write("'%d'.'%s' 1\n"%(ind, rid))
            fid.write("\n")
    fid.close()
    return 1


def runAnalysis(resultDict, numATP, outputFilePath, imgFormat='png', shift_pathway_id_by=0,\
                sourceSubstrateID='C00031', endSubstrateID='C00022', darkBackgroundMode=False):
    logging.info("Analyzing results... \n")

    pathway_objects = []

    outputFileName='OptStoic_gams_{0}ATP'.format(numATP)

    f = open(os.path.join(outputFilePath, outputFileName + '_KeggModel.txt'), 'w+')

    for ind, res in sorted(resultDict.items()):
        logging.debug("Pathway %s"%ind)
        if "pathway" not in res:
            logging.info("OptStoic terminated with infeasible solution.")
            continue
        elif "num_reaction" not in res:
            logging.info("Pathway is incomplete...")
            continue

        p = Pathway(id=int(ind)+shift_pathway_id_by, name='OptStoic_gams', reaction_ids=list(res['pathway'].keys()),
                    fluxes=list(res['pathway'].values()), sourceSubstrateID=sourceSubstrateID, endSubstrateID=endSubstrateID,
                    total_flux_no_exchange=res['total_flux_no_exchange'],
                    note={'modelstat': res.get("modelstat"), 'solvestat': res.get("solvestat")})
        p.rearrange_reaction_order()
        pathway_objects.append(p)
        generate_kegg_model(p, filehandle=f)

        graph_title = "{0}_{1}ATP_P{2}".format(p.name, p.nATP, p.id)
        if res['modelstat'] != 1:
            graph_title += '; Modelstat={0}'.format(res['modelstat'])
        if imgFormat:
            draw_pathway(p, os.path.join(outputFilePath+'/pathway_{0:03d}'.format(p.id)),
                        imageFormat=imgFormat, graphTitle=graph_title, darkBackgroundMode=darkBackgroundMode)


    #pickle.dump(resultDict, open(outputFilePath+outputFileName+'_pathways_dict.pkl', 'w+'))
    pickle.dump(pathway_objects, open(outputFilePath + outputFileName + '_pathways_obj.pkl', 'w+'))

    logging.info("\nDone!\n")
    logging.info("Check your output folder: %s"%outputFilePath)

    return pathway_objects




if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    """
    numATP = 2.5
    fname = 'loopless_minflux_glycolysis_2_5ATP_20160703_zlb.json'

    fix_incomplete_json(fname, fname)

    data = json.load(open(fname, 'r'))

    data2 = {int(k): v for k, v in data.items()}

    make_integer_cut(data2, 'integer_cut/20160703_2_5ATP_integer_cut.txt')

    # outputFilePath = 'result/{0}ATP/'.format(numATP)

    # if not os.path.exists(outputFilePath):
    #     os.makedirs(outputFilePath)
    """

    #runAnalysis(data2, numATP, outputFilePath, imgFormat='png')

    ###-----------------------
    ###Compare pathways with existing pathways
    all_pathways = pickle.load(open('result/1ATPOptStoic_gams_1ATP_pathways_obj.pkl', 'r'))

    ED = {
            'R00200': -1.00000000,
            'R00299': 1.00000000,
            'R00658': 1.00000000,
            'R00835': 1.00000000,
            'R01061': 1.00000000,
            'R01512': -1.00000000,
            'R01518': -1.00000000,
            'R02035': 1.00000000,
            'R02036': 1.00000000,
            'R05605': 1.00000000,
            'EX_glc': -1.00000000,
            'EX_nad': -1.00000000,
            'EX_adp': -1.00000000,
            'EX_phosphate'  :  -1.00000000,
            'EX_pyruvate' :   2.00000000,
            'EX_nadh' :   1.00000000,
            'EX_atp' :   1.00000000,
            'EX_h2o' :   1.00000000,
            'EX_h+':   2.00000000,
            'EX_nadp'  :  -1.00000000,
            'EX_nadph' :   1.00000000,
            }

    EMP = {
            'R00200': -2.00000000,
            'R00299': 1.00000000,
            'R00658': 2.00000000,
            'R00756': 1.00000000,
            'R00771': 1.00000000,
            'R01015': -1.00000000,
            'R01061': 2.00000000,
            'R01068': 1.00000000,
            'R01512': -2.00000000,
            'R01518': -2.00000000,
            'EX_glc': -1.00000000,
            'EX_nad': -2.00000000,
            'EX_adp': -2.00000000,
            'EX_phosphate':-2.00000000,
            'EX_pyruvate': 2.00000000,
            'EX_nadh': 2.00000000,
            'EX_atp': 2.00000000,
            'EX_h2o': 2.00000000,
            'EX_h+': 2.00000000
    }

    EDpath = Pathway(id='ED', name='ED', reaction_ids=list(ED.keys()), fluxes=list(ED.values()))
    EDpath.rearrange_reaction_order()
    res = []
    for cpath in all_pathways:
        if EDpath.is_same_pathway_with(cpath):
            print(cpath)
        res.append(EDpath.get_pathway_similarity_index(cpath))