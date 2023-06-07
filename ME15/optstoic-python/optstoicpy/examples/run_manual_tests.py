import sys
import os
from optstoicpy.script.utils import create_logger
from optstoicpy.script.database_preprocessing import (
    test_blocked_reactions_analysis,
    test_internal_loop_analysis)
import optstoicpy.script.optstoic_glycolysis as optstoic_gly
import optstoicpy.script.optstoic as optstoic
from optstoicpy.core import drawpathway
from optstoicpy.core.pathway import Pathway, generate_kegg_model

current_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.normpath(os.path.join(current_dir,'../'))
# sys.path.append(parent_dir)


def test_all_optimization_scripts():
    """TODO: To be deprecated or convert to example in Jupyter notebook.
    """
    logger = create_logger(name="optstoicpy.test.testAll")

    logger.info("Test optstoic")
    lp_prob1, pathways1 = optstoic.test_optstoic()

    logger.info("Test optstoic glycolysis. Runtime depends on the solver used.")
    lp_prob2, pathways2 = optstoic_gly.test_optstoic_glycolysis()

    logger.info("Generate kegg_model and draw pathway")
    res_dir = os.path.join(current_dir, 'result/')
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    f = open(os.path.join(res_dir, 'test_KeggModel.txt'), 'w+')

    for ind, p in pathways2.items():
        p.rearrange_reaction_order()
        generate_kegg_model(p, filehandle=f)
        graph_title = "{0}_{1}ATP_P{2}".format(p.name, p.nATP, p.id)
        drawpathway.draw_pathway(p,
                                 os.path.join(res_dir + '/pathway_{0:03d}'.format(p.id)),
                                 imageFormat='png',
                                 graphTitle=graph_title)
    f.close()

    logger.info("optstoicpy.test.testAll completed!")


if __name__ == "__main__":
    test_all_optimization_scripts()
