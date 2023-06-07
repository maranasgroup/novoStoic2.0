from __future__ import division
from builtins import str
from past.utils import old_div
from .pathway import Pathway
from .config import cofactorsList, kegg_compound, color_configs
import graphviz as gv
import os
import logging
import math

# ##################CONSTANTS######################
kegg_compound['C00009'] = 'Pi'
kegg_compound['C00013'] = 'PPi'
kegg_compound['C00236'] = '1,3-Bisphospho-D-glycerate'
kegg_compound['C00111'] = 'dihydroxyacetone phosphate'

REACTION_FONT_SIZE = '20'

# ########################################


def load_global_styles(colorConfig):
    """
    Create a global styles dictionary for all Graphviz graph
    """
    colorMapping = colorConfig['colorMapping']
    for c in cofactorsList:
        if c not in colorMapping:
            colorMapping[c] = colorConfig['OTHER_COFACTOR_COLOR']

    global_styles = {
        'graph': {
            'fontsize': '20',
            'fontname': 'Helvetica',
            'bgcolor': colorConfig['BACKGROUND_COLOR'],
            # 'rankdir': 'BT',
        },
        'nodes': {
            'fontname': 'Helvetica',
            'fontsize': '25',
            # 'shape': 'hexagon',
            'fontcolor': colorConfig['ALL_FONT_COLOR'],
            # 'color': 'white',
            # 'style': 'filled',
            # 'fillcolor': '#006699',
        },
        'edges': {
            # 'style': 'dashed',
            # 'color': 'white',
            # 'arrowhead': 'open',
            'fontname': 'Helvetica',
            'fontsize': '14',
            # 'fontcolor': 'white',
        }
    }
    return global_styles, colorMapping


def apply_styles(graph, styles):
    """
    Apply the styles to a Graphviz graph.
    """
    graph.graph_attr.update(
        ('graph' in styles and styles['graph']) or {}
    )
    graph.node_attr.update(
        ('nodes' in styles and styles['nodes']) or {}
    )
    graph.edge_attr.update(
        ('edges' in styles and styles['edges']) or {}
    )
    return graph


def draw_pathway(
        Pathway,
        imageFileName=None,
        imageFormat='png',
        graphTitle='',
        scaleLineWidth=False,
        scalingFactor=200.0,
        cleanup=True,
        engine='dot',
        darkBackgroundMode=False,
        width=5,
        height=5):
    """
    Draw a digraph for a Pathway objects and render it as
    the given imageFormat using Graphviz.

    Args:
        Pathway (TYPE): A Pathway object (pathway.py)
        imageFileName (None, optional): Name of the output file (default Pathway.name)
        imageFormat (str, optional): Any format that Graphviz can support (default 'png')
        graphTitle (str, optional): Title of the output graph
        scaleLineWidth (bool, optional): If true, scale the penwidth of an edge
            to a value between 1 and 10. This is useful when
            fluxes are too large.
            Else, the penwidth of an edge is absolute value
            of the flux value.  (default False)
        scalingFactor (float, optional): If scaleLineWidth is true,
            penwidth = (abs(flux)/scalingFactor) * 10 + 1.
            (E.g. Use the maximum flux values of a
            pathway as the scaling Factor)
        cleanup (bool, optional): delete the ".dot" file after drawing
        engine (str, optional): Graphviz layout engine used to render the graph.
            Layout engines = {'circo', 'dot', 'fdp', 'neato', 'nop1', 'nop2',
                             'osage', 'patchwork', 'sfdp', 'twopi'}
        darkBackgroundMode (bool, optional): change all color settings to make graph
            for dark background.
        width (int, optional): Graphviz graph width
        height (int, optional): Graphviz graph height

    Returns:
        TYPE: Description
    """
    if darkBackgroundMode:
        colorConfig = color_configs['dark']
    else:
        colorConfig = color_configs['light']

    global_styles, colorMapping = load_global_styles(colorConfig)

    g = gv.Digraph('G', format=imageFormat, engine=engine)

    # g.graph_attr['ratio']='fill'
    g.graph_attr['size'] = "{},{}".format(width, height)
    if imageFormat == 'png':
        g.graph_attr['dpi'] = '300'
    elif imageFormat == 'svg':
        g.graph_attr['dpi'] = '72'
    g.graph_attr['forcelabels'] = 'true'
    g.graph_attr['labelloc'] = 't'  # top or 'b'
    g.graph_attr['label'] = graphTitle
    g = apply_styles(g, global_styles)
    r_counter = 1

    # Automatically use scaling factor if flux > 10
    # (scaling factor is set to the nearest 100)
    all_f = [abs(f) for f in Pathway.fluxes]

    if max(all_f) > 10:
        scaleLineWidth = True
        scalingFactor = 10**(math.ceil(math.log(max(all_f), 10)))

    for rxn in Pathway.reactions:

        g.node(rxn.rid, shape='point',
               color=colorConfig['RXN_NODE_COLOR'],
               xlabel=rxn.rid + '; ' + str(abs(rxn.flux)),
               fontsize=REACTION_FONT_SIZE,
               fontcolor=colorConfig['REACTION_COLOR'])

        if scaleLineWidth:
            lineW = '%i' % (old_div(10 * abs(rxn.flux), scalingFactor) + 1)

        else:
            if rxn.flux >= 1 or rxn.flux <= -1:
                lineW = '%s' % (abs(rxn.flux) * 2)

            # penwidth for any flux between -1 < v < 1 is 1.
            else:
                lineW = '1'

        for met in rxn.reactants:
            if met in cofactorsList:
                cf = abs(rxn.metabolites[met])
                if cf > 1:
                    # show stoichiometric coefficient only when it is >= 2
                    clabel = '%i %s' % (cf, kegg_compound[met])
                else:
                    clabel = kegg_compound[met]
                g.node(met + '_' + str(r_counter),
                       shape=colorConfig['COFACTOR_SHAPE'],
                       color=colorMapping[met], style="filled", label=clabel)
                g.edge(met + '_' + str(r_counter), rxn.rid, penwidth=lineW,
                       weight='1', arrowhead="none",
                       color=colorConfig['EDGE_COLOR'])
            else:
                g.node(met, shape=colorConfig['NONCOFACTOR_SHAPE'],
                       label=kegg_compound[met], fontname='helvetica bold',
                       style="filled", color=colorConfig['NONCOFACTOR_COLOR'])
                g.edge(met, rxn.rid, penwidth=lineW, weight='2',
                       arrowhead="none", color=colorConfig['EDGE_COLOR'])

        for met in rxn.products:
            if met in cofactorsList:
                cf = abs(rxn.metabolites[met])
                if cf > 1:
                    clabel = '%i %s' % (cf, kegg_compound[met])
                else:
                    clabel = kegg_compound[met]
                g.node(met + '_' + str(r_counter),
                       shape=colorConfig['COFACTOR_SHAPE'],
                       color=colorMapping[met], style="filled", label=clabel)
                g.edge(rxn.rid, met + '_' + str(r_counter), weight='1',
                       penwidth=lineW, color=colorConfig['EDGE_COLOR'])
            else:
                g.node(met, shape=colorConfig['NONCOFACTOR_SHAPE'],
                       label=kegg_compound[met], fontname='helvetica bold',
                       style="filled", color=colorConfig['NONCOFACTOR_COLOR'])
                g.edge(rxn.rid, met, penwidth=lineW, weight='2',
                       color=colorConfig['EDGE_COLOR'])
        r_counter += 1

    if imageFileName is None:
        imageFileName = Pathway.name
    g.render(imageFileName, cleanup=cleanup)

    return g


def test_drawpathway():
    """
    Test for drawing pathways. Create two version of the pathways.
    """
    testpathway = {'flux': [-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0,
                            1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 2.0, 1.0,
                            1.0, 1.0, -1.0, 1.0],
                   'iteration': 1,
                   'reaction_id': ['R00200', 'R00300', 'R00658', 'R01059',
                                   'R01063', 'R01512', 'R01518', 'R01519',
                                   'R01538', 'R08570', 'EX_glc', 'EX_nad',
                                   'EX_adp', 'EX_phosphate', 'EX_pyruvate',
                                   'EX_nadh', 'EX_atp', 'EX_h2o', 'EX_nadp',
                                   'EX_nadph']}

    p1 = Pathway(id=1, name='OptStoic',
                 reaction_ids=testpathway['reaction_id'],
                 fluxes=testpathway['flux'])
    # outputFilename = 'OptStoic'
    # current_dir = dirname(abspath(__file__))
    # outputFilepath = normpath(join(current_dir,'../result', outputFilename))

    logging.info("Creating 'res' folder in current directory if not exist...")
    outputFilepath = 'res'
    try:
        os.makedirs(outputFilepath)
    except OSError:
        if not os.path.isdir(outputFilepath):
            raise

    draw_pathway(p1, os.path.join(outputFilepath, 'test_drawpath_light'),
                 imageFormat='png', graphTitle='test_pathway_light',
                 scaleLineWidth=False, scalingFactor=200.0,
                 darkBackgroundMode=False)

    draw_pathway(p1, os.path.join(outputFilepath, 'test_drawpath_dark'),
                 imageFormat='png', graphTitle='test_pathway_dark',
                 scaleLineWidth=False, scalingFactor=200.0,
                 darkBackgroundMode=True)

    logging.info("Testing drawpathway.py: Pass\n")
    return None
# ------------------------------------------------------


if __name__ == "__main__":

    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    test()
