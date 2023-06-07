from __future__ import division
from builtins import map
from builtins import zip
from builtins import object
from .reaction import Reaction
from .config import cofactors, default_params, rxnSji
import os
from collections import OrderedDict
from optstoicpy.script.utils import create_logger


class Pathway(object):
    """OptStoic Pathway class"""

    def __init__(self,
                 id=None,
                 name=None,
                 reaction_ids=[],
                 fluxes=[],
                 reactions=None,
                 sourceSubstrateID='C00031',
                 endSubstrateID='C00022',
                 total_flux_no_exchange=None,
                 note={},
                 logger=None):
        """
        A Pathway instance can be initialized by either
            (a) List of reaction_ids and fluxes (Let reactions = None)
            (b) List of reaction instances (reaction_ids and fluxes
                will be populated)
        Args:
            id (None, optional): Pathway id
            name (None, optional): Pathway name
            reaction_ids (list, optional): A list of reaction IDs (kegg_id) in the pathway
            fluxes (list, optional): A list of reaction fluxes corresponding to the reaction_ids
            reactions (None, optional): Description
            sourceSubstrateID (str, optional): Kegg compound ID of the source metabolite of
                the pathway
            endSubstrateID (str, optional): Kegg compound ID of the end metabolite of the pathway
            total_flux_no_exchange (None, optional): Sum of absolute flux through the pathway (Exclude export reactions)
            note (dict, optional): (For debugging purpose) modelstat and solvestat can be added
        """
        if logger is None:
            self.logger = create_logger('core.Pathway')
        else:
            self.logger = logger
        self.id = id
        self.name = name
        self.note = note

        # Iniatilize pathway object using list of reaction_ids and fluxes
        if reactions is None:
            # Check if both reaction_ids and fluxes are list and
            # contain the same number of item
            assert (isinstance(reaction_ids, list) == 1)
            assert (isinstance(fluxes, list) == 1)
            assert len(reaction_ids) == len(
                fluxes), "number of reactions must equal number of fluxes!"

            # Change EX_h+ to EX_hplus as optstoic pulp fail to read "+" symbol
            self.reaction_ids = ["EX_hplus" if x ==
                                 "EX_h+" else x for x in reaction_ids]

            self.fluxes = fluxes
            # Create list of reaction objects
            self.reactions = Reaction.create_Reaction_list_from_dict(
                {'reaction_id': self.reaction_ids, 'flux': self.fluxes},
                excludeExchangeRxn=True)

        # Iniatilize pathway object using list of reaction objects
        else:
            self.reactions = reactions
            self.fluxes = [r.flux for r in reactions]
            self.reaction_ids = [r.rid for r in reactions]

        self.reaction_ids_no_exchange = [
            r for r in reaction_ids if 'EX_' not in r]

        if not total_flux_no_exchange:
            self.total_flux_no_exchange = sum(map(
                abs, [r.flux for r in self.reactions]))
        else:
            self.total_flux_no_exchange = total_flux_no_exchange

        self.rxn_flux_dict = dict(list(zip(self.reaction_ids, self.fluxes)))

        try:
            self.nATP = self.rxn_flux_dict['EX_atp']
        except BaseException:
            self.nATP = None
        self.sourceSubstrateID = sourceSubstrateID
        self.endSubstrateID = endSubstrateID

    def get_pathway_dict(self):
        """
        return a dictionary of the {reaction:flux}
        """
        return dict(list(zip(self.reaction_ids, self.fluxes)))

    def update_nATP(self):

        try:
            self.nATP = self.rxn_flux_dict['EX_atp']
        except BaseException:
            self.nATP = None

    def get_total_flux(self):
        """
        return total flux through the pathway
        """
        return sum(map(abs, self.fluxes))

    def get_total_flux_no_exchange(self):
        """
        return total flux through the pathway excluding exchange reactions
        """
        return self.total_flux_no_exchange

    def get_modelstat(self):
        if 'modelstat' in self.note:
            return self.note['modelstat']
        else:
            return None

    def get_solvestat(self):
        if 'solvestat' in self.note:
            return self.note['solvestat']
        else:
            return None

    def get_time(self):
        if 'time' in self.note:
            return self.note['time']
        else:
            return None

    def get_reaction_involving_reactant(self, substrate_ID):
        """get reactions that involve reactant substrate_ID"""
        output = []
        for rxn in self.reactions:
            if substrate_ID in rxn.reactants:
                output.append(rxn)
        return output

    def rearrange_reaction_order(self):
        """
        Try to implement a topological sorting of the pathway.
        This was done in a very early stage of the project.
        (Todo: use a different algorithm, e.g., graph-based breadth-first search algorithm)

        """
        # Exclude exchange reaction from being rearranged
        sortedRxn = []
        flag = 1

        next_substrate = [self.sourceSubstrateID]
        # Find substrateID in list of reactants
        while flag == 1:
            rstore = []
            for subs in next_substrate:
                if subs == self.endSubstrateID:
                    flag = 0
                    break
                r1 = self.get_reaction_involving_reactant(subs)

                if len(r1) == 0:
                    # print "\nWarning: No reaction found using %s" %subs
                    if len(next_substrate) == 1:
                        flag = 0
                        break
                else:
                    sortedRxn.extend([r.rid for r in r1])
                rstore += r1
            next_substrate = []
            for rxn in rstore:
                next_substrate.extend(rxn.products)
            next_substrate = list(set(next_substrate) - cofactors)

        # Add all reactions to the sorted reactions
        # (as the duplicates will be removed in the following command)
        sortedRxn += self.reaction_ids
        # Make unique ordered list
        sortedRxnUnique = list(OrderedDict.fromkeys(sortedRxn))

        # Raise error if the number of reactions changes after processing
        assert len(sortedRxnUnique) == len(self.reaction_ids), "Error: \
        the number of reactions does not match after processing"

        # Sort all the flux according to the order of the reaction ID
        sortedRxnFlux = [self.rxn_flux_dict[rxn] for rxn in sortedRxnUnique]

        self.reaction_ids = sortedRxnUnique
        self.fluxes = sortedRxnFlux
        self.reactions = Reaction.create_Reaction_list_from_dict(
            {'reaction_id': sortedRxnUnique, 'flux': sortedRxnFlux})
        return self

    # @staticmethod
    def get_pathway_similarity_index(self, pathway2):
        """Calculate the jaccard index of two pathways"""
        a = set(self.reaction_ids)
        b = set(pathway2.reaction_ids)
        idscore = len(a & b) / len(a | b)
        return idscore

    def get_pathway_similarity_index_no_exchange(self, pathway2):
        """Calculate the jaccard index of two pathways"""
        a = set(self.reaction_ids_no_exchange)
        b = set(pathway2.reaction_ids_no_exchange)
        idscore = len(a & b) / len(a | b)
        return idscore

    def is_same_pathway_with(self, another_pathway):
        idscore = self.get_pathway_similarity_index(another_pathway)
        if idscore == 1:
            return 1
        else:
            return 0

    def to_dict(self):
        return dict(pathway=self.get_pathway_dict(),
                    num_reaction=len(self.reaction_ids),
                    total_flux_no_exchange=self.get_total_flux_no_exchange(),
                    modelstat=self.get_modelstat(),
                    solvestat=self.get_solvestat(),
                    time=self.get_time())

    @staticmethod
    def pathways_to_dict(list_of_pathways):
        """
        Output list of Pathway instances as dictionary
        """
        output = {}
        for p in list_of_pathways:
            output[p.id] = p.to_dict()

        return output

    def __repr__(self):
        return "<OptStoicPathway(id='%s', numRxn='%s', nATP='%s')>" % (
            self.id, len(self.reaction_ids), self.nATP)

# ----------------------------------------------------------------------------


def generate_kegg_model(pathway,
                        params=default_params,
                        filehandle=None,
                        add_ratio_constraints=False):
    """
    Convert the pathway to KEGG model format
    (as input for Component Contribution/MDF)

    Args:
        pathway (TYPE): A pathway instance
        params (TYPE, optional): KEGG model parameters (default parameters are given)
        filehandle (None, optional): If a text file handle is provided,
            it writes the model text to the file (default None)
        add_ratio_constraints (bool, optional): Description

    Returns:
        TYPE: Description

    """
    params['ENTRY'] = "{0}_{1}ATP_P{2}".format(pathway.name,
                                               pathway.nATP, pathway.id)
    params['NAME'] = "{0}_{1}ATP_P{2}".format(pathway.name,
                                              pathway.nATP, pathway.id)

    modeltext = """\
ENTRY\t\t{ENTRY}
SKIP\t\t{SKIP}
NAME\t\t{NAME}
PH\t\t\t{PH}
I\t\t\t{I}
T\t\t\t{T}
C_RANGE\t\t{C_RANGE[0]:.0e} {C_RANGE[1]:.0e}\n""".format(**params)

    all_bounds = params['BOUND']

    if add_ratio_constraints:
        all_bounds = params['RATIO_BOUND']

        # write the ratios
        for i, (cids, ratios) in enumerate(sorted(params['RATIO'].items())):
            if i == 0:
                modeltext += "RATIO\t\t{C[0]} {C[1]} {B[0]:e} {B[1]:e}\n".format(
                    C=cids, B=ratios)
            else:
                modeltext += "\t\t\t{C[0]} {C[1]} {B[0]:e} {B[1]:e}\n".format(
                    C=cids, B=ratios)

    # write the bounds
    for i, (cid, bounds) in enumerate(sorted(all_bounds.items())):
        if i == 0:
            modeltext += "BOUND\t\t{0} {1[0]:e} {1[1]:e}\n".format(cid, bounds)
        else:
            modeltext += "\t\t\t{0} {1[0]:e} {1[1]:e}\n".format(cid, bounds)

    # write the reactions (in the direction of the flux)
    for i, rxn in enumerate(pathway.reactions):
        rxn.set_equation()
        if i == 0:
            modeltext += "REACTION\t{0} {1} (x{2:1.2f})\n".format(
                rxn.rid, rxn.equation, abs(rxn.flux))
        else:
            modeltext += "\t\t\t{0} {1} (x{2:1.2f})\n".format(
                rxn.rid, rxn.equation, abs(rxn.flux))

    modeltext += "///\n"
    if filehandle:
        filehandle.write(modeltext)

    return modeltext
