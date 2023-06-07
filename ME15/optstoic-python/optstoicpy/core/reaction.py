from builtins import range
from builtins import object
from .config import rxnSji
from optstoicpy.script.utils import create_logger


class Reaction(object):
    """Reaction class

    Attributes:
        equation (TYPE): Description
        flux (TYPE): Description
        metabolites (TYPE): Description
        reversible (TYPE): Description
        rid (TYPE): Description
    """

    def __init__(self,
                 rid=None,
                 flux=1,
                 metabolites={},
                 equation='',
                 reversible=True,
                 logger=None):

        if logger is None:
            self.logger = create_logger('core.Reaction')
        else:
            self.logger = logger

        self.rid = rid
        self.flux = flux
        self.metabolites = metabolites
        self.equation = equation
        self.reversible = reversible

    def autoset_metabolites(self):
        if len(self.metabolites) > 0:
            self.logger.warning("Metabolites exists!")
        else:
            self.logger.info("Retrieving metabolites from default database")
            self.metabolites = rxnSji[self.rid]
        return self.metabolites

    @property
    def reactants(self):
        if self.flux > 0:
            return [k for k, v in list(self.metabolites.items()) if v < 0]
        else:
            return [k for k, v in list(self.metabolites.items()) if v > 0]

    @property
    def products(self):
        if self.flux > 0:
            return [k for k, v in list(self.metabolites.items()) if v > 0]
        else:
            return [k for k, v in list(self.metabolites.items()) if v < 0]

    def set_equation(self):
        """Write equation in the direction of the flux.
        The main purpose is to simplify downstream MDF/protein cost analysis.

        Returns:
            TYPE: Description
        """
        if len(self.equation) != 0:
            self.logger.warning("Equation exists!")
        else:
            if len(self.metabolites) == 0:
                self.logger.info(
                    "Metabolites are not available! Auto-updating metabolites...")
                self.autoset_metabolites()

            temp_list = []
            for cpd in sorted(self.reactants):
                coeff = abs(self.metabolites[cpd])
                if coeff == 1:
                    temp_list.append(cpd)
                else:
                    temp_list.append('%1.0f %s' % (coeff, cpd))

            eqStr = ' + '.join(temp_list)
            eqStr += ' <=> '

            temp_list = []
            for cpd in sorted(self.products):
                coeff = abs(self.metabolites[cpd])
                if coeff == 1:
                    temp_list.append(cpd)
                else:
                    temp_list.append('%1.0f %s' % (coeff, cpd))
            eqStr += ' + '.join(temp_list)
            self.equation = eqStr
        return self.equation

    @classmethod
    def create_Reaction_list_from_dict(cls, dataDict, excludeExchangeRxn=True):
        """
        Make a list of Reaction object from dataDict, excluding exchange reactions.

        E.g.
        dataDict = {'reaction_id': ['R00001', 'R00002'], 'flux': [-1, 1]}
        output = [Reaction('R00001'), Reaction('R00002')]

        Args:
            dataDict (TYPE): ictionary with reaction_id and flux
            excludeExchangeRxn (bool, optional): Exclude all exchange reactions in the list. Default to True.

        Returns:
            TYPE: Description
        """
        RxnObjList = []
        for i in range(len(dataDict['reaction_id'])):
            if excludeExchangeRxn:
                if 'EX_' in dataDict['reaction_id'][i]:
                    continue
            tempRxn = cls(dataDict['reaction_id'][i], dataDict['flux'][i])
            # Get the metabolites dictionary {'C00001': -1, ...} for each
            # reaction
            tempRxn.metabolites = rxnSji[tempRxn.rid]
            RxnObjList.append(tempRxn)
        return RxnObjList

    def __str__(self):
        return "Reaction('%s')" % self.rid

    def __repr__(self):
        return "Reaction('%s')" % self.rid
