from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import object
import os
import json
import copy
import pandas as pd
from optstoicpy.script import gams_parser
from optstoicpy.script.utils import create_logger

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(
    os.path.join(CURRENT_DIR, '../data/', 'optstoic_db_v3')
)


class BaseReactionDatabase(object):
    """The initial reaction database to be used for pre-processing.

    This converts GAMS model files to Python while retaining the
    model structure that GAMS users are familiar with.
    """
    REACTION_TYPE = {0: 'Forward irreversible', 1: 'Reversible',
                     2: 'Reverse irreverisble', 4: 'Export reaction'}

    def __init__(
            self,
            data_filepath='../data/',
            dbdict_json=None,
            dbdict_gams=None,
            logger=None):

        if logger is None:
            self.logger = create_logger('core.BaseDatabase')
        else:
            self.logger = logger

        self.data_filepath = data_filepath
        self.dbdict_json = dbdict_json
        self.dbdict_gams = dbdict_gams

        # initalize
        self.reactions = []
        self.metabolites = []
        self.S = {}
        self.Sji = {}
        self.rxntype = []
        self.user_defined_export_rxns = []
        self._S_df = None

    def load(self):
        # Method 1: JSON approach
        if self.dbdict_json is not None:
            self.logger.debug('Reading S matrix from JSON...')
            self.Sji = json.load(
                open(
                    os.path.join(
                        self.data_filepath,
                        self.dbdict_json['Sji']),
                    'r+'))
            self.S = self.transpose_S(self.Sji)
            self.reactions = sorted(self.Sji.keys())
            self.internal_rxns = copy.deepcopy(self.reactions)
            self.metabolites = sorted(self.S.keys())
            self.logger.debug('Reading reaction type file...')
            self.rxntype = json.load(
                open(
                    os.path.join(
                        self.data_filepath,
                        self.dbdict_json['reactiontype']),
                    'r+'))

        # Method 2: Standard GAMS input file
        else:
            self.logger.debug('Reading S matrix from txt...')
            self.S = gams_parser.convert_parameter_table_to_dict(
                os.path.join(self.data_filepath,
                             self.dbdict_gams['Sji'])
            )
            self.Sji = self.transpose_S(self.S)

            # Load reactions
            self.logger.debug('Reading reaction file...')
            self.reactions = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict_gams['reaction'])
            )
            # Create internal reactions list
            self.internal_rxns = copy.deepcopy(self.reactions)

            self.logger.debug('Reading metabolite file...')
            self.metabolites = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict_gams['metabolite'])
            )

            self.logger.debug('Reading reaction type file...')
            self.rxntype = gams_parser.convert_parameter_list_to_dict(
                os.path.join(self.data_filepath, self.dbdict_gams['reactiontype']),
                datadict=None
            )
        self.validate()

    def validate(self):
        """Validate the S matrix, the reaction vector and the reaction type vector.

        Raises:
            Exception: Description
        """
        self.logger.info("Validating database")

        if set(self.Sji.keys()) != set(self.reactions):
            raise Exception("The number of reactions do not match!")

        if set(self.S.keys()) != set(self.metabolites):
            raise Exception("The number of metabolites do not match!")

        for rxn in self.reactions:
            assert self.rxntype.get(
                rxn, None) is not None, "%s does not have rxntype assigned!" % rxn

        if None in set(self.rxntype.values()):
            raise Exception("Some reaction type is not assigned!")

    @staticmethod
    def transpose_S(Sji):
        """Tranpose Sji into Sij and also Sij to Sji dictionary."""
        # Update to pandas 0.19 (using sparse dataframe)
        df_Sji = pd.DataFrame(Sji).T
        Sij = dict(
            (k, v.dropna().to_dict()) for k, v in df_Sji.items()
        )
        return Sij

    def create_or_update_S_df(self):
        """Create Pandas.DataFrame of the Sij matrix
        """
        self._S_df = pd.DataFrame(self.Sji).fillna(0)

    def refresh_database(self, previous_operations_on='Sji'):
        """Afer loading the database, if any operations were
        performed on S(ij) or Sji, refresh the database to reflect
        the changes on all related attribute.
        TODO: This can be simplified

        Args:
            previous_operations_on (str, optional): Description
        """
        if previous_operations_on == 'Sji':
            # re-create Sij from Sji
            self.S = self.transpose_S(self.Sji)
        elif previous_operations_on == "S":
            # re-create Sji from S
            self.Sji = self.transpose_S(self.S)
        else:
            raise Exception(
                "The previous_operations_on argument must be from the list: ['S', 'Sji']")

        # update S df
        self.create_or_update_S_df()
        # update metabolites
        self.metabolites = sorted(self.S.keys())
        # update reactions
        self.reactions = sorted(self.Sji.keys())
        self.validate()

    @property
    def S_df(self):
        """Return a Pandas.DataFrame of the S matrix
        metabolites = self.S_df.index.tolist()
        reactions = self.db.S_df.columns.tolist()
        Smat = self.S_df.as_matrix()
        Returns:
            `Pandas.DataFrame`: Description
        """
        if self._S_df is None:
            self.create_or_update_S_df()
        return self._S_df

    @staticmethod
    def to_json(Sdict, filepath):
        with open(filepath, 'w+') as fp:
            json.dump(Sdict, fp, sort_keys=True, indent=4)

    def to_mat_file():
        """write to matlab file"""
        pass

    def get_reaction_type(self, rid, verbose=True):
        try:
            self.rxntype[rid]
        except BaseException:
            self.logger.warning("Reaction %s not in database!" % rid)
            return None
        else:
            if verbose:
                print(
                    "Reaction: {0} is ({1}) {2}".format(
                        rid,
                        self.rxntype[rid],
                        self.REACTION_TYPE.get(
                            self.rxntype[rid])))
            return self.rxntype[rid]

    def set_reaction_type(self, rid, rxntype):
        try:
            t0 = self.rxntype[rid]
            self.rxntype[rid] = rxntype
        except KeyError:
            self.logger.error('Reaction %s not in database!' % rid)
        else:
            self.logger.debug(
                'Reaction %s has been updated from %s to %s.' %
                (rid, self.REACTION_TYPE.get(t0), self.REACTION_TYPE.get(rxntype)))

    def extend_S_from_gams_inputfile(self, filename=None):
        """Extend S matrix using the gams S matrix format.
        Args:
            filename (None, optional): The name of the inputfile.
        """
        self.S = gams_parser.convert_parameter_table_to_dict(
            os.path.join(self.data_filepath, filename),
            Sdict=self.S)

        self.refresh_database(previous_operations_on='S')

    def update_S(self, extension_dict, default_reactiontype=None):
        """Add new reactions using dictionary
        {'EX_glc': {'C00031': -1.0},
         'EX_nad': {'C00003': -1.0}}
        Args:
            extension_dict (TYPE): Description

        Returns:
            TYPE: Description
        """
        temp_rxn = []
        for met, entries in extension_dict.items():
            if met not in self.S:
                self.S[met] = {}
                self.metabolites.append(met)
            for rxn, coeff in entries.items():
                self.S[met][rxn] = float(coeff)
                if rxn not in self.reactions:
                    self.reactions.append(rxn)
                    temp_rxn.append(rxn)
                    self.rxntype[rxn] = default_reactiontype

        self.refresh_database(previous_operations_on='S')
        return self.S, temp_rxn

    def set_database_export_reaction(self, export_reactions_Sij_dict):

        _, temp_rxn = self.update_S(
            export_reactions_Sij_dict, default_reactiontype=4)
        if len(self.user_defined_export_rxns) != 0:
            self.logger.warning("Warning: The current list of export reactions\
                will be replaced! %s" % str(self.user_defined_export_rxns))
        self.user_defined_export_rxns = temp_rxn
        # for rxn in self.user_defined_export_rxns:
        #     self.set_reaction_type(rxn, 4)

    def update_rxntype(self, new_reaction_type_dict):
        for (r, rtype) in new_reaction_type_dict.items():
            self.set_reaction_type(r, rtype)
        return self.rxntype

    def remove_reaction(self, reaction_id, refresh_database=True):
        """Remove a reaction by reaction_id from the database.

        Args:
            reaction_id (TYPE): The reaction_id
            refresh_database (bool, optional):  If True, the S matrix will be updated,
                but this is slow when removing a large number of reactions. In that case,
                perform the refresh_database after the last iteration.
        """
        assert reaction_id in self.Sji, "The reaction_id provided do not present in the database!"
        # remove from S matrix
        self.Sji.pop(reaction_id, None)
        # remove reactions
        self.reactions.remove(reaction_id)
        # remove from internal reactions, hard coded for the moment
        if not reaction_id.startswith("EX_"):
            self.internal_rxns.remove(reaction_id)
        # remove from reaction type
        self.rxntype.pop(reaction_id, None)

        if refresh_database:
            self.refresh_database(previous_operations_on='Sji')

        self.logger.debug(
            "Reaction %s removed from the database." %
            reaction_id)

    def __repr__(self):
        return "BaseReactionDatabase"


class Database(BaseReactionDatabase):
    """The optstoic Database class: loading database from GAMS input files.
    TODO: Use cobrapy Model/interconvert between different modes of input.
    """

    def __init__(
            self,
            description='',
            data_filepath='../data/',
            dbdict_json=None,
            dbdict_gams=None,
            blocked_rxns=None,
            excluded_reactions=None,
            reduce_model_size=True,
            logger=None):
        """Summary

        Args:
            description (str, optional): Description of the Database
            data_filepath (str, optional): Description
            dbdict_json (None, optional): filename for json
            dbdict_gams (None, optional): filename for gams input
            blocked_rxns (None, optional): A list of blocked reactions
            excluded_reactions (None, optional): A list of reactions to be excluded
                from optstoic result
            reduce_model_size (bool, optional): If True, remove blocked_rxns from
                the S matrix to speed up optstoic analysis
            logger (None, optional): Description
        """
        if logger is None:
            logger = create_logger('core.Database')

        self.description = description

        super(Database, self).__init__(
            data_filepath=data_filepath,
            dbdict_json=dbdict_json,
            dbdict_gams=dbdict_gams,
            logger=logger)

        # Initalize child-only attributes
        self.loops = []
        self.Ninternal = {}
        self.all_excluded_reactions = []
        self.excluded_reactions = excluded_reactions

        if blocked_rxns is not None:
            assert isinstance(
                blocked_rxns, list), "blocked_rxns must be a list"
        self.blocked_rxns = blocked_rxns
        self.reduce_model_size = reduce_model_size

    def load(self):

        # Method 1: JSON approach
        if self.dbdict_json is not None:
            self.logger.debug('Reading S matrix from JSON...')
            self.Sji = json.load(
                open(
                    os.path.join(
                        self.data_filepath,
                        self.dbdict_json['Sji']),
                    'r+'))
            self.S = self.transpose_S(self.Sji)
            self.reactions = sorted(self.Sji.keys())
            self.internal_rxns = copy.deepcopy(self.reactions)
            self.metabolites = sorted(self.S.keys())
            self.logger.debug('Reading reaction type file...')
            self.rxntype = json.load(
                open(
                    os.path.join(
                        self.data_filepath,
                        self.dbdict_json['reactiontype']),
                    'r+'))

            self.logger.debug('Reading Nint(loop, j) from JSON...')
            self.Ninternal = json.load(
                open(
                    os.path.join(
                        self.data_filepath,
                        self.dbdict_json['Nint']),
                    'r+'))

            self.loops = sorted(self.Ninternal.keys())
        # Method 2: Standard GAMS input file
        else:
            self.logger.debug('Reading S matrix from txt...')
            self.S = gams_parser.convert_parameter_table_to_dict(
                os.path.join(self.data_filepath,
                             self.dbdict_gams['Sji'])
            )
            self.Sji = self.transpose_S(self.S)

            self.logger.debug('Reading metabolite file...')
            self.metabolites = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict_gams['metabolite'])
            )
            # Load reactions
            self.logger.debug('Reading reaction file...')
            self.reactions = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict_gams['reaction'])
            )
            self.internal_rxns = copy.deepcopy(self.reactions)

            self.logger.debug('Reading reaction type file...')
            self.rxntype = gams_parser.convert_parameter_list_to_dict(
                os.path.join(self.data_filepath, self.dbdict_gams['reactiontype']),
                datadict=None
            )
            self.logger.debug('Reading Nint(loop, j) from txt...')
            self.Ninternal = gams_parser.convert_parameter_table_to_dict(
                os.path.join(self.data_filepath, self.dbdict_gams['Nint']))

            self.logger.debug('Reading loop file...')
            self.loops = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict_gams['loops']))

        if self.excluded_reactions is not None:
            self.all_excluded_reactions = list(
                set(self.excluded_reactions + self.blocked_rxns)
            )
        else:
            self.all_excluded_reactions = self.blocked_rxns

        self.validate()

        if self.reduce_model_size:
            self.remove_blocked_reactions()
            self.validate()

    def remove_blocked_reactions(self):
        self.logger.warning("Removing blocked reactions to reduce model size!")

        loop_rxns = [list(v.keys()) for v in list(self.Ninternal.values())]
        loop_rxns = set([rid for sublist in loop_rxns for rid in sublist])
        assert len(loop_rxns & set(self.blocked_rxns)
                   ) == 0, "Blocked reactions must not present in loops"

        for rxn in self.blocked_rxns:
            self.remove_reaction(rxn, refresh_database=False)

        self.refresh_database(previous_operations_on='Sji')

    def __repr__(self):
        return "OptStoic Database(Description='%s')" % self.description


def load_custom_reactions_to_be_excluded():
    """A list of undesirable reactions that are specific
    to the glycolysis study

    Returns:
        TYPE: Description
    """
    NTP_involving_rxns = gams_parser.convert_set_to_list(
        os.path.join(DATA_DIR, 'NTP_and_AMP_reactions.txt')
    )
    cofactor_only_rxns = gams_parser.convert_set_to_list(
        os.path.join(DATA_DIR, 'cofactor_only_reactions.txt')
    )
    cofactor_only_rxns.append('R10092')

    methylglyoxal_rxns = [
        'R00203',
        'R00205',
        'R01016',
        'R02260',
        'R02527',
        'R02528',
        'R02529',
        'R02530',
        'R02531',
        'R07183',
        'R09796',
        'R10049',
        'R10050'
    ]

    other_undesirable_rxns = [
        # Bicarbonate and pyrrole cycle
        'R09794',
        'R09795',
        # Undesirable glucose uptake loop
        'R00305',
        'R00874',
        'R07359',
        'R00837',
        'R09749',
        'R03075',
        'R02985',
        'R02558',
        'R01555',
        'R02727',
        'R00010',
        'R02778',
        'R08946',
        'R00306'
    ]

    excluded_reactions = (NTP_involving_rxns +
                          methylglyoxal_rxns +
                          cofactor_only_rxns +
                          other_undesirable_rxns)

    all_excluded_reactions = list(set(excluded_reactions))
    return all_excluded_reactions


def load_base_reaction_db(
        user_defined_export_rxns_Sji=None,
        logger=None):
    """Load the base reaction database with all reactions
    (i.e., no blocked reactions and no loops)

    Args:
        user_defined_export_rxns_Sji (None, optional): Description

    Returns:
        :obj:`BaseReactionDatabase`:
    """
    if logger is None:
        logger = create_logger(
            name="optstoicpy.core.database.load_base_reaction_db")

    # get reactions that are manually curated to be excluded for glycolysis
    # study
    excluded_reactions = load_custom_reactions_to_be_excluded()

    dbdict_json = {
        'Sji': 'optstoic_v3_Sji_dict.json',
        'reactiontype': 'optstoic_v3_reactiontype.json'
    }
    dbdict_gams = {
        'Sji': 'optstoic_v3_Sij.txt',
        'reaction': 'optstoic_v3_reactions.txt',
        'metabolite': 'optstoic_v3_metabolites.txt',
        'reactiontype': 'optstoic_v3_reactiontype.txt'
    }

    DB = BaseReactionDatabase(
        data_filepath=DATA_DIR,
        dbdict_json=dbdict_json,
        dbdict_gams=dbdict_gams)

    DB.load()

    # Update reaction type
    # Update reaction type  = 0
    irreversible_fwd_rxns = gams_parser.convert_set_to_list(os.path.join(
        DATA_DIR, 'optstoic_v3_ATP_irreversible_forward_rxns.txt')
    )

    new_reaction_type_dict = dict(list(zip(
        irreversible_fwd_rxns, [0] * len(irreversible_fwd_rxns)))
    )
    # Update reaction type  =  2
    irreversible_bwd_rxns = gams_parser.convert_set_to_list(os.path.join(
        DATA_DIR, 'optstoic_v3_ATP_irreversible_backward_rxns.txt')
    )

    new_reaction_type_dict.update(dict(
        list(zip(irreversible_bwd_rxns, [2] * len(irreversible_bwd_rxns))))
    )

    DB.update_rxntype(new_reaction_type_dict)

    # Add a list of export reactions and the metabolites
    if user_defined_export_rxns_Sji is not None:
        user_defined_export_rxns_Sij = Database.transpose_S(
            user_defined_export_rxns_Sji
        )

        DB.set_database_export_reaction(user_defined_export_rxns_Sij)
    return DB


def load_db_v3(
    reduce_model_size=True,
    user_defined_export_rxns_Sji={
        'EX_glc': {'C00031': -1.0},
        'EX_nad': {'C00003': -1.0},
        'EX_adp': {'C00008': -1.0},
        'EX_phosphate': {'C00009': -1.0},
        'EX_pyruvate': {'C00022': -1.0},
        'EX_nadh': {'C00004': -1.0},
        'EX_atp': {'C00002': -1.0},
        'EX_h2o': {'C00001': -1.0},
        'EX_hplus': {'C00080': -1.0},
        'EX_nadp': {'C00006': -1.0},
        'EX_nadph': {'C00005': -1.0}
    },
    logger=None
):
    """Load OptStoic database v3

    Returns:
        TYPE: Description

    Args:
        reduce_model_size (bool, optional): True if you want to reduce the size of the
            model by removing blocked reactions from the S matrix.
        user_defined_export_rxns_Sji (dict, optional): The list of export reactions that
            need to be added to the model for metabolite exchange (i.e., any metabolite
            that participate in the design equation)
    """
    if logger is None:
        logger = create_logger(name="optstoicpy.core.database.load_db_v3")

    # get reactions that are manually curated to be excluded for glycolysis
    # study
    excluded_reactions = load_custom_reactions_to_be_excluded()

    dbdict_json = {
        'Sji': 'optstoic_v3_Sji_dict.json',
        'Nint': 'optstoic_v3_Nint.json',
        'reactiontype': 'optstoic_v3_reactiontype.json'
    }

    dbdict_gams = {
        'Sji': 'optstoic_v3_Sij.txt',
        'reaction': 'optstoic_v3_reactions.txt',
        'metabolite': 'optstoic_v3_metabolites.txt',
        'reactiontype': 'optstoic_v3_reactiontype.txt',
        'loops': 'optstoic_v3_loops_nocofactor.txt',
        'Nint': 'optstoic_v3_null_sij_nocofactor.txt',
        'blocked_rxns': 'optstoic_v3_blocked_reactions_0to5ATP.txt',
    }

    logger.debug('Reading blocked reactions file...')
    if 'blocked_rxns' in dbdict_gams:
        blocked_rxns = gams_parser.convert_set_to_list(
            os.path.join(DATA_DIR, dbdict_gams['blocked_rxns'])
        )
    else:
        blocked_rxns = None

    DB = Database(
        description='v3',
        data_filepath=DATA_DIR,
        dbdict_json=dbdict_json,
        dbdict_gams=dbdict_gams,
        blocked_rxns=blocked_rxns,
        excluded_reactions=excluded_reactions,
        reduce_model_size=reduce_model_size)

    DB.load()

    # Update reaction type
    # Update reaction type  = 0
    irreversible_fwd_rxns = gams_parser.convert_set_to_list(os.path.join(
        DATA_DIR, 'optstoic_v3_ATP_irreversible_forward_rxns.txt')
    )

    new_reaction_type_dict = dict(list(zip(
        irreversible_fwd_rxns, [0] * len(irreversible_fwd_rxns)))
    )
    # Update reaction type  =  2
    irreversible_bwd_rxns = gams_parser.convert_set_to_list(os.path.join(
        DATA_DIR, 'optstoic_v3_ATP_irreversible_backward_rxns.txt')
    )

    new_reaction_type_dict.update(dict(
        list(zip(irreversible_bwd_rxns, [2] * len(irreversible_bwd_rxns))))
    )

    DB.update_rxntype(new_reaction_type_dict)

    # user_defined_export_rxns = ['EX_glc', 'EX_nad', 'EX_adp',
    #                             'EX_phosphate', 'EX_pyruvate', 'EX_nadh',
    #                             'EX_atp', 'EX_h2o', 'EX_hplus', 'EX_nadp',
    #                             'EX_nadph']

    # Add a list of export reactions and the metabolites
    if user_defined_export_rxns_Sji is not None:
        user_defined_export_rxns_Sij = Database.transpose_S(
            user_defined_export_rxns_Sji
        )

        DB.set_database_export_reaction(user_defined_export_rxns_Sij)

    return DB


if __name__ == "__main__":
    DB = load_db_v3()
