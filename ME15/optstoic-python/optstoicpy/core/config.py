import os
from os.path import dirname, abspath, normpath
import pandas as pd
import json
# import pdb

current_dir = dirname(abspath(__file__))
data_dir = normpath(os.path.join(current_dir, '../data'))

# Load cofactors list
cofactors_df = pd.read_csv(os.path.join(data_dir, 'cofactors.csv'))
cofactorsList = cofactors_df['KEGG_ID'].tolist()
cofactors = set(cofactorsList)

# reaction Sij
# rxnSij = json.load(open(
#    os.path.join(data_dir, 'optstoic_db_v2',
#                 '20160616_optstoic_Sji_dict.json')
#    ,'r+'))
rxnSji = json.load(open(os.path.join(data_dir,
                                     'optstoic_db_v3',
                                     'optstoic_v3_Sji_dict.json'), 'r+'))

kegg_compound = json.load(open(os.path.join(data_dir,
                                            'kegg_compound.json'), 'r+'))

# Kegg_model default argument
default_bound = {
    'C00001': [1, 1],
    'C00002': [5e-3, 5e-3],
    'C00008': [5e-4, 5e-4],
    'C00009': [1e-2, 1e-2],
    'C00020': [5e-4, 5e-4],
    'C00003': [5e-3, 5e-3],
    'C00004': [5e-5, 5e-5],
    'C00005': [5e-4, 5e-4],
    'C00006': [5e-5, 5e-5],
    'C00011': [1e-5, 1e-5]
}


ratio_bound = {
    'C00001': [1, 1],
    'C00009': [5e-3, 5e-3],
    'C00011': [1e-5, 1e-5]
}


ratio = {
    ('C00004', 'C00003'): [5e-4, 5e-1],
    ('C00002', 'C00008'): [2e-1, 10],
    ('C00005', 'C00006'): [23e-2, 1e2],
    ('C00008', 'C00020'): [1, 1]
}

default_params = {
    'ENTRY': 'PathwayID',
    'SKIP': 'False',
    'NAME': 'PathwayName',
    'PH': 7.0,
    'I': 0.1,
    'T': 298.15,
    'C_RANGE': [1e-6, 1e-2],
    'BOUND': default_bound,
    'RATIO_BOUND': ratio_bound,
    'RATIO': ratio
}


# Color configurations for drawpathway.py
color_configs = {}
color_configs['light'] = dict(COFACTOR_SHAPE="ellipse",  # "diamond"
                              OTHER_COFACTOR_COLOR='#B6B6B6',
                              NONCOFACTOR_SHAPE="plaintext",  # "box"
                              NONCOFACTOR_COLOR='transparent',  # "#D2EBEB"
                              REACTION_COLOR="#512DA8",
                              RXN_NODE_COLOR="#323232",
                              EDGE_COLOR="#323232",  # "#505050"
                              BACKGROUND_COLOR="transparent",
                              ALL_FONT_COLOR="black")
# color_configs['light']['colorMapping'] = {
#     'C00002': '#F05456', 'C00008': '#F05456',
#     'C00003': '#FFEB3B', 'C00004': '#FFEB3B',
#     'C00005': '#9dd3ee', 'C00006': '#9dd3ee'
# }

color_configs['light']['colorMapping'] = {
    'C00002': '#F05456', 'C00008': '#FFC000',
    'C00003': '#149B76', 'C00004': '#149B76',
    'C00005': '#2393CB', 'C00006': '#2393CB'
}

# For blue background
color_configs['dark'] = dict(COFACTOR_SHAPE="ellipse",  # "diamond"
                             OTHER_COFACTOR_COLOR="#7F7F7F",
                             NONCOFACTOR_SHAPE="plaintext",  # "box"
                             NONCOFACTOR_COLOR="transparent",  # "#CCFF33"
                             REACTION_COLOR="#FFFF00",
                             EDGE_COLOR="#E5E5E5",  # "#505050"
                             RXN_NODE_COLOR="#E5E5E5",
                             BACKGROUND_COLOR="transparent",
                             ALL_FONT_COLOR="white")

color_configs['dark']['colorMapping'] = {
    'C00002': '#F05456', 'C00008': '#FFC000',
    'C00003': '#149B76', 'C00004': '#149B76',
    'C00005': '#2393CB', 'C00006': '#2393CB'
}
