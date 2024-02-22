import streamlit as st
import pandas as pd
import numpy as np
import re
from PIL import Image
import webbrowser

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions as Reactions

import tensorflow as tf
from tensorflow import keras
from keras.preprocessing import sequence
from keras.utils import pad_sequences
import keras
from keras import backend as K
from keras.models import load_model
import argparse
import h5py
import pdb


seq_rdic = ['A', 'I', 'L', 'V', 'F', 'W', 'Y', 'N', 'C', 'Q', 'M',
    'S', 'T', 'D', 'E', 'R', 'H', 'K', 'G', 'P', 'O', 'U', 'X', 'B', 'Z']
seq_dic = {w: i+1 for i, w in enumerate(seq_rdic)}


@st.cache_data
def encodeSeq(seq, seq_dic):
    if pd.isnull(seq):
        return [0]
    else:
        return [seq_dic[aa] for aa in seq]


@st.cache_resource
def load_modelfile(model_string):
	loaded_model = tf.keras.models.load_model(model_string)
	return loaded_model


@st.cache_resource
def prot_feature_gen_from_str_input(prot_input_str, prot_len=2500):
    Prot_ID = prot_input_str.split(':')[0]
    Prot_seq = prot_input_str.split(':')[1]
    prot_dataframe = pd.DataFrame(
        {'Protein_ID': Prot_ID, 'Sequence': Prot_seq}, index=[0])
    prot_dataframe.set_index('Protein_ID')

    prot_dataframe["encoded_sequence"] = prot_dataframe.Sequence.map(
        lambda a: encodeSeq(a, seq_dic))
    prot_feature = pad_sequences(
        prot_dataframe["encoded_sequence"].values, prot_len)

    return prot_feature, Prot_ID


@st.cache_resource
def mol_feature_gen_from_str_input(mol_str, kegg_id_flag, kegg_df):

	if kegg_id_flag == 1:
		KEGG_ID = mol_str
		kegg_id_loc = kegg_df.index[kegg_df.Compound_ID == KEGG_ID][0]
		KEGG_ID_info = kegg_df.loc[kegg_id_loc]
		KEGG_ID_info_df = KEGG_ID_info.to_frame().T.set_index('Compound_ID')

		final_return = KEGG_ID_info_df
		final_id = KEGG_ID

	else:
		try:
			mol_ID = mol_str.split(':')[0]
			mol_smiles = mol_str.split(':')[1]
			mol = Chem.MolFromSmiles(mol_smiles)
			fp1 = AllChem.GetMorganFingerprintAsBitVect(
			    mol, useChirality=True, radius=2, nBits=2048)
			fp_list = list(np.array(fp1).astype(float))
			fp_str = list(map(str, fp_list))
			mol_fp = '\t'.join(fp_str)

			mol_dict = {}
			mol_dict['Compound_ID'] = mol_ID
			mol_dict['Smiles'] = mol_smiles
			mol_dict['morgan_fp_r2'] = mol_fp

			mol_info_df = pd.DataFrame(mol_dict, index=[0])
			mol_info_df = mol_info_df.set_index('Compound_ID')

			final_return = mol_info_df
			final_id = mol_ID

		except Exception as error:
			print('Something wrong with molecule input string...' + repr(error))

	return final_return, final_id


@st.cache_resource
def act_df_gen_mol_feature(mol_id, prot_id):
	act_df = pd.DataFrame(
	    {'Protein_ID': prot_id, 'Compound_ID': mol_id}, index=[0])

	return act_df


@st.cache_resource
def compound_feature_gen_df_input(act_df, comp_df, comp_len=2048, comp_vec='morgan_fp_r2'):
	act_df = pd.merge(act_df, comp_df, left_on='Compound_ID', right_index=True)
	comp_feature = np.stack(act_df[comp_vec].map(lambda fp: fp.split("\t")))
	comp_feature = comp_feature.astype('float')
	return comp_feature


@st.cache_resource
def model_prediction(compound_feature, enz_feature, model):
    prediction_vals = model.predict([compound_feature, enz_feature])

    return prediction_vals[0][0]


def main():
	graph = tf.compat.v1.get_default_graph()
	ld_model = tf.keras.models.load_model('./data/models/EnzRank/Final_model.model')

	KEGG_compound_read = pd.read_csv('./data/kegg_compound.csv', index_col = 'Compound_ID')
	kegg_df = KEGG_compound_read.reset_index()
	
	st.image('./data/header_image/EnzRank_header.png', use_column_width=True)

	st.subheader('Enzyme-Substrate Activity Predictor ')

	st.subheader('Enzyme sequence')
	st.caption('Please follow the input format show in the text box--> id:Sequence')

	enz_str = st.text_input('', value="A0A4P8WFA8:MTKRVLVTGGAGFLGSHLCERLLSEGHEVICLDNFGSGRRKNIKEFEDHPSFKVNDRDVRISESLPSVDRIYHLASRASPADFTQFPVNIALANTQGTRRLLDQARACDARMVFASTSEVYGDPKVHPQPETYTGNVNIRGARGCYDESKRFGETLTVAYQRKYDVDARTVRIFNTYGPRMRPDDGRVVPTFVTQALRGDDLTIYGDGEQTRSFCYVDDLIEGLISLMRVDNPEHNVYNIGKENERTIKELAYEVLGLTDTESDIVYEPLPEDDPGQRRPDITRAKTELDWEPKISLREGLEDTITYFDN")

    # url = 'https://www.genome.jp/dbget-bin/www_bget?rn:R00801'
    # if st.button('KEformat example'):
    #     webbrowser.open_new_tab(url)

	st.subheader('Substrate ')
	st.caption('Please follow the input format show in the text box--> KEGG id or click the checkbox')

	comp_str = st.text_input('', value="C00149")
	if st.checkbox('If you are entering smiles string along with KEGG ID'):
		add_info = st.text_area('Additional information (id: Smiles):', "C00149:O[C@@H](CC([O-])=O)C([O-])=O")
	else:
		add_info = ''

	if st.button("Predict"):
    # if session_state.button_search:
# 		st.subheader('Enzyme-Substrate activity score')
		with st.spinner('Calculating...'):
			try:
# 				st.write('I am inside')
				prot_feature, prot_id = prot_feature_gen_from_str_input(enz_str)
				if len(add_info) == 0:
					kegg_id_flag = 1
					comp_feature, comp_id = mol_feature_gen_from_str_input(comp_str, kegg_id_flag, kegg_df)
				else:
					kegg_id_flag = 0
					comp_feature, comp_id = mol_feature_gen_from_str_input(add_info, kegg_id_flag, kegg_df)
				
				act_dataframe = act_df_gen_mol_feature(comp_id, prot_id)
# 				st.write(act_dataframe)
				compound_feature = compound_feature_gen_df_input(act_dataframe, comp_feature)
# 				st.write(compound_feature)
            
			except Exception as e:
				st.write('Error somewhere...' + repr(e))
				
# 			st.write(compound_feature)
# 			st.write(prot_feature)
# 			keras.backend.clear_session()

			y = ld_model.predict([compound_feature, prot_feature])
			
			subheaderstring = 'EnzRank Score for '+ prot_id +  '-' + comp_id + ' pair:'
			st.subheader(subheaderstring)
			st.write(str(y[0][0]))
			
if __name__ == '__main__':
    main()
