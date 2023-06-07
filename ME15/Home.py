import streamlit as st

# st.set_page_config(
    # page_title="Hello",
    # page_icon="=K",
# )

#st.set_page_config(layout = "wide")

st.write("# Welcome to novoStoic 2.0")

#st.sidebar.success("Select a demo above.")

st.markdown(
    """
    novoStoic2.0 is an integrated pathway design tool that accounts for 
    thermodynamic considerations and enzyme selection for novel reactions in pathway design 
    
    novoStoic2.0 combines 4 different tools into a single platform
    
    1. optStoic    :- Finds the optimal overall stoichiometry of conversion maximizing yield
    2. novoStoic   :- Pathway design using optimal overall stoichiometry obtained from optStoic including both known and novel transformations learned from known reactions
    3. dGPredictor :- Thermodynamic analysis to safegaurd against infeasible reaction directionalities, especially for reactions with novel molecules
    4. EnzRank     :- Enzyme selection tool for novel reactions
    
    
    Citations: 
    - optStoic : Chowdhury, A., Maranas, C. Designing overall stoichiometric conversions and intervening metabolic reactions. Sci Rep 5, 16009 (2015). https://doi.org/10.1038/srep16009
    - novoStoic: Kumar, A., Wang, L., Ng, C.Y. et al. Pathway design using de novo steps through uncharted biochemical spaces. Nat Commun 9, 184 (2018). https://doi.org/10.1038/s41467-017-02362-x
    - dGPredictor: Wang L, Upadhyay V, Maranas CD (2021) dGPredictor: Automated fragmentation method for metabolic reaction free energy prediction and de novo pathway design. PLoS Comput Biol 17(9): e1009448. https://doi.org/10.1371/journal.pcbi.1009448
    - EnzRank: Upadhyay V, Maranas CD (2023) Rank-ordering of known enzymes as starting points for re-engineering novel substrate activity using a convolutional neural network, Metab Eng (Accepted)    

"""
)
