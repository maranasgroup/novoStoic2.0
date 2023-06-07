# EnzRank

CNN based model for enzyme-substrate activity prediction using enzyme sequence and substrate structural information. EnzRank is used for rank-ordering novel enzyme-substrate activities identified in de novo biosynthesis pathway design and enzyme re-engineering projects. 

User interface is available at "https://huggingface.co/spaces/vuu10/EnzRank"

## Requirements
- Tensorflow 2
- rdkit
- pandas
- tqdm
- numpy 
- keras 
- scikit-learn
- streamlit

## Steps to run the tool locally

- create a conda environment using: `conda create --prefix MLenv`
- activate the created environment using: `conda activate MLenv`
- install packages using requirement.txt file: `pip install -r requirement.txt`
- install rdkit using: `pip install rdkit` 
- install streamlit using: `pip install streamlit`

## The EnzRank tool is tested on Linux-based system

Three data inputs are required- proteins sequence, substrate smiles, and enzyme-substrate pairs (Look at CNN_data or CNN_split_data folder for the required input files)- input data must follow the same format

For performance prediction on different splits of data (CNN_data_split folder) run. 

`python predict_with_model.py ./CNN_results/CNN_results_1/model.model -n predict -i ./CNN_data_split/CNN_data_1/test_dataset/test_act.csv -d ./CNN_data_split/CNN_data_1/test_dataset/test_compound.csv -t ./CNN_data_split/CNN_data_1/test_dataset/test_protein.csv -v Convolution -l 2500 -V morgan_fp_r2 -L 2048 -W -o ./CNN_results/CNN_results_1/output_predictions.csv`

-- Please use "predict_with_model_less.py" if you have less than 50 enzyme-substrate 

This will output "output_prediction.csv" file in the CNN_results folder

To launch the streamlit based graphical user interface: 

move to 'Streamlit' directory using- `cd Streamlit` 
launch interface using command: `streamlit run main.py` 


## Usage 
```
  usage: EnzRank.py [-h] [--test-name [TEST_NAME [TEST_NAME ...]]]
                     [--test-act-dir [TEST_act_DIR [TEST_act_DIR ...]]]
                     [--test-mol-dir [TEST_mol_DIR [TEST_mol_DIR ...]]]
                     [--test-enz-dir [TEST_enz_DIR [TEST_enz_DIR ...]]]
                     [--with-label WITH_LABEL]
                     [--window-sizes [WINDOW_SIZES [WINDOW_SIZES ...]]]
                     [--enz-layers [enz_LAYERS [enz_LAYERS ...]]]
                     [--mol-layers [mol_LAYERS [mol_LAYERS ...]]]
                     [--fc-layers [FC_LAYERS [FC_LAYERS ...]]]
                     [--learning-rate LEARNING_RATE] [--n-epoch N_EPOCH]
                     [--prot-vec PROT_VEC] [--prot-len PROT_LEN]
                     [--mol-vec mol_VEC] [--mol-len mol_LEN]
                     [--activation ACTIVATION] [--dropout DROPOUT]
                     [--n-filters N_FILTERS] [--batch-size BATCH_SIZE]
                     [--decay DECAY] [--validation] [--predict]
                     [--save-model SAVE_MODEL] [--output OUTPUT]
                     act_dir mol_dir enz_dir
```

## Parameter specification

### Positional arguments
```
    act_dir               Training act information [mol, target, label]
    mol_dir               Training mol information [mol, SMILES,[feature_name,..]]
    enz_dir               Training enz information [enz, seq]
```
### Optional arguments

```
    --validation          Excute validation with independent data, will give AUC and AUPR (No prediction result)
    --predict             Predict interactions of independent test set
```

EnzRank script has two mode, validation and predict mode.

In validation mode, performances (AUC, AUPR, threshold for AUC and AUPR) on each step and selected hyperparameters are recorded.

In test step, prediction results for given test dataset will be reported after training.

```
    --test-name [TEST_NAME [TEST_NAME ...]], -n [TEST_NAME [TEST_NAME ...]]                     Name of test data sets
    --test-act-dir [TEST_act_DIR [TEST_act_DIR ...]], -i [TEST_act_DIR [TEST_act_DIR ...]]      Test act [mol, target, [label]]
    --test-mol-dir [TEST_mol_DIR [TEST_mol_DIR ...]], -d [TEST_mol_DIR [TEST_mol_DIR ...]]      Test mol information [mol, SMILES,[feature_name,..]]
    --test-enz-dir [TEST_enz_DIR [TEST_enz_DIR ...]], -t [TEST_enz_DIR [TEST_enz_DIR ...]]      Test enz information [enz, seq]
    --with-label WITH_LABEL, -W WITH_LABEL                                                      Existence of label information in test act
```
You can input multiple datasets for validation or test with argument specifier.

In addition to act information file, mol information file and target enz information file, you need to name of validation or test datasets.

For test dataset, you can inform that test dataset has label or not with `-W` value

```
    --n-epoch N_EPOCH, -e N_EPOCH     The number of epochs for training or validation
    
    --prot-vec      PROT_VEC, -v Convolution     Convolution, can change this to use difference protein feature instead of convolution
    --prot-len      PROT_LEN, -l PROT_LEN        enz vector length
    --mol-vec       mol_VEC,  -V mol_VEC         Type of mol feature (morgan_fp_r2)
    --mol-len       mol_LEN,  -L mol_LEN         mol vector length
    --window-sizes  [WINDOW_SIZES [WINDOW_SIZES ...]], -w [WINDOW_SIZES [WINDOW_SIZES ...]]   Window sizes for model (only works for Convolution)
    --enz-layers    [enz_LAYERS [enz_LAYERS ...]], -p [enz_LAYERS [enz_LAYERS ...]]           Dense layers for enz
    --mol-layers    [mol_LAYERS [mol_LAYERS ...]], -c [mol_LAYERS [mol_LAYERS ...]]           Dense layers for mols
    --fc-layers     [FC_LAYERS [FC_LAYERS ...]], -f [FC_LAYERS [FC_LAYERS ...]]               Dense layers for concatenated layers of mol and enzyme layer
    --n-filters     N_FILTERS, -F N_FILTERS                                                   Number of filters for convolution layer, only works for Convolution

    --activation    ACTIVATION, -a ACTIVATION           Activation function of model
    --dropout       DROPOUT, -D DROPOUT                 Dropout ratio
    --batch-size    BATCH_SIZE, -b BATCH_SIZE           Batch size
    --learning-rate LEARNING_RATE, -r LEARNING_RATE     Learning late for training
    --decay         DECAY, -y DECAY                     Learning rate decay
    --save-model    SAVE_MODEL, -m SAVE_MODEL           save model
    --output        OUTPUT, -o OUTPUT                   Prediction output
```

### Training and validation of model

For an example, if you have training dataset, `CNN_data_split/split_0/training_dataset/training_act.csv`, `CNN_data_split/split_0/training_dataset/training_compound.csv` and `CNN_data_split/split_0/training_dataset/training_enz.csv` with right specification.
You can validate model with validation dataset, `CNN_data_split/split_0/validation_dataset/validation_act.csv`, `CNN_data_split/split_0/validation_dataset/validation_compound.csv` and `CNN_data_split/split_0/validation_dataset/validation_enz.csv` by using this command line.

`python EnzRank.py ./CNN_data_split/split_0/training_dataset/training_act.csv ./CNN_data_split/split_0/training_dataset/training_compound.csv ./CNN_data_split/split_0/training_dataset/training_enz.csv --validation -n validation_dataset -i ./CNN_data_split/split_0/validation_dataset/validation_act.csv -d ./CNN_data_split/split_0/validation_dataset/validation_compound.csv -t ./CNN_data_split/split_0/validation_dataset/validation_enz.csv -W -c 512 128 -w 10 15 20 25 30 -p 128 -f 128 -r 0.0001 -n 30 -v Convolution -l 2500 -V morgan_fp_r2 -L 2048 -D 0 -a elu -F 128 -b 32 -y 0.0001 -o ./validation_output.csv -m ./model.model -e 1`

This command will train model with given hyper-parameters for 1 epoch. (because of -e 1).
And resulting in validation result `./validation_output.csv`, and corresponding model `./model.model`

