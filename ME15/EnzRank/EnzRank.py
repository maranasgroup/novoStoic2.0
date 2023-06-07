# import numpy and pandas
import numpy as np
import pandas as pd
import pdb


# import keras modules
import tensorflow as tf
import keras.backend as K
from keras.models import Model
from keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Embedding, Lambda
from keras.layers import Convolution1D, GlobalMaxPooling1D, SpatialDropout1D
from keras.layers import Concatenate
# from tensorflow.keras.optimizers import Adam
from keras.optimizers import Adam
from keras.regularizers import l2,l1
from keras.preprocessing import sequence
from keras.utils import pad_sequences

from sklearn.metrics import precision_recall_curve, auc, roc_curve

physical_devices = tf.config.list_physical_devices('GPU')
## select how many GPUs you have here and uncomment
# try:
#     tf.config.experimental.set_memory_growth(physical_devices[0], True)
#     tf.config.experimental.set_memory_growth(physical_devices[1], True)
#     # tf.config.experimental.set_memory_growth(physical_devices[2], True)
#     # tf.config.experimental.set_memory_growth(physical_devices[3], True)
# except:
#   # Invalid device or cannot modify virtual devices once initialized.
#   pass


strategy = tf.distribute.MirroredStrategy()


# init = tf.compat.v1.global_variables_initializer()
seq_rdic = ['A','I','L','V','F','W','Y','N','C','Q','M','S','T','D','E','R','H','K','G','P','O','U','X','B','Z']
seq_dic = {w: i+1 for i,w in enumerate(seq_rdic)}


def encodeSeq(seq, seq_dic):
    if pd.isnull(seq):
        return [0]
    else:
        return [seq_dic[aa] for aa in seq]

def parse_data(act_dir, mol_dir, protein_dir, with_label=True,
               prot_len=2500, prot_vec="Convolution",
               mol_vec="Convolution", mol_len=2048):

    print("Parsing {0} , {1}, {2} with length {3}, type {4}".format(*[act_dir ,mol_dir, protein_dir, prot_len, prot_vec]))

    protein_col = "Protein_ID"
    mol_col = "Compound_ID"
    col_names = [protein_col, mol_col]
    if with_label: 
        label_col = "Label"
        col_names += [label_col]
    act_df = pd.read_csv(act_dir)
    mol_df = pd.read_csv(mol_dir, index_col="Compound_ID")
    protein_df = pd.read_csv(protein_dir, index_col="Protein_ID")


    if prot_vec == "Convolution":
        protein_df["encoded_sequence"] = protein_df.Sequence.map(lambda a: encodeSeq(a, seq_dic))
    act_df = pd.merge(act_df, protein_df, left_on=protein_col, right_index=True)
    act_df = pd.merge(act_df, mol_df, left_on=mol_col, right_index=True)
    mol_feature = np.stack(act_df[mol_vec].map(lambda fp: fp.split("\t")))
    mol_feature = mol_feature.astype('float')
    # pdb.set_trace()
    if prot_vec=="Convolution":
        protein_feature = pad_sequences(act_df["encoded_sequence"].values, prot_len)
        # pdb.set_trace()
    else:
        protein_feature = np.stack(act_df[prot_vec].map(lambda fp: fp.split("\t")))
    if with_label:
        label = act_df[label_col].values
        print("\tPositive data : %d" %(sum(act_df[label_col])))
        print("\tNegative data : %d" %(act_df.shape[0] - sum(act_df[label_col])))
        return {"protein_feature": protein_feature, "mol_feature": mol_feature, "label": label}
    else:
        return {"protein_feature": protein_feature, "mol_feature": mol_feature}


class mol_Enz_Prediction(object):
    
    def PLayer(self, size, filters, activation, initializer, regularizer_param):
        def f(input):
            model_p = Convolution1D(filters=filters, kernel_size=size, padding='same', kernel_initializer=initializer, kernel_regularizer=l2(regularizer_param))(input)
            model_p = BatchNormalization()(model_p)
            model_p = Activation(activation)(model_p)
            return GlobalMaxPooling1D()(model_p)
        return f

    def modelv(self, dropout, mol_layers, protein_strides, filters, fc_layers, prot_vec=False, prot_len=2500,
               activation='relu', protein_layers=None, initializer="glorot_normal", mol_len=2048, mol_vec="ECFP4"):
        def return_tuple(value):
            if type(value) is int:
               return [value]
            else:
               return tuple(value)

        regularizer_param = 0.001
        input_d = Input(shape=(mol_len,))
        input_p = Input(shape=(prot_len,))
        params_dic = {"kernel_initializer": initializer,
                      # "activity_regularizer": l2(regularizer_param),
                      "kernel_regularizer": l2(regularizer_param),
        }
        input_layer_d = input_d
        if mol_layers is not None:
            mol_layers = return_tuple(mol_layers)
            for layer_size in mol_layers:
                model_d = Dense(layer_size, **params_dic)(input_layer_d)
                model_d = BatchNormalization()(model_d)
                model_d = Activation(activation)(model_d)
                model_d = Dropout(dropout)(model_d)
                input_layer_d = model_d

        if prot_vec == "Convolution":
            model_p = Embedding(26,20, embeddings_initializer=initializer,embeddings_regularizer=l2(regularizer_param))(input_p)
            model_p = SpatialDropout1D(0.2)(model_p)
            model_ps = [self.PLayer(stride_size, filters, activation, initializer, regularizer_param)(model_p) for stride_size in protein_strides]
            if len(model_ps)!=1:
                model_p = Concatenate(axis=1)(model_ps)
            else:
                model_p = model_ps[0]
        else:
            model_p = input_p
        if protein_layers:
            input_layer_p = model_p
            protein_layers = return_tuple(protein_layers)
            for protein_layer in protein_layers:
                model_p = Dense(protein_layer, **params_dic)(input_layer_p)
                model_p = BatchNormalization()(model_p)
                model_p = Activation(activation)(model_p)
                model_p = Dropout(dropout)(model_p)
                input_layer_p = model_p

        model_t = Concatenate(axis=1)([model_d,model_p])
        #input_dim = filters*len(protein_strides) + mol_layers[-1]
        if fc_layers is not None:
            fc_layers = return_tuple(fc_layers)
            for fc_layer in fc_layers:
                model_t = Dense(units=fc_layer,#, input_dim=input_dim,
                                **params_dic)(model_t)
                model_t = BatchNormalization()(model_t)
                model_t = Activation(activation)(model_t)
                # model_t = Dropout(dropout)(model_t)
                input_dim = fc_layer
        model_t = Dense(1, activation='tanh', activity_regularizer=l2(regularizer_param),**params_dic)(model_t)
        model_t = Lambda(lambda x: (x+1.)/2.)(model_t)

        model_f = Model(inputs=[input_d, input_p], outputs = model_t)

        return model_f

    def __init__(self, dropout=0.2, mol_layers=(1024,512), protein_strides = (10,15,20,25), filters=64,
                 learning_rate=1e-3, decay=0.0, fc_layers=None, prot_vec=None, prot_len=2500, activation="relu",
                 mol_len=2048, mol_vec="ECFP4", protein_layers=None):
        self.__dropout = dropout
        self.__mols_layer = mol_layers
        self.__protein_strides = protein_strides
        self.__filters = filters
        self.__fc_layers = fc_layers
        self.__learning_rate = learning_rate
        self.__prot_vec = prot_vec
        self.__prot_len = prot_len
        self.__mol_vec = mol_vec
        self.__mol_len = mol_len
        self.__activation = activation
        self.__protein_layers = protein_layers
        self.__decay = decay
        self.__model_t = self.modelv(self.__dropout, self.__mols_layer, self.__protein_strides,
                                     self.__filters, self.__fc_layers, prot_vec=self.__prot_vec,
                                     prot_len=self.__prot_len, activation=self.__activation,
                                     protein_layers=self.__protein_layers, mol_vec=self.__mol_vec,
                                     mol_len=self.__mol_len)

        opt = Adam(lr=learning_rate, decay=self.__decay)
        self.__model_t.compile(optimizer=opt, loss='binary_crossentropy', metrics=['accuracy'])
        # K.get_session().run(tf.compat.v1.global_variables_initializer())
#         tf.compat.v1.global_variables_initializer()
#         K.get_session().run(tf.initializer.global_variables())
#         K.get_session().run(tf.global_variables_initializer())
#         K.get_session().run(tf.initialize_all_variables())
        

        

    def fit(self, mol_feature, protein_feature, label, n_epoch=10, batch_size=32):
        for _ in range(n_epoch):
            history = self.__model_t.fit([mol_feature,protein_feature],label, epochs=_+1, batch_size=batch_size, shuffle=True, verbose=1,initial_epoch=_)
        return self.__model_t
    
    def summary(self):
        self.__model_t.summary()
    
    def validation(self, mol_feature, protein_feature, label, output_file=None, n_epoch=10, batch_size=32, **kwargs):

        if output_file:
            param_tuple = pd.MultiIndex.from_tuples([("parameter", param) for param in ["window_sizes", "mol_layers", "fc_layers", "learning_rate"]])
            result_df = pd.DataFrame(data = [[self.__protein_strides, self.__mols_layer, self.__fc_layers, self.__learning_rate]]*n_epoch, columns=param_tuple)
            result_df["epoch"] = range(1,n_epoch+1)
        result_dic = {dataset: {"AUC":[], "AUPR": [], "opt_threshold(AUPR)":[], "opt_threshold(AUC)":[] }for dataset in kwargs}

        for _ in range(n_epoch):
            history = self.__model_t.fit([mol_feature,protein_feature],label,
                                         epochs=_+1, batch_size=batch_size, shuffle=True, verbose=1, initial_epoch=_)
            for dataset in kwargs:
                print("\tPredction of " + dataset)
                test_p = kwargs[dataset]["protein_feature"]
                test_d = kwargs[dataset]["mol_feature"]
                test_label = kwargs[dataset]["label"]
                prediction = self.__model_t.predict([test_d,test_p])
                fpr, tpr, thresholds_AUC = roc_curve(test_label, prediction)
                AUC = auc(fpr, tpr)
                precision, recall, thresholds = precision_recall_curve(test_label,prediction)
                distance = (1-fpr)**2+(1-tpr)**2
                EERs = (1-recall)/(1-precision)
                positive = sum(test_label)
                negative = test_label.shape[0]-positive
                ratio = negative/positive
                opt_t_AUC = thresholds_AUC[np.argmin(distance)]
                opt_t_AUPR = thresholds[np.argmin(np.abs(EERs-ratio))]
                AUPR = auc(recall,precision)
                print("\tArea Under ROC Curve(AUC): %0.3f" % AUC)
                print("\tArea Under PR Curve(AUPR): %0.3f" % AUPR)
                print("\tOptimal threshold(AUC)   : %0.3f " % opt_t_AUC)
                print("\tOptimal threshold(AUPR)  : %0.3f" % opt_t_AUPR)
                print("=================================================")
                result_dic[dataset]["AUC"].append(AUC)
                result_dic[dataset]["AUPR"].append(AUPR)
                result_dic[dataset]["opt_threshold(AUC)"].append(opt_t_AUC)
                result_dic[dataset]["opt_threshold(AUPR)"].append(opt_t_AUPR)
        if output_file:
            for dataset in kwargs:
                result_df[dataset, "AUC"] = result_dic[dataset]["AUC"]
                result_df[dataset, "AUPR"] = result_dic[dataset]["AUPR"]
                result_df[dataset, "opt_threshold(AUC)"] = result_dic[dataset]["opt_threshold(AUC)"]
                result_df[dataset, "opt_threshold(AUPR)"] = result_dic[dataset]["opt_threshold(AUPR)"]
            print("save to " + output_file)
            print(result_df)
            result_df.to_csv(output_file, index=False)

    def predict(self, **kwargs):
        results_dic = {}
        for dataset in kwargs:
            result_dic = {}
            test_p = kwargs[dataset]["protein_feature"]
            test_d = kwargs[dataset]["mol_feature"]
            result_dic["label"] = kwargs[dataset]["label"]
            result_dic["predicted"] = self.__model_t.predict([test_d, test_p])
            results_dic[dataset] = result_dic
        return results_dic
    
    def save(self, output_file):
        self.__model_t.save(output_file)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="""
    This Python script is used to train, validate, test deep learning model for prediction of mol-Enz interaction (act)\n
    Deep learning model will be built by Keras with tensorflow.\n
    You can set almost hyper-parameters as you want, See below parameter description\n
    act, mol and protein data must be written as csv file format. And feature should be tab-delimited format for script to parse data.\n
    Basically, this script builds convolutional neural network on sequence.\n
    If you don't want convolutional neural network but traditional dense layers on provide protein feature, specify type of feature and feature length.\n
    \n
    requirement\n
    ============================\n
    tensorflow > 1.0\n
    keras > 2.0\n
    numpy\n
    pandas\n
    scikit-learn\n
    ============================\n
    \n
    contact : dlsrnsladlek@gist.ac.kr\n
    """)
    # train_params
    parser.add_argument("act_dir", help="Training act information [mol, Enz, label]")
    parser.add_argument("mol_dir", help="Training mol information [mol, SMILES,[feature_name, ..]]")
    parser.add_argument("protein_dir", help="Training protein information [protein, seq, [feature_name]]")
    # test_params
    parser.add_argument("--test-name", '-n', help="Name of test data sets", nargs="*")
    parser.add_argument("--test-act-dir", "-i", help="Test act [mol, Enz, [label]]", nargs="*")
    parser.add_argument("--test-mol-dir", "-d", help="Test mol information [mol, SMILES,[feature_name, ..]]", nargs="*")
    parser.add_argument("--test-protein-dir", '-t', help="Test Protein information [protein, seq, [feature_name]]", nargs="*")
    parser.add_argument("--with-label", "-W", help="Existence of label information in test act", action="store_true")
    # structure_params
    parser.add_argument("--window-sizes", '-w', help="Window sizes for model (only works for Convolution)", default=None, nargs="*", type=int)
    parser.add_argument("--protein-layers","-p", help="Dense layers for protein", default=None, nargs="*", type=int)
    parser.add_argument("--mol-layers", '-c', help="Dense layers for mols", default=None, nargs="*", type=int)
    parser.add_argument("--fc-layers", '-f', help="Dense layers for concatenated layers of mol and Enz layer", default=None, nargs="*", type=int)
    # training_params
    parser.add_argument("--learning-rate", '-r', help="Learning late for training", default=1e-4, type=float)
    parser.add_argument("--n-epoch", '-e', help="The number of epochs for training or validation", type=int, default=10)
    # type_params
    parser.add_argument("--prot-vec", "-v", help="Type of protein feature, if Convolution, it will execute conlvolution on sequeunce", type=str, default="Convolution")
    parser.add_argument("--prot-len", "-l", help="Protein vector length", default=2500, type=int)
    parser.add_argument("--mol-vec", "-V", help="Type of mol feature", type=str, default="morgan_fp")
    parser.add_argument("--mol-len", "-L", help="mol vector length", default=2048, type=int)
    # the other hyper-parameters
    parser.add_argument("--activation", "-a", help='Activation function of model', type=str)
    parser.add_argument("--dropout", "-D", help="Dropout ratio", default=0.2, type=float)
    parser.add_argument("--n-filters", "-F", help="Number of filters for convolution layer, only works for Convolution", default=64, type=int)
    parser.add_argument("--batch-size", "-b", help="Batch size", default=32, type=int)
    parser.add_argument("--decay", "-y", help="Learning rate decay", default=0.0, type=float)
    # mode_params
    parser.add_argument("--validation", help="Excute validation with independent data, will give AUC and AUPR (No prediction result)", action="store_true")
    parser.add_argument("--predict", help="Predict interactions of independent test set", action="store_true")
    # output_params
    parser.add_argument("--save-model", "-m", help="save model", type=str)
    parser.add_argument("--output", "-o", help="Prediction output", type=str)

    args = parser.parse_args()
    # train data
    train_dic = {
        "act_dir": args.act_dir,
        "mol_dir": args.mol_dir,
        "protein_dir": args.protein_dir,
        "with_label": True
    }
    # create dictionary of test_data
    test_names = args.test_name
    tests = args.test_act_dir
    test_proteins = args.test_protein_dir
    test_mols = args.test_mol_dir
    test_sets = zip(test_names, tests, test_mols, test_proteins)
    output_file = args.output
    # model_structure variables
    mol_layers = args.mol_layers
    window_sizes = args.window_sizes
    if window_sizes==0:
        window_sizes = None
    protein_layers = args.protein_layers
    fc_layers = args.fc_layers
    # training parameter
    train_params = {
        "n_epoch": args.n_epoch,
        "batch_size": args.batch_size,
    }
    # type parameter
    type_params = {
        "prot_vec": args.prot_vec,
        "prot_len": args.prot_len,
        "mol_vec": args.mol_vec,
        "mol_len": args.mol_len,
    }
    # model parameter
    model_params = {
        "mol_layers": mol_layers,
        "protein_strides": window_sizes,
        "protein_layers": protein_layers,
        "fc_layers": fc_layers,
        "learning_rate": args.learning_rate,
        "decay": args.decay,
        "activation": args.activation,
        "filters": args.n_filters,
        "dropout": args.dropout
    }

    model_params.update(type_params)
    print("\tmodel parameters summary\t")
    print("=====================================================")
    for key in model_params.keys():
        print("{:20s} : {:10s}".format(key, str(model_params[key])))
    print("=====================================================")

    act_prediction_model = mol_Enz_Prediction(**model_params)
    act_prediction_model.summary()

    # read and parse training and test data
    train_dic.update(type_params)
    train_dic = parse_data(**train_dic)
    test_dic = {test_name: parse_data(test_act, test_mol, test_protein, with_label=True, **type_params)
                for test_name, test_act, test_mol, test_protein in test_sets}

    # validation mode
    if args.validation:
        validation_params = {}
        validation_params.update(train_params)
        validation_params["output_file"] = output_file
        print("\tvalidation summary\t")
        print("=====================================================")
        for key in validation_params.keys():
            print("{:20s} : {:10s}".format(key, str(validation_params[key])))
        print("=====================================================")
        validation_params.update(train_dic)
        validation_params.update(test_dic)
        act_prediction_model.validation(**validation_params)
    # prediction mode
    elif args.predict:
        print("prediction")
        train_dic.update(train_params)
        act_prediction_model.fit(**train_dic)
        test_predicted = act_prediction_model.predict(**test_dic)
        result_df = pd.DataFrame()
        result_columns = []
        for dataset in test_predicted:
            temp_df = pd.DataFrame()
            value = test_predicted[dataset]["predicted"]
            value = np.squeeze(value)
            print(dataset+str(value.shape))
            temp_df[dataset,'predicted'] = value
            temp_df[dataset, 'label'] = np.squeeze(test_predicted[dataset]['label'])
            result_df = pd.concat([result_df, temp_df], ignore_index=True, axis=1)
            result_columns.append((dataset, "predicted"))
            result_columns.append((dataset, "label"))
        result_df.columns = pd.MultiIndex.from_tuples(result_columns)
        print("save to %s"%output_file)
        result_df.to_csv(output_file, index=False)

    # save trained model
    if args.save_model:
        act_prediction_model.save(args.save_model)
    exit()
