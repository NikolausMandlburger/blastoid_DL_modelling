
##################
### Load libraries
##################

import datetime
import tensorflow as tf
import tensorflow.keras.layers as kl
from tensorflow.keras.layers import Conv1D, Conv2D, MaxPooling1D, MaxPooling2D
from tensorflow.keras.layers import Dropout, Reshape, Dense, Activation, Flatten
from tensorflow.keras.layers import BatchNormalization, InputLayer, Input, GlobalAvgPool1D, GlobalMaxPooling1D, LSTM
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam, SGD, RMSprop
from tensorflow.keras.callbacks import EarlyStopping, History, ModelCheckpoint
from tensorflow.keras.regularizers import l1_l2
from tensorflow.keras.utils import plot_model

import pandas as pd
import numpy as np

from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt 
from matplotlib.lines import Line2D
import seaborn as sns

from scipy import stats
from scipy.stats import chisquare
import os

import sys
sys.path.append('/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/model_training/accessory_functions/')
from helper import IOHelper, SequenceHelper # from https://github.com/bernardo-de-almeida/Neural_Network_DNA_Demo.git

import random
random.seed(1234)

import re

import tensorflow as tf
from Final_model_architectures import *


import h5py
from deeplift.dinuc_shuffle import dinuc_shuffle

import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns


from DeepExplainer_accessory_functions import  load_model, prepare_input
from random_sequence_generation import  plot_CG_vs_accsessibility




def run_shuffled_control_predictions(model_name,results_path, datapath, target_variable="ATAC_pTE_log1p", cmap_name = "Set1"):
    
    """
    Run shuffled lable predictions on a certain model and dataset.
    Model name is expected to contain foldnumber (foldNN). This number is extracted and used to load the corresponding Test data set (<fold>_sequences_activity_Test.txt, <fold>_sequences_Test.fa) from the datapath.
    Predictions are run on the test set as well as on the dinucleotide and nucleotide shuffled versions of the test set. 
    The results are plotted and a figure object is returned
    
    Parameters:
    model_name (str): Name of the model. <model_name>.h5 and .json expected to be found in the results dir.
    results_path (str): Path to directory where model h5 and json files are.
    datapath (str): Path to directory where sequence fasta files are located.
    target_variable (str): name of the predicted quantity.
    
    Returns:
    fig (matplotlib figure): The figure object containing the plotted results
    """

    #things that come afterwards
    fold_pattern = r'fold\d{2}'
    fold = re.findall(fold_pattern, model_name)[0]

    #testset fasta/txt path
    test_labels_txt_file = os.path.join(datapath, f"{fold}_sequences_activity_Test.txt")
    test_fasta_file = os.path.join(datapath, f"{fold}_sequences_Test.fa")

    #load model
    model, weights, json = load_model(results_path,model_name)

    #load test data
    X_test = prepare_input("", test_fasta_file, 1001)
    y_test_allvar = pd.read_csv(test_labels_txt_file, sep = "\t")
    y_test = y_test_allvar[target_variable]

    #predict test set
    y_test_pred = model.predict(X_test)

    #test set shuffled dinuc
    X_test_shuf_dinuc = np.array([dinuc_shuffle(X_test[i,:,:]) for i in range(X_test.shape[0])])
    y_test_shuf_dinuc_pred = model.predict(X_test_shuf_dinuc)


    #test set shuffled random
    X_test_shuf_random = np.array([X_test[i,:,:][np.random.choice(1001, 1001, replace=False)] for i in range(X_test.shape[0])])
    y_test_shuf_random_pred = model.predict(X_test_shuf_random)



    #visualise 
    fig, ax = plt.subplots(nrows = 3, ncols = 3, figsize = (22,22))
    
    #create coloring according to class
    cmap = plt.get_cmap(cmap_name)
    
    class_entries = y_test_allvar["class"].unique()
    color_dict = {cl:cmap(i) for i,cl in enumerate(class_entries)}
    class_colors = y_test_allvar["class"].map(color_dict)
    scattercolor_handles = [Line2D([0],[0],marker="o",markerfacecolor= color_dict[cl] ,color="w",label=cl) for cl in class_entries]
    

    #original test set sequences

    plot_CG_vs_accsessibility(X_test,y_test_allvar,ax = ax[0,0], seqset_name="Test set sequences", target_variable = target_variable)
    pcc = stats.pearsonr(y_test, y_test_pred)[0][0]

    ax[0,1].scatter(y_test, y_test_pred , c=class_colors ,s=0.1)
    ax[0,1].set_title(f"true seq performance: {pcc:.2f} pcc")
    ax[0,1].set_xlabel("y test obs")
    ax[0,1].set_ylabel("y test pred")
    ax[0,1].legend(handles = scattercolor_handles)

    sns.kdeplot(y_test_pred, ax = ax[0,2])
    sns.kdeplot(y_test, ax = ax[0,2])
    ax[0,2].set_title("Original test set sequences" )
    ax[0,2].set_xlabel("Distribution of predicted accessibility values" )
    ax[0,2].set_xlim(0,2)

    handles = [
        Line2D([0],[0],color="tab:blue",label="predicted accessibility values"),
        Line2D([0],[0],color="tab:orange",label="original accessibility values")

    ]
    ax[0,2].legend(handles=handles, loc = "upper right")

    #dinucleotide shuffled

    plot_CG_vs_accsessibility(X_test_shuf_dinuc,y_test_allvar,ax = ax[1,0], seqset_name="Test set sequences (dinucleotide shuffled)")
    pcc = stats.pearsonr(y_test, y_test_shuf_dinuc_pred)[0][0]

    ax[1,1].scatter(y_test, y_test_shuf_dinuc_pred,c=class_colors ,s=0.1)
    ax[1,1].set_title(f"dinucleotide shuffled seq performance: {pcc:.2f} pcc")
    ax[1,1].set_xlabel("y test obs")
    ax[1,1].set_ylabel("y test dinuc shuffled pred")
    ax[1,1].legend(handles = scattercolor_handles)

    sns.kdeplot(y_test_shuf_dinuc_pred, ax = ax[1,2])
    sns.kdeplot(y_test, ax = ax[1,2])
    ax[1,2].set_title("Test set sequences dinucleotide shuffled" )
    ax[1,2].set_xlabel("Distribution of predicted accessibility values" )
    ax[1,2].set_xlim(0,2)


    ax[1,2].legend(handles=handles, loc = "upper right")


    #randomly shuffled seqs

    plot_CG_vs_accsessibility(X_test_shuf_random,y_test_allvar,ax = ax[2,0], seqset_name="Test set sequences (randomly shuffled)")
    pcc = stats.pearsonr(y_test, y_test_shuf_random_pred)[0][0]

    ax[2,1].scatter(y_test, y_test_shuf_random_pred,c=class_colors ,s=0.1)
    ax[2,1].set_title(f"randomly shuffled seq performance: {pcc:.2f} pcc")
    ax[2,1].set_xlabel("y test obs")
    ax[2,1].set_ylabel("y test randomly shuffled pred")
    ax[2,1].legend(handles = scattercolor_handles)

    sns.kdeplot(y_test_shuf_random_pred, ax = ax[2,2])
    sns.kdeplot(y_test, ax = ax[2,2])
    ax[2,2].set_title("Test set sequences randomly shuffled" )
    ax[2,2].set_xlabel("Distribution of predicted accessibility values" )
    ax[2,2].set_xlim(0,2)


    ax[2,2].legend(handles=handles, loc = "upper right")
    
    return fig