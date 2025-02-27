
#!/groups/stark/nikolaus.mandlburger/.conda/envs/tf2_clone_shenzhi/bin/python3.7

import sys, getopt
import tensorflow as tf
#tf.config.threading.set_inter_op_parallelism_threads(1)
#tf.config.threading.set_intra_op_parallelism_threads(1)
#tf.config.threading.enable_op_determinism()
#tf.compat.v1.config.threading.enable_op_determinism()
#from tensorflow.compat.v1.keras.backend import get_session
tf.compat.v1.disable_v2_behavior()
tf.compat.v1.disable_eager_execution()
import keras
import keras.layers as kl
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers import BatchNormalization, InputLayer, Input
from keras import models
from keras.models import Sequential, Model
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History, ModelCheckpoint
from keras import backend as K

import pandas as pd
import numpy as np
import os
import random
import matplotlib.pyplot as plt 
import json
import pkbar
import datetime

from scipy.stats.mstats import spearmanr
from scipy.stats.mstats import pearsonr
from scipy import stats
from sklearn.metrics import mean_squared_error
import deeplift
from keras.models import model_from_json
from deeplift.dinuc_shuffle import dinuc_shuffle
import shap # forked from https://github.com/AvantiShri/shap/blob/master/shap/explainers/deep/deep_tf.py
from deeplift.visualization import viz_sequence
import h5py

sys.path.append('/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/model_training/accessory_functions/')
from helper import IOHelper, SequenceHelper # from https://github.com/const-ae/Neural_Network_DNA_Demo

### other parameters
# number of background sequences to take an expectation over
ns=1000
# number of dinucleotide shuffled sequences per sequence as background
dinuc_shuffle_n=100
input_length=1001


def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1


def prepare_input(set, data_input_path, target_length):
    print(set)
    # Convert sequences to one-hot encoding matrix
    file_seq = str(data_input_path)
    input_fasta_data_A = IOHelper.get_fastas_from_file(file_seq, uppercase=True)

    target_length = int(target_length)
    if len(input_fasta_data_A.sequence.iloc[0]) < target_length:
        print("Pad sequence at 3' end because train data is augmented with insertions")
        print("before padding length = " + str(len(input_fasta_data_A.sequence.iloc[0])))
        padded_sequences = []
        for seq in input_fasta_data_A.sequence:
            diff = target_length - len(seq)
            if diff != 0:
                n_seq = seq + ''.join(np.random.choice(('C','G','T','A'), diff ))
            else:
                seq = seq
            padded_sequences.append(n_seq)
        input_fasta_data_A.sequence = padded_sequences
        print("after padding length = " + str(len(input_fasta_data_A.sequence.iloc[0])))

    # length of first sequence
    sequence_length = len(input_fasta_data_A.sequence.iloc[0])

    # Convert sequence to one hot encoding matrix
    seq_matrix_A = SequenceHelper.do_one_hot_encoding(input_fasta_data_A.sequence, sequence_length, SequenceHelper.parse_alpha_to_seq)
    print(seq_matrix_A.shape)
    
    X = np.nan_to_num(seq_matrix_A) # Replace NaN with zero and infinity with large finite numbers
    #X_reshaped = X.reshape((X.shape[0], X.shape[1], X.shape[2]))
    print("X", X.shape)
    return X #X_reshaped


def load_model(model_path, model):
    #import deeplift
    #from tensorflow.keras.models import model_from_json
    keras_model_weights = model_path + model + '.h5'
    keras_model_json = model_path + model + '.json'
    keras_model = model_from_json(open(keras_model_json).read())
    keras_model.load_weights(keras_model_weights)
    print(keras_model.summary())
    return keras_model, keras_model_weights, keras_model_json

#from deeplift.dinuc_shuffle import dinuc_shuffle
#import numpy as np


### deepExplainer functions
def dinuc_shuffle_several_times(list_containing_input_modes_for_an_example,seed=1234):
  assert len(list_containing_input_modes_for_an_example)==1
  onehot_seq = list_containing_input_modes_for_an_example[0]
  rng = np.random.RandomState(seed)
  to_return = np.array([dinuc_shuffle(onehot_seq, rng=rng) for i in range(dinuc_shuffle_n)])
  return [to_return] #wrap in list for compatibility with multiple modes

# get hypothetical scores also
def combine_mult_and_diffref(mult, orig_inp, bg_data):
    assert len(orig_inp)==1
    projected_hypothetical_contribs = np.zeros_like(bg_data[0]).astype("float")
    assert len(orig_inp[0].shape)==2
    #At each position in the input sequence, we iterate over the one-hot encoding
    # possibilities (eg: for genomic sequence, this is ACGT i.e.
    # 1000, 0100, 0010 and 0001) and compute the hypothetical 
    # difference-from-reference in each case. We then multiply the hypothetical
    # differences-from-reference with the multipliers to get the hypothetical contributions.
    #For each of the one-hot encoding possibilities,
    # the hypothetical contributions are then summed across the ACGT axis to estimate
    # the total hypothetical contribution of each position. This per-position hypothetical
    # contribution is then assigned ("projected") onto whichever base was present in the
    # hypothetical sequence.
    #The reason this is a fast estimate of what the importance scores *would* look
    # like if different bases were present in the underlying sequence is that
    # the multipliers are computed once using the original sequence, and are not
    # computed again for each hypothetical sequence.
    for i in range(orig_inp[0].shape[-1]):
        hypothetical_input = np.zeros_like(orig_inp[0]).astype("float")
        hypothetical_input[:,i] = 1.0
        hypothetical_difference_from_reference = (hypothetical_input[None,:,:]-bg_data[0])
        hypothetical_contribs = hypothetical_difference_from_reference*mult[0]
        projected_hypothetical_contribs[:,:,i] = np.sum(hypothetical_contribs,axis=-1) 
    return [np.mean(projected_hypothetical_contribs,axis=0)]


def my_deepExplainer(model, one_hot, bg):
    #import shap # forked from Avanti https://github.com/AvantiShri/shap/blob/master/shap/explainers/deep/deep_tf.py
    #import numpy as np
    
    # output layer
    out_layer=-1

    # output layer
    if bg=="random":
        np.random.seed(seed=1111)
        background = one_hot[np.random.choice(one_hot.shape[0], ns, replace=False)]
        explainer = shap.DeepExplainer((model.layers[0].input, model.layers[out_layer].output),
          data=background,
          combine_mult_and_diffref=combine_mult_and_diffref)
    if bg=="dinuc_shuffle":
        explainer = shap.DeepExplainer((model.layers[0].input, model.layers[out_layer].output),
          data=dinuc_shuffle_several_times,
          combine_mult_and_diffref=combine_mult_and_diffref)
        
    # running on all sequences
    shap_values_hypothetical = explainer.shap_values(one_hot)
    
    # normalising contribution scores
    # sum the deeplift importance scores across the ACGT axis (different nucleotides at the same position)
    # and “project” that summed importance onto whichever base is actually present at that position
    shap_values_contribution=shap_values_hypothetical[0]*one_hot

    return shap_values_hypothetical[0], shap_values_contribution