#!/groups/stark/nikolaus.mandlburger/.conda/envs/tf2_clone_shenzhi/bin/python3.7

### Load arguments
#########

import sys, getopt

def main(argv):
   conv_n = '4'
   dense_n = '2'
   conv_filt = '256'
   filt_size = '7'
   outdir = ''
   try:
      opts, args = getopt.getopt(argv,"hi:v:a:d:o:n:")
   except getopt.GetoptError:
      print('Run_deepSTARR_model.py -i <fold> -v <output variable> -a <architecture> -o <output model ID>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('Run_deepSTARR_model.py -i <fold> -v <output variable> -a <architecture> -o <output model ID>')
         sys.exit()
      elif opt in ("-i"):
         fold = arg
      elif opt in ("-v"):
         output = arg
      elif opt in ("-a"):
         architecture = arg
      elif opt in ("-d"):
         datadir = arg
      elif opt in ("-o"):
         outdir = arg
      elif opt in ("-n"):
         modelname = arg
   if fold=='': sys.exit("fold not found")
   if output=='': sys.exit("variable output not found")
   if architecture=='': sys.exit("architecture not found")
   if outdir=='': sys.exit("datadir not found")
   if outdir=='': sys.exit("outdir not found")
   print('fold ', fold)
   print('variable output ', output)
   print('Model architecture ', architecture)
   print('datadir ', datadir)
   print('outdir ', outdir)
   print("modelname", modelname)
   return fold, output, architecture,datadir, outdir, modelname

if __name__ == "__main__":
   fold, output, architecture,datadir, outdir,modelname  = main(sys.argv[1:])

#########
### Load libraries
#########
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
from sklearn.metrics import roc_curve, auc, average_precision_score
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.linear_model import LinearRegression
from sklearn.utils import class_weight
from matplotlib import pyplot as plt 
import seaborn as sns

import os
import json
import sys

sys.path.append('/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/model_training/accessory_functions/')
from helper import IOHelper, SequenceHelper # from https://github.com/bernardo-de-almeida/Neural_Network_DNA_Demo.git
from Final_model_architectures import *

import random
random.seed(1234)

import tensorflow as tf
from shuffled_controls_functions import run_shuffled_control_predictions

#########
### Load sequences and activity
#########

def prepare_input(fold, set, output):
	# Convert sequences to one-hot encoding matrix
    file_seq = str(fold + "_sequences_" + set + ".fa")
    input_fasta_data_A = IOHelper.get_fastas_from_file(file_seq, uppercase=True)

    # length of first sequence
    sequence_length = len(input_fasta_data_A.sequence.iloc[0])

    # Convert sequence to one hot encoding matrix
    seq_matrix_A = SequenceHelper.do_one_hot_encoding(input_fasta_data_A.sequence, sequence_length,
                                                      SequenceHelper.parse_alpha_to_seq)
    print(seq_matrix_A.shape)
    
    X = np.nan_to_num(seq_matrix_A) # Replace NaN with zero and infinity with large finite numbers
    X_reshaped = X.reshape((X.shape[0], X.shape[1], X.shape[2]))

    Activity = pd.read_table(fold + "_sequences_activity_" + set + ".txt")
    Y = Activity[output]

    weights_dict = {'closed': 1,
   'intermediate': 1,
   'CG_rich_closed': 1,
   'open': 6,
   'pTE_spec_closed': 6,
   'pTE_specopen': 60}

    sample_weights = class_weight.compute_sample_weight(class_weight=weights_dict, y=Activity["class"])
    #sample_weights = class_weight.compute_sample_weight(class_weight='balanced', y=Activity["class"])
    
    print(set)

    return X_reshaped, Y, sample_weights

print("\nLoad sequences\n")

os.chdir(datadir)

X_train, Y_train, importance_weights_train = prepare_input(fold, "Train", output)
X_valid, Y_valid, _ = prepare_input(fold, "Val", output)
X_test, Y_test, _ = prepare_input(fold, "Test", output)

from scipy.stats import spearmanr
def Spearman(y_true, y_pred):
     return ( tf.py_function(spearmanr, [tf.cast(y_pred, tf.float32),
                       tf.cast(y_true, tf.float32)], Tout = tf.float32) )

#########
### Load model architecture
#########

print("\nBuild model: " + architecture +"\n")


if architecture == 'DeepSTARR2':
   get_model=DeepSTARR2

get_model()[0].summary()
get_model()[1] # dictionary

#########
### Model training
#########

print("\nModel training")
#log_dir = model_ID_output + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
#tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)

def train(selected_model, X_train, Y_train, X_valid, Y_valid, params):

	my_history=selected_model.fit(X_train, Y_train,
                                validation_data=(X_valid, Y_valid),
                                sample_weight=importance_weights_train,
                                batch_size=params['batch_size'], epochs=params['epochs'],
                                callbacks=[EarlyStopping(patience=params['early_stop'], monitor="val_loss", restore_best_weights=True),
                                History()])

	return selected_model, my_history


# Model fit
main_model, main_params = get_model()
#main_params["early_stop"]=20
main_model, my_history = train(main_model, X_train, Y_train, X_valid, Y_valid, main_params)

#########
### Save model
#########

print("\nSaving model ...\n")


model_json = main_model.to_json()
with open(os.path.join(outdir,f"{modelname}.json"), "w") as json_file:
    json_file.write(model_json)

main_model.save_weights(os.path.join(outdir,f"{modelname}.h5"))

#history

#with open(os.path.join(outdir,f'{modelname}_history.json'),"w") as historyfile:
#    json.dump(my_history.history, historyfile)



#########
### Model evaluation
#########

print("\nEvaluating model ...\n")

from scipy import stats
from sklearn.metrics import mean_squared_error



y_test_pred = main_model.predict(X_test, batch_size=main_params['batch_size'])
try:
    PCC = stats.pearsonr(Y_test, y_test_pred)[0][0]
except:
    PCC = 5


#plot performance history
fig,ax = plt.subplots(ncols=2,figsize=(18,8))
fig.suptitle(f"Model performance: {modelname}")

ax[0].plot(my_history.history[str('loss')])
ax[0].plot(my_history.history[str('val_loss')])
ax[0].set_title("ATAC_log1p")
ax[0].set_ylabel('Loss MSE')
ax[0].set_xlabel('Epoch')
ax[0].legend(['Train', 'Validation'], loc='upper left')
min_val_loss = min(my_history.history['val_loss'])
ax[0].axvline(x=my_history.history['val_loss'].index(min_val_loss), color='red', linestyle='--')

ax[1].scatter(Y_test, y_test_pred, alpha=0.05, s=1)
ax[1].set_xlabel("ATAC_log1p observed")
ax[1].set_ylabel("ATAC_log1p predicted")
ax[1].set_title(f"PCC: {PCC:.2f}")

# Fit linear regression model
linreg = LinearRegression()
linreg.fit(np.array(Y_test).reshape(-1, 1), y_test_pred)

# Get coefficients and intercept
slope = linreg.coef_[0]
intercept = linreg.intercept_

# Plot the regression line
ax[1].plot(Y_test, slope * np.array(Y_test) + intercept, color='red', label='Regression Line')

# Add legend
ax[1].legend()
#save figure
fig.savefig(os.path.join(outdir,f'{modelname}_performance.png'))

#density plot
g = sns.jointplot(x=Y_test, y=y_test_pred.squeeze(), kind="kde", fill=True, color="red")
#g.ax_marg_x.remove() # remove marginal densities
#g.ax_marg_y.remove() # remove marginal densities
sns.regplot(x=Y_test, y=y_test_pred.squeeze(), scatter=False, color='black', ax=g.ax_joint)
plt.xlabel('Measured accessibility [log1p]')
plt.ylabel('Predicted accessibility [log1p]')
g.savefig(os.path.join(outdir,f'{modelname}_densityplot.png'))

shuffled_seq_predictions_fig = run_shuffled_control_predictions(modelname, outdir, datadir, target_variable=output)
shuffled_seq_predictions_fig.savefig(os.path.join(outdir, f"{modelname}_shuffled_controls.png"))
