#!/groups/stark/nikolaus.mandlburger/.conda/envs/tf2_clone_shenzhi/bin/python3.7

#!/groups/stark/nikolaus.mandlburger/.conda/envs/tf2_clone_shenzhi_tangermeme/bin/python3.7



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
from matplotlib import pyplot as plt 
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from keras.models import model_from_json

import os
import glob
import json
import sys
import re
import argparse

import random
random.seed(1234)

import tensorflow as tf


import re
sys.path.append('/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/model_training/accessory_functions/')
from Final_model_architectures import *
from DeepExplainer_accessory_functions import  load_model, prepare_input


sys.path.append('/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/primate_sequence_variant_prediction/')
from candidate_predictions_accessory_functions import species_dict, species_names_trivial, species_names_scientific, species_info_full, color_dict



#parse arguments (argpars): chromosome, model_dir, res_dir, ROI_name
parser = argparse.ArgumentParser(prog='ProgramName', description="Make predictions for a candidate enhancer in different species")
parser.add_argument("-c","--chromosome",type=str, help="The chromosome of the candidate enhancer")
parser.add_argument("-m","--model_dir",type=str, help="The directory where the model is stored")
parser.add_argument("-r","--res_dir",type=str, help="The directory where the results are stored")
parser.add_argument("-n","--ROI_name",type=str, help="The name of the candidate enhancer")
args=parser.parse_args()

chromosome = args.chromosome
model_dir = args.model_dir
res_dir = args.res_dir
ROI_name = args.ROI_name


chrom_to_fold_dict = {"chr6":"fold06",
                      "chr2":"fold02",
                      "chr10":"fold10"}

fold=chrom_to_fold_dict[chromosome]




ROI_name_fa_list = glob.glob(os.path.join(res_dir, f"{ROI_name}*_shift_*.fa"))
X_info_all_species_all_shifts_list = []
modelname_list = [os.path.basename(p).replace(".h5","") for p in glob.glob(model_dir+f"pTE_{fold}_rep_*.h5")]

for i,model_name in enumerate(modelname_list):
    #load_model 
    rep = re.search(r'rep_(\d+)', model_name).group(1)
    model, weights, json = load_model(model_dir,model_name)
    X_info_all_species = species_info_full.copy()
    
    for fa_name in ROI_name_fa_list:
        print(fa_name)
        shift_match = re.search(r"shift_(-?\d+)", fa_name)
        shift_value = int(shift_match.group(1))
        
        #load data corresponding to current 
        X = prepare_input("",fa_name,1001)
        X_info = pd.read_csv(fa_name.replace(".fa",".tsv"), sep = "\t")
        y_pred = model.predict(X)

        #adding strand and prediction information to the dataframe
        X_info["pred_pTE_accessibility"]=y_pred
        X_info["strand"] = [re.search(r'_(fwd|rv)_', e).group(1) for e in X_info["seqname"]]


        #merge predictions into full species df
        species_info_pred = pd.merge(species_info_full, X_info, on=['species','species_trivial', 'strand'], how='left')
        species_info_pred["strand"].loc[species_info_pred["pred_pTE_accessibility"].isna()]="no_seq"
        #species_info_pred["pred_pTE_accessibility"].loc[species_info_pred["pred_pTE_accessibility"].isna()]=0

        #save predictions for current replicate
        X_info_all_species[f"pred_pTE_accessibility_rep_{rep}_shift_{shift_value}"]=np.array(species_info_pred["pred_pTE_accessibility"])
        print(f"rep {rep} shift {shift_value} finished")
    X_info_all_species_all_shifts_list.append(X_info_all_species)
    #if i==2:
        #break

X_info_all_species_all_shifts_df = pd.concat([df.set_index(["species","species_trivial","strand"]) for df in X_info_all_species_all_shifts_list], axis = 1).reset_index()



#create df for plotting
X_info_all_species_all_shifts_df_fwd = X_info_all_species_all_shifts_df.loc[X_info_all_species_all_shifts_df["strand"]=="fwd"]
species_names_df=X_info_all_species_all_shifts_df_fwd.iloc[:,:2].copy()
species_names_df.reset_index(drop=True, inplace=True)

X_info_all_species_all_shifts_df_fwd.drop(columns=["species","species_trivial","strand"], inplace=True)
X_info_all_species_all_shifts_df_fwd.columns = [colname+"_fwd" for colname in X_info_all_species_all_shifts_df_fwd.columns]
X_info_all_species_all_shifts_df_fwd.reset_index(drop=True, inplace=True)


X_info_all_species_all_shifts_df_rv = X_info_all_species_all_shifts_df.loc[X_info_all_species_all_shifts_df["strand"]=="rv"]
X_info_all_species_all_shifts_df_rv.drop(columns=["species","species_trivial","strand"], inplace=True)
X_info_all_species_all_shifts_df_rv.columns = [colname+"_rv" for colname in X_info_all_species_all_shifts_df_rv.columns]
X_info_all_species_all_shifts_df_rv.reset_index(drop=True, inplace=True)

X_info_all_species_merged = pd.concat([species_names_df, X_info_all_species_all_shifts_df_fwd, X_info_all_species_all_shifts_df_rv], axis=1)

X_info_all_species_merged_long = X_info_all_species_merged.melt(
    id_vars="species_trivial", 
    value_vars=X_info_all_species_merged.columns[2:],
    var_name="Prediction",
    value_name="Value"
)

human = ['Human']
apes = ['Human', 'Chimp', 'Bonobo', 'Gorilla']
nonhuman_apes = ['Chimp', 'Bonobo', 'Gorilla']
old_world_monkeys = ['Orangutan', 'Gibbon', 'Rhesus', 'Crab-eating macaque', 'Pig-tailed macaque', 'Sooty_mangabey', 'Baboon', 'Green_monkey', 'Drill', 'Proboscis_monkey', 'Angolan_colobus', 'Golden_snub-nosed monkey', 'Black_snub-nosed monkey']
new_world_monkeys = ['Marmoset', 'Squirrel_monkey', 'White-faced_sapajou', "Ma's_night monkey"]
other_primates = ['Tarsier', 'Mouse_lemur', "Coquerel's_sifaka", 'Black_lemur', "Sclater's_lemur", 'Bushbaby']
other_mammals = ['Mouse', 'Dog', 'Armadillo']


#X_info_all_species_merged.loc[X_info_all_species_merged["species_trivial"].isin(apes)]["pred_pTE_accessibility_median"].median()
#median human
median_human = np.nanmedian(X_info_all_species_merged.loc[X_info_all_species_merged["species_trivial"].isin(human)].iloc[:,2:].to_numpy())
median_apes = np.nanmedian(X_info_all_species_merged.loc[X_info_all_species_merged["species_trivial"].isin(apes)].iloc[:,2:].to_numpy())
median_nonhuman_apes = np.nanmedian(X_info_all_species_merged.loc[X_info_all_species_merged["species_trivial"].isin(nonhuman_apes)].iloc[:,2:].to_numpy())
median_old_world_monkeys = np.nanmedian(X_info_all_species_merged.loc[X_info_all_species_merged["species_trivial"].isin(old_world_monkeys)].iloc[:,2:].to_numpy())
median_new_world_monkeys = np.nanmedian(X_info_all_species_merged.loc[X_info_all_species_merged["species_trivial"].isin(new_world_monkeys)].iloc[:,2:].to_numpy())

fc_human_apes = median_human/median_nonhuman_apes
fc_apes_old_world_monkeys = median_apes/median_old_world_monkeys
fc_apes_new_world_monkeys = median_apes/median_new_world_monkeys

# Create a boxplot
fig, ax = plt.subplots(figsize=(15, 7))

sns.boxplot(x="species_trivial", y="Value",whis=[25, 75], data=X_info_all_species_merged_long,fliersize =0, ax=ax)
sns.swarmplot(x="species_trivial", y="Value", data=X_info_all_species_merged_long,ax=ax, alpha = 0.5, c="tab:grey", s = 2)

ax.text(
    0.95, 0.95,  # Coordinates near the upper-right corner
    f"Foldchange Human/other Great apes: {fc_human_apes:.2f}\nFoldchange Great apes/Old World Monkeys: {fc_apes_old_world_monkeys:.2f}\nFoldchange Great apes/New World Monkeys: {fc_apes_new_world_monkeys:.2f}",
    fontsize=10,
    ha='right',  # Align text to the right
    va='top',    # Align text to the top
    transform=ax.transAxes  # Relative to the axes (0 to 1)
)

ax.tick_params(axis='x', labelrotation=90)
ax.set_ylabel("Predicted accessibility (ATAC log1p)")
ax.set_xlabel("Species")
coordinates_unshifted_window = "_".join([faname for faname in  ROI_name_fa_list if "shift_0" in faname][0].split("_")[-5:-2])
ax.set_title(f"Cross species variant prediction - {ROI_name} {coordinates_unshifted_window}")

plt.tight_layout()
plt.show()

#save foldchange df to csv, with ROI name
os.makedirs(os.path.join(res_dir,"predictions"), exist_ok=True)
fc_df = pd.DataFrame({"candidate_enhancer_name":[ROI_name],"fc_human_apes":[fc_human_apes], "fc_apes_old_world_monkeys":[fc_apes_old_world_monkeys], "fc_apes_new_world_monkeys":[fc_apes_new_world_monkeys]})
fc_df.to_csv(os.path.join(res_dir,f"predictions/{ROI_name}_{coordinates_unshifted_window}_species_foldchange.csv"), index=False)

#savefig
fig.savefig(os.path.join(res_dir,f"predictions/cross_species_predictions_{ROI_name}_{coordinates_unshifted_window}.pdf") ,bbox_inches='tight', dpi=600)