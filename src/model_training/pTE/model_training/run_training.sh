#!/bin/bash

#training with importance weights

SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/model_training/pTE/model_training/Train_model_with_weighing_TESTING.py
VARIABLE=ATAC_pTE_log1p

#OUTDIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/model_training/trained_models/
OUTDIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/model_training/trained_models_test/
DATADIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/data_selection/fold_fa_txt_te_va_tr/

FOLD=fold01
REP=999
modelname=pTE_${FOLD}_rep_${REP}
$SCRIPT -i $FOLD -v $VARIABLE -a DeepSTARR2 -d $DATADIR -o $OUTDIR -n $modelname