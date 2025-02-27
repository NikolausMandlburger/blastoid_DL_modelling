#!/bin/bash

SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/motif_discovery/pTE/run_DeepSHAP_DeepExplainer_for_modiscolite.py
model_path=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/model_training/trained_models/
modelname=pTE_fold06_rep_1
outdir=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_open_regions/attribution_profiles
open_fasta_dir=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/data_selection/folds_open_fw_vs_rv/

background=dinuc_shuffle
setname=Test
input_len=1001

for FOLD in fold01 fold02 fold03 fold04 fold05 fold07 fold08 fold09 fold10;do
    FASTA_FW=${open_fasta_dir}/${FOLD}_sequences_fw${setname}.fa
    FASTA_RV=${open_fasta_dir}/${FOLD}_sequences_rv${setname}.fa
    $MYBSUB -P g -G "gpu:1" -c "g2|g3|g1" -m 50 -o ${outdir}/log/ -n attr_fw_${FOLD} -T 03:30:00 "$SCRIPT -i $FASTA_FW -p $model_path -o $outdir -m $modelname -s $setname -w $input_len -b $background"
    $MYBSUB -P g -G "gpu:1" -c "g2|g3|g1" -m 50 -o ${outdir}/log/ -n attr_rv_${FOLD} -T 03:30:00 "$SCRIPT -i $FASTA_RV -p $model_path -o $outdir -m $modelname -s $setname -w $input_len -b $background"
done
