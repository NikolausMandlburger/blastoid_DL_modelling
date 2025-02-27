#!/bin/bash

SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/motif_discovery/pTE/run_DeepSHAP_DeepExplainer_for_modiscolite.py
model_path=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/model_training/trained_models/
modelname=pTE_fold06_rep_1
outdir=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_specopen_vs_mTE/attribution_profiles/

background=dinuc_shuffle
setname=Test
input_len=1001

FASTA=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/model_training/pTE/data_selection/pTE_specopen_vs_mTE.fa
$MYBSUB -P g -G "gpu:1" -c "g2|g3" -m 50 -o ${outdir}/log/ -n attr -T 03:30:00 "$SCRIPT -i $FASTA -p $model_path -o $outdir -m $modelname -s $setname -w $input_len -b $background"