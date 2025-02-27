#!/bin/bash

#Run using env with modisco-lite e.g. env <root_dir>/blastoid_DL_modelling/src/envs/modisco_lite_env.yml

SEQS=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_specopen_vs_mTE/attribution_profiles/pTE_specopen_vs_mTE_onehot.npz
CONTRIBS=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_specopen_vs_mTE/attribution_profiles/pTE_specopen_vs_mTE_contrib.npz
OUT_H5=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_specopen_vs_mTE/modisco/modisco_results_pTEspec.h5
$MYBSUB -m 60 -o ./log -n modisco_pTEspec -T 4:20:00 "modisco motifs -s ${SEQS} -a ${CONTRIBS} -n 50000 -o ${OUT_H5}"