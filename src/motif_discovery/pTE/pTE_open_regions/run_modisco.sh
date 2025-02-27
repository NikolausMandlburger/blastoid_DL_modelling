#!/bin/bash

#Run using env with modisco-lite e.g. env <root_dir>/blastoid_DL_modelling/src/envs/modisco_lite_env.yml

SEQS=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_open_regions/attribution_profiles/fold06_sequences_onehot.npz
CONTRIBS=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_open_regions/attribution_profiles/fold06_sequences_contrib.npz
OUT_H5=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_open_regions/modisco/modisco_results_fold06.h5
$MYBSUB -m 60 -o ./log -n modisco_fold06 -T 14:20:00 "modisco motifs -s ${SEQS} -a ${CONTRIBS} -n 50000 -o ${OUT_H5}"