#!/bin/bash

source activate tf_modisco_light 
module load meme/5.1.1-foss-2018b-python-3.6.6

JASPAR_VERT_FILE=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/ext_data/JASPAR_vertebrate_non_redundant_motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
MODISCO_h5=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_open_regions/modisco/modisco_results_fold06.h5
OUTDIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/motif_discovery/pTE_open_regions/modisco/
$MYBSUB -m 30 -n modisco_fold06 -o ./log -T 01:00:00  "modisco report -i ${MODISCO_h5} -o ${OUTDIR} -m ${JASPAR_VERT_FILE}"