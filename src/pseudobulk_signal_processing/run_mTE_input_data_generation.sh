#!/bin/bash

SCRIPT=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/pseudobulk_signal_processing/Prepare_input_data_from_pseudobulk.sh
BW=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/mTE/mural_TE_cell_96h_rpm_norm.bw
HG38=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/ext_data/reference_genomes/hg38/hg38.chrom.sizes
OUTDIR=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/mTE

bash $SCRIPT -i $BW -o $OUTDIR -c $HG38

#bedtools v2.31.1
#bigWigAverageOverBed v2