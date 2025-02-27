#!/bin/bash

#############################################################################
#Generate rpm normalised ATAC track of 96h blastoid polar trophectoderm cells
#############################################################################

frag2bw=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/src/pseudobulk_signal_processing/fragment2bw_cov.sh
FRAGMETSTSV=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/data/pTE/fragments.tsv.gz
OUT_BW=/groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/res/pseudobulks/pTE/polar_cell_96h_rpm_norm.bw
bash $frag2bw -i $FRAGMETSTSV -o $OUT_BW  -c /groups/stark/nikolaus.mandlburger/Projects/blastoid_project/blastoid_DL_modelling/ext_data/reference_genomes/hg38/hg38.chrom.sizes

#Run with environment that has:
#bedtools v2.31.1
#bedGraphToBigWig v 2.10
#e.g. <root_dir>/blastoid_DL_modelling/src/envs/bedtools_env.yml