Description of files:

mural_TE_cell_96h_rpm_norm.bg, mural_TE_cell_96h_rpm_norm.bw:
rmp normalised ATAC track of 96h blastoid polar cells, produced by running <rootdir>/blastoid_DL_modelling/src/pseudobulk_signal_processing/run_mTE_pseudobulk_processing.sh

hg38.1001bp_50s.windows.bed, Regions_coverage_mural_TE_cell_96h_rpm_norm.bed:
Produced by running <root_dir>/blastoid_DL_modelling/src/pseudobulk_signal_processing/run_mTE_input_data_generation.sh

hg38.1001bp_50s.windows.bed:
1001 bp long windows of hg38 genome, tiled with stride 50

Regions_coverage_mural_TE_cell_96h_rpm_norm.bed:
The regions of hg38.1001bp_50s.windows.bed with the averaged (rpm normalised) ATAC signal of the pseudobulk of 96h mural trophecoderm cells over the central 200 bp parts of the windows.