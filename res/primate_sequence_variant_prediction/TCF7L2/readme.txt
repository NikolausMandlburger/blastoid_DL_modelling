Description of files:

Enhancer candidates of TCF7L2: Human versions and homologous sequences of different primate and mammal species
.maf files: multiple sequence alignments of loci as obtained from primate cons 30 conservation track from UCSC genome browser
.fa files: Continuous form of the DNA sequences from maf files, trimmed to 1001 bp lenght, for input to the models

predictions subfolder:
Results of the model predictions for the different loci

Explanation "shifts":
For each candidate peak, sequence from the the center part has been used (shift_0) as well as sequences of windows that are shifted n basepairs relative to peak center (shift_n)

All files in this folder and the predictions subfolder have been generated by running the script <root_dir>/blastoid_DL_modelling/src/primate_sequence_variant_prediction/run_pipeline_TCF7L2.sh