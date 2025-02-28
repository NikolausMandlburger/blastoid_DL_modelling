## Computational modeling of sequence encoding of chromatin accessibility in blastoid cell types

This repository contains the code for the deep learning based analyses of the paper <b>"PLACEHOLDER"</b>,<br>
This includes preprocessing for and training of DeepSTARR type models for chromatin accessibility prediction from sequence as well as downstream application of<br> 
these models for de novo motif discovery as well as variant effect predictions of enhancer candidates across different primate species.<br>
<br>
<b>Structure:</b>
<br>
<br>
<u>src:</u><br>
Contains all scripts for performing data preprocessing (starting from ATAC pseudobulks), 
model training, motif discovery and cross species variant effect prediction.<br>
<br>
<u>res:</u><br>
Contains...<br>
- trained pTE accessibility models (keras models, weights in h5 format, structure in json format) 
    - /res/model_training/pTE/model_training/trained_models<br>
<br>
- Final results of motif discovery
    - /res/motif_discovery/pTE_specopen_vs_mTE/modisco
    - /res/motif_discovery/pTE_open_regions/modisco<br>
<br>
- Final results of cross species variant effect prediction
    - /res/primate_sequence_variant_prediction/GCM1/predictions
    - /res/primate_sequence_variant_prediction/GRHL1/predictions
    - /res/primate_sequence_variant_prediction/TCF7L2/predictions

<br><b>Raw data and intermediary data files are not contained in the repository for space reasons</b><br>however readme.txt files distributed throughout the different directories of the repository hold information on which intermediary files are expected, how they are generated and what they are.
