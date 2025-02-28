Description of files:

Files represent data for training, validating and testing pTE accessibility models.
Files have been produced with the following notebook: <root_dir>/blastoid_DL_modelling/src/model_training/pTE/data_selection/dataselection_pTE.Rmd

fold<n>_sequences_activity_<Train,TestVal>.txt:
Accessibility values and "class" labels for train, validation and test data for pTE models.

fold<n>_sequences_<Train,TestVal>.fa:
Accessibility values and "class" labels for train, validation and test data for pTE models.
Input sequences for train, validation and test data for pTE models.

NOTE: The files in this directory contain redundant information, the Train, test, Val files for any single fold contain the information for all of the dataset.
Files named according to different folds differ in terms of which sequences have been assigned to Train, val or test set. 