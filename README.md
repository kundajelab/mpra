# MPRA

This project applies convolutional neural networks to predict output from massively parallel reporter assays (MPRAs), with the aim of systematically decoding regulatory sequence patterns and identifying disease-causing noncodign variants.

Most of the code involved with training and interpreting the models can be found at https://github.com/kundajelab/mpra/tree/master/deeplearn/scripts. The parameters of the final model described in the manuscript are here: https://github.com/kundajelab/mpra/blob/master/deeplearn/yamls/sharpr_znormed_jul23_config/hyperparameter_configs_list.yaml.

Some code written to process the Sharpr-MPRA and SuRE-seq datasets is available in the folders "sharpr" and "sureseq" respectively. These scripts process the data and produce structured input/output matrices that are used to train the deep learning models.

Feel free to contact Rajiv Movva (rmovva at mit dot edu) with any questions.
