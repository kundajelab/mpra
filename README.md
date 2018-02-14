# MPRA

This project applies convolutional neural networks to predict output from massively parallel reporter assays (MPRAs), with the aim of systematically decoding regulatory sequence patterns.

Code written by Rajiv Movva, with advice from Peyton Greenside and Anshul Kundaje.

Most of the code involved with training and interpreting the models can be found at https://github.com/kundajelab/mpra/tree/master/deeplearn/scripts.

Some code written to process the Sharpr-MPRA and SuRE-seq datasets is available in the folders "sharpr" and "sureseq" respectively. These scripts process the data and produce structured input/output matrices that are used to train the deep learning models.
