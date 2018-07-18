# MPRA

This project applies convolutional neural networks to predict output from massively parallel reporter assays (MPRAs), with the aim of systematically decoding regulatory sequence patterns and identifying disease-causing noncodign variants.

Most of the code involved with training and interpreting the models can be found in Jupyter notebooks in `deeplearn/scripts/`. The final model was trained in Keras 1.2.2 and can be loaded with the `load_model` or `model_from_json` methods using the weights (`deeplearn/model_files/record_13_model_bgGhy_modelWeights.h5`) and/or architecture files (`deeplearn/model_files/record_13_model_bgGhy_modelJson.json`).

Dependencies:
* numpy, scipy, pandas, seaborn, etc.
* theano (0.9.0)
* keras (1.2.2)
* deeplift (0.5.2)

Here are descriptions of the relevant Jupyter notebooks in `deeplearn/scripts/`:

* `Sharpr Model Interpretation.ipynb`: Scatter plots for replicates (Figure 1C), prediction performances (Figure 2A), performances for specific chromatin states (Figure 2B), exploration of regulatory motifs and grammars learned by the model.
* `GBM Performance Benchmarking.ipynb`: Training and performance testing of gradient boosting tree models (Figure S2).
* `Sharpr DeepLIFT Scoring Validation`: CENTIPEDE TFBS validation (Figure 2C-D), DeepLIFT motif scores, comparison to Sharpr (Fig. 3).
* `Regulatory Grammar Discovery with Sharpr Models.ipynb`: Exploration of predictive TF motif PWMs by comparing to DeepLIFT score profiles (Fig. 4).
* `Variant Scoring.ipynb`: Evaluation of SNPpet ISM scores for variant prioritization (Figures 5 and 6, Supp. Figures).
* `ggplot2 visualizations.ipynb`: R notebook to generate several of the manuscript's figures.

Some code written to process the Sharpr-MPRA and SuRE-seq datasets is available in the folders "sharpr" and "sureseq" respectively. These scripts process the data and produce structured input/output matrices that are used to train the deep learning models.

Feel free to contact Rajiv Movva (rmovva at mit dot edu) with any questions.
