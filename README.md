# Reproducibility Repository for "Prediction-Powered Inference with Imputed Covariates and Nonuniform Sampling"

This Github repository gives the reproducibility code for the paper entitled "Prediction-Powered Inference with Imputed Covariates and Nonuniform Sampling" by Dan M. Kluger, Kerri Lu, Tijana Zrnic, Sherrie Wang, and Stephen Bates. It also provides a function titled "PTDBootModularized" that implements various variants of the Predict-Then-Debias Bootstrap simultaneously. Descriptions for the input argument to the function "PTDBootModularized.R" can be found at the top of that function's file.

Please note that there a couple of notational differences between the reproducibility code and the corresonding paper. Most notably, the code uses the term "calibration" or "calib" to refer to the small sample where all variables of interest are measured while the paper refers to this sample as the "complete" sample. In addition the code referes to the larger sample with missing variables as the "main" or "noncalibration" sample whereas the paper refers to it as the "incomplete" sample. Also note that the paper on arXiv used an unspecified seed, whereas the reproducibility code involved preset seed, so the results from running the reproducibility code differ slightly from those in the paper.

Please direct any questions or bug reports to dkluger@mit.edu. 
