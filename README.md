# Confidence_Sets_GLM

Adapted from https://github.com/AlexBowring/Confidence_Sets_Manuscript.

Compute Cohen's d confidence sets in a GLM setting, at the group level.

The algorithm is adapted from Algorithm 2 in https://doi.org/10.1016/j.neuroimage.2020.117477.

## Inputs
- contrast4d (string):        path to a 4d nii file of individual contrast estimates (created for example with fslmerge). gz is not supported. 
- brainmask (string):         path to a binary brain mask. Dimensions must be identical to contrast4d.
- design (string):            path to the group-level design matrix in a csv format, with named columns. Should include an intercept: a column named "CONSTANT" consisting of 1s.
- contrastvector (vector):    contrast representing which covariate to test (e.g. [1,0,0]). The order should be the same as in the design matrix.
- thr (float, optional):      cohen's d threshold, default 0.5.
- outputfolder (string, optional):    path to a folder where the ouputs will be written.

