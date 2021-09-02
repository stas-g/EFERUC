Code for the paper:
Extended Federated Ensemble Regression using Classification (EFERUC).

To reproduce the EFERUC results having copied the correct codebase in the EFERUC folder to a local directory:
1. Download the `experiments datasets.zip` folder from http://dx.doi.org/10.17632/mpvwnhv4vb.2
2. There are four folders in the `experiments datasets` folder downloaded in (1). To reproduce the gene expression experiments, copy the `gene_expression` dataset folder into the `gene expression` folder in this repository and rename it to `datasets`. It should be in the same folder level as the `code` folder. To run the experiments, go into the `code` folder and execute the `eferuc.R` script. It will automatically create an `output` folder in the same level as the `code` and `datasets` folders.
The same process should be followed for the other dataset folders, however, the `general` folder in this repository should be used instead of the `gene_expression` folder. Some code folders have `_nominal` attached to them. This indicates that the datasets with nominal attributes should be performed using the code in that folder.

The output folder has the following structure:

```
+-- output
|   +-- dataset_name_A
|   |   +-- weight_details
|   |   +-- best_bin_size
|   |   +-- performance
|   |   |   +-- regression
|   |   |   +-- classification
|   |   |   |   +-- RDS
|   +-- dataset_name_B
|   ...
|   +-- dataset_name_C
```

To reproduce the resampling results use `resampler_script.R` in the `Resampling` folder. Perform the following:

1. Download `datasets_and_splits.zip` from [https://data.mendeley.com/datasets/mpvwnhv4vb/2](https://data.mendeley.com/datasets/mpvwnhv4vb/2) and unpack in your working directory. These are the datasets and training/test splits used in the paper.
2. Download the `Code and Data` folder from [https://www.dcc.fc.up.pt/~ltorgo/ExpertSystems/] and unpack in your working directory. This is code accompanying the paper by [Torgo et al](https://onlinelibrary.wiley.com/doi/abs/10.1111/exsy.12081).
3. You will need to install the `uba` [package](https://rdrr.io/github/rpribeiro/uba/) in order to use the relevance function. Note that one of the dependencies of the package has been taken of CRAN so you will need to install it by hand by first downloading the freshest version [here](https://cran.r-project.org/src/contrib/Archive/DMwR/) and then running `install.packages( "Path/To/DMwR_0.4.1.tar.gz", repos=NULL, type="source" )`.
4. Make sure `smoter_helper.R`, `nominal_var_info.rds` and `relevance.rds` are in your working directory, or amend the script (lines 9, 20, 22) so that it knows where it is.
5. To run an experiment (for all the learners: ranger, xgboost, lasso and ridge) for data collection `d`, dataset `a`, resampling method Smoter (`s=TRUE`) or undersampler (`s=FALSE`), undersampling level `u` and, for SmoteR, oversampling level `o`, execute:
`Rscript --vanilla resampler_script.R d=d a=a s=s u=u o=o`

- dataset collections are: `Yeast`, `QSAR`, `OpenML`, `gene_expression` and `PaoBrancoImbalanced`, corresponding to the four sub-folders in the `datasets_and_splits/datasets` folder.
- dataset names is the first portion of the file names in `datasets_and_splits/datasets/` subfolders, e.g. for `yeast_x5.Fluorocytosine_X.csv`, it is `yeast_x5.Fluorocytosine`.
- to run SmoteR, choose `s=TRUE`, else, to run undersampler, choose `s=FALSE`.
- undersampling and oversampling levels can be freely chosen by the user.

Example use:
`Rscript --vanilla resampler_script.R d=Yeast a=yeast_x4NQO s=TRUE u=50 o=100`

5. Results will be saved in `output` in the format `d_a_l_s_u_o_mod.rds` (model fitted on a resampled dataset), `d_a_l_s_u_o_pred.rds` (prediction obtained on the test set) and `d_a_l_s_u_o_rsq.rds` (R-squared).

To reproduce analysis performed in the paper, user will need to run `resampler_script.R` for all dataset collections, resampling methods and combinations of under- and oversampling detailed in our paper.