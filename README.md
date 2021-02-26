Code for the paper:

Extended Federated Ensemble Regression using Classification.
To reproduce the results having copied to repository to a local directory:
1. Download the `experiments datasets.zip` folder from http://dx.doi.org/10.17632/mpvwnhv4vb.1
2. There are four folders in the `experiments datasets` folder downloaded in (1). To reproduce the gene expression experiments, copy the `gene_expression` dataset folder into the `gene expression` folder in this repository and rename it to `datasets`. It should be in the same folder level as the `code` folder. To run the experiments, go into the `code` folder and execute the `eferuc.R` script. It will automatically create an `output` folder in the same level as the `code` and `datasets` folders.
The same process should be followed for the other dataset folders, however, the `general` folder in this repository should be used instead of the `gene_expression` folder. 

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