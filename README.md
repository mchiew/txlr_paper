# Results and Figures 

## Accelerated calibrationless parallel transmit mapping using joint transmit and receive low-rank tensor completion
Aaron T. Hess, Iulius Dragonu, Mark Chiew

Within each folder (except for Figure 1), there are `fig_XX.m` files which generate the data/results for that figure or sub-figure, and `plt_XX.m` files which generate the MATLAB figure. Generally speaking, `fig_XX.m` will save output to a `res_XX.mat` mat-file, which is read by `plt_XX.m`

Requires the following data (which can be downloaded from <https://doi.org/10.5287/bodleian:YQpGNevaa>):
- `data/masks.mat` poisson-disc sampling masks
- `data/syn_data.mat` synthetic body dataset (8 Rx, 8 Tx)
- `data/body_data.mat` in-vivo body datasets with Rx noise scans (8 Rx, 8 Tx)
- `data/brain_data.mat` in-vivo brain dataset with Rx noise scan (32 Rx, 8 Tx)
