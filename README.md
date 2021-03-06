# lcms

This repository contains the data, MATLAB code, and Stan code used in the following paper "Toward the general mechanistic model of liquid chromatographic retention"

The data can also be found in: Kubik, Ł.; Jacyna, J.; Struck-Lewicka, W.; Markuszewski, M.; Wiczling, P. LC-TOF-MS Data Collected for 300 Small Molecules. XBridge-C18 Column. 2022. https://doi.org/10.17605/OSF.IO/ZQTJ7.

The files:
1) hplc-gra-redsum-qsrr-L.m is the main script.
2) boxplot_pwhisker.m, hplc_gra_sim.m, plot_data.m, plot_sim.m, plot_uncertainity_chromatogram.m are helper functions.
3) hplc-gra-redsum-qsrr-L.stan (main), hplc-gra-redsum-qsrr-L-priors.stan (for prior predictions), hplc-gra-redsum-qsrr-L-fixed.stan (with fixed population parameters) are Stan programs 
4) The full summary of the distributions of estimated parameters and derived quantities is in hplc-gra-redsum-qsrr-L.txt

The R code is available here: https://github.com/akamedulska/lc-ms