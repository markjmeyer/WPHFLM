# WPHFLM
Code to implement Bayesian Wavelet-packet Historical Functional Linear Models as described in Meyer, Malloy, and Coull (2021).

sample_script.m contains a sample script with a single simulated setting for illustration.

The Simulations folder contains scripts to run simulations from the manuscript, sim_n50 and sim_n200 for the cumulative (c), delayed time-specific (d), lagged (l), and time-specific (t) effects. The peak effect is in sim_n20p. The comparison simulations, comp_ref_mr_, call R from MATLAB using flat-file communication. This requires the file refFDB.R to be located in the MATLAB directory. It also requires the user to install the refund and FDBoost packages in R as well as the R.matlab package.
