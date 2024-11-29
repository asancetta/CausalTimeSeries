# CausalTimeSeries
This code replicates the empirical results of the paper Causal Inference for High-Dimensional Time Series by Cordoni and Sancetta.
However, it is written so that it can be modfiied and applied to general data-sets. 
The file CausalTimeSeries_fx.R contains generic functions loaded by the two scripts script_MACRO_application.R and script_HFT_application.R. 

script_MACRO_application.R replicates the results of the macro application. It uses the datafile data.csv.

script_HFT_application.R replicates the high frequency results. It uses the datafile Xall_volume_time.csv. The additional files HFT_application_PC_result.RData and lambda_tau_cv_hft_empirical.RData contain intermediate results. They will be automatically generated if flag_first_run in the script script_HFT_application.R is set to true. 

