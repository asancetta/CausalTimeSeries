####################
# Cordoni Francesco and Alessio Sancetta 19 May 2023
# script for implementation of the procedure of the paper 
# "Consistent Causal Inference for High Dimensional Time Series" Cordoni Francesco and Alessio Sancetta
#####################


#-------------------------------------------- Main Algorithm and Function
source("CausalTimeSeries_fx.R")
#######################################################################################
#----------------------------------------------------------------------------------#
#                                                                   
# Load DATA and setting parameters                         
#----------------------------------------------------------------------------------#

# As example, we use the dataset in Kanzig 2021
# Example. VAR(12) with  
#    Real oil price"              "World oil production"       
#   "World oil inventories"       "World industrial production"
#   "U.S. industrial production"  "U.S. CPI"
# DATA set from Kanzig 2021 AER "The macroeconomic effects of oil supply news: Evidence
# from OPEC announcements" Diego R. Känzig ∗ We consider his Baseline VAR model

#Changes to be made to use the algorithms in more general context are discussed below.

#----------working directory - data dir

dataMacro = "data.csv"#use full path if needed

#-----------------estimation settings
p_lag = 12# number of lags in latent VAR
# number of smoothing parameters to construct an initial grid of lambdas in Algorithm 2
nlambda = 20
# ratio between smallest over largest lambda parameter
lambda_min_ratio = 0.01 
# set of values for tau divided by lambda, e.g. 2 means tau=2*lambda, in the paper 
tau_rel_set = c(1.2,2)# a grid size of tau relative to lambda (tau = tau_rel_set*lambda) 
fold_cv     = 5# number of folds for CV to choose smoothing parameters in Algorithm 2
seed_num_cv = 1# seed for exact replication

# number of lambdas to consider below and above lambda_aic
#(lambda_aic is an estimated smoothing parameter around which we build a grid of lambdas for cross validation )
n_lambdas_cv_interval = 3 
flag_minimise_false_neg = TRUE# selects a smaller penalty to reduce false negatives in Algorithm 2


alpha_pc    = 0.05 # the alpha parameter in Algorithm 4

alpha_irf_ci      = 0.1# confidence interval for irfs of size 1-alpha_irf_ci, 10% as in Kanzig 

flag_cumsum       = TRUE#compute the cumulative irfs
#In the paper, in Section 6.1 we carry out the same analysis twice, 
#where the second time we add the instrument. 
#When adding the instrument use the penalty obtained in the analysis with no instrument 
#Set the blow to true if want to use the instrument. 
#Best run with FALSE, obtain penalty by CV and hard code that one in lines 952 953 
#when this is set to TRUE
flag_use_instrument = FALSE

#M-B lasso graph calculation via Lasso or CLIME
flag_MB_CLIME = "MB"# Select CLIME for CLIME and MB for LASSO
flag_perform_CV = FALSE # set TRUE if you want to perform CV 
#(when we include the instrument we use the same penalty obtained by CV whitout the penalty, so set this to FALSE)

#------ load the data
data = as.matrix(read.csv(file = dataMacro,header = TRUE))
#remove first column because it is a date in test data and possibly remove the instrument
if( flag_use_instrument==FALSE){
  cols2remove = c(1,2)
  
}else{
  cols2remove = c(1)
  
}


X           = data[,-cols2remove]
col_names_X =as.character(colnames(X))

#--------Algorithm 1
K = ncol(X) 
n = nrow(X)

W         = get_W(X,p_lag)
K_W       = ncol(W)
Sigma_hat = huge.npn(W, npn.func="skeptic",verbose = FALSE)

#--------Algorithm 2 and 3 using cross-validation

if (flag_MB_CLIME == "MB"){
	mb_res  = huge(Sigma_hat,nlambda = nlambda,lambda.min.ratio =lambda_min_ratio,method = "mb",verbose = FALSE)
	# penalty chosen by the algorithm when we use the baseline model (i.e. no instrument)
	#this is for like for like comparison
	lambda_hard_coded =0.03779779810576886639 
	# To use CLIME replace the above function with 
	# mb_res  = sugm(Sigma_hat, nlambda = nlambda,lambda.min.ratio =lambda_min_ratio,method = "clime",verbose = FALSE)
	# it is slow
}else{

	# it is a little bit slow compared to LASSO
	# for CLIME we can use the sugm function which has the same output of the huge, no need to change any other functions
	mb_res  = sugm(Sigma_hat, nlambda = nlambda,lambda.min.ratio =lambda_min_ratio,method = "clime",verbose = FALSE)
	lambda_hard_coded = 0.0423346642385447058387626384501345455647
}




lambdas = mb_res$lambda# the set of lambdas generated within huge for M-B procedure
res_aic = get_mb_aic(Sigma_hat,mb_res,n, c(0))# thresh=0 for initial value of lambda and a set of lambdas
#the index in lambdas that minimises AIC
ind_lambda_aic = which.min(res_aic)
#a subset of lambdas around the AIC minimiser to be used as grid in cross validation

lambdas4cv     = lambdas[max(ind_lambda_aic-n_lambdas_cv_interval,1)
                         :min(ind_lambda_aic+n_lambdas_cv_interval,length(lambdas))]


#the set of all allowed combinations of lambda and tau for cross validation
lambda_tau_set = get_lambda_tau_set(lambdas4cv, tau_rel_set)

# lambda and tau derived from cross validation

if (flag_perform_CV){
  lambda_tau_cv = get_lambda_tau_cv(seed_num_cv, fold_cv, W, lambda_tau_set)
  lambda_cv    = lambda_tau_cv$lambda
  tau_cv       = lambda_tau_cv$tau
  obj_cv       = lambda_tau_cv$obj
  #this plot shows as many curves as length(tau_rel_set) one on the right of the other
  plot(obj_cv)
  
  # lambda and tau, possibly adjusted to reduce false positives
  # Choose smaller penalty to avoid false negatives 
  # if not adjusted, it just returns the input plus the index in lambdas and tau_rel_set above
  lambda_tau_res = get_ind4lambda_tau_sets(lambdas, tau_rel_set, lambda_cv,tau_cv,
                                           flag_minimise_false_neg = flag_minimise_false_neg)
  lambda_cv_adj  = lambda_tau_res$lambda
  tau_cv_adj     = lambda_tau_res$tau  
}else{
  # To replicate Figure 1 panel b of the paper 
  # select a penalty equal to the penalty selected by CV when using no instrument  
  # i.e. Baseline Model of Kanzig removing the OilShock variable.
  lambda_cv_adj  = lambda_hard_coded
  tau_cv_adj     = lambda_cv_adj*1.2
  
}



#reuse the results from mb_est rather than recomputing
Omega_hat      = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
Theta_hat      = get_Theta(Sigma_hat, Omega_hat)
Sigma_eps_hat  = get_Sigma_eps(Theta_hat,K)



#----------------- Algorithm 5

#let PC skip edges that are zero in Theta11, set to NULL otherwise
zero_rest      = Theta_hat[1:K,1:K]==0
pc_res         = get_pc_alg(Sigma_eps_hat,n, alpha=alpha_pc, 
                            zero_restr = zero_rest)

#---------------------- Algorithm 6
#extract the adjacency matrix from pc_res  
adj_mat     = pc_res$adj_mat


colnames(adj_mat) =  col_names_X
row.names(adj_mat) = col_names_X

g = graph_from_adjacency_matrix(adj_mat)

is.DAG(g)
plot(g, edge.arrow.size=0.5)


if (!is.DAG(g)){
  print("PC has returned a CPDAG instead of a DAG. This is could happen for two reasons:")
  print("1- The penalty is too small, so please select greater lambda/tau parameters.")
  print("You can also try different values of the penalties by investigating different penalty grids")
  print(" in the cross-validation step, i.e., by varying 'lambda_min_ratio' and 'nlambda'.")
  print("2- The assumptions of the PC could be violated, e.g., that the underlying system is not a DAG.")
  print("In this case you can continue to analyse the part of the graph which is identified, (partial identification).")
  print("For instance you can consider a shock, related to identified graph part, and consider the relative columns in the SVAR model")
  print("to study its effect on the other variable,")
  print("see the discussion in the Manuscript. This eventuality is not yet implemented in the code.")
  print("---------")
  print("Below you can find those edges which are not identified (path with 2 nodes)")
  print("and those paths (with more than 2 nodes) that violate the recursiveness assumption (acyclicity).")
  print("---------")
  print("A solution is to remove from the graph these paths (or at least one edges of the cycle) and check if the resulted.")
  print("graph is a DAG and continue the analysis.")
  
  find_directed_cycles = function(g) {
    Cycles = NULL
    for(v1 in V(g)) {
      if(degree(g, v1, mode="in") == 0) { next }
      GoodNeighbors = neighbors(g, v1, mode="out")
      GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
      for(v2 in GoodNeighbors) {
        TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
        TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
        Cycles  = c(Cycles, TempCyc)
      }
    }
    Cycles
  }
  
  
  cycles = find_directed_cycles(g)
  
  
  if (length(cycles)>=1){
    for (i in cycles){
      print(col_names_X[i])
    }
  }
}else{
  
  # 
  ####################################################################################
  #     THE FOLLOWING COULD BE IMPROVED
  #----------------------------------------------------------------------------------#
  #                                                                   
  # GET IRFS                              
  #----------------------------------------------------------------------------------#
  #This uses cumulative IRFs. Modify accordingly for general use
  max_horizon = 50#if fast PC use 50
  # number of sampling to compute the IRFs
  m = 500#if fast PC use 500
  delta_shock = 1#1/0.29541030
  # SELECT THE VARIABLE TO SHOCK
  k_shock_variable = c(1)  # we consider a shock on the REAL OIL PRICE variable 
  # Select the variable to be shocked
  k_shocked_variable = c(1:K)
  
  
  # YOU CAN CHOOSE THE GRAPH YOU WANT BY SELECTING THE RELATED Adjacency matrix
  adj_mat     = pc_res$adj_mat
  #find the SVAR parameters
  svar_params = get_SVAR_params(Sigma_eps_hat, adj_mat)
  PI          = svar_params$PI
  D           = svar_params$D
  H           = svar_params$H
  # (eye(K)-D) == solve(H)
  
  # ?? start
  # permute the H to the original order provided in input (i.e., the nominal order)
  H_permuted = t(PI)  %*% H %*% (PI)  
  colnames(H_permuted) = col_names_X
  row.names(H_permuted) = col_names_X
  # following Lemma 2 of the paper A^s H_permuted will provide IRFs for the latente VAR  
  # end
  # ------------------------------------------------------------ get the IRFs 
  # 
  f_k_inverse = get_f_k_inverse(X)
  
  # bootstrap for confidence interval 
  svar_params_comp = get_svar_params_companion(Sigma_eps_hat, Theta_hat, H, PI,  K)
  A_comp = svar_params_comp$A
  
  # on average it takes 262.158 sec , but you can choose different mboot and m parameters
  
  # PI_comp = svar_params_comp$PI
  # H_comp = svar_params_comp$H
  
  # IRF_latent_svar = get_IRF_Z(H_comp,PI_comp,A_comp, max_horizon, Lag=1)
  # then try to recovere these IRFs using bootstrap``
  
  #===================================
  # IRFs_X_point_est <- get_IRF_X(m = m, adj_mat, Sigma_eps_hat,
  #                               Theta_hat, f_k_inverse, max_horizon = max_horizon,
  #                               k_shock_variable = k_shock_variable, 
  #                               delta_shock = delta_shock, seed_num=0,
  #                               flag_shock_Z=FALSE, flag_cum = TRUE)
  # start
  IRFs_X_point_est_all = NULL
  for (ii_kk in 1:length(k_shock_variable)){
    
    i_k_shock_variable = k_shock_variable[ii_kk]
    IRFs_X_point_est_all[[ii_kk]] <-get_IRF_X(m = m, adj_mat, Sigma_eps_hat,
                                              Theta_hat, f_k_inverse, max_horizon = max_horizon,
                                              k_shock_variable = i_k_shock_variable,
                                              delta_shock = delta_shock, seed_num=0,
                                              flag_shock_Z=FALSE, flag_cum = TRUE)
  }
  # restrict the IRF to only the k_shocked_variable to be consistent with IRFs_X_boot
  
  IRFs_X_point_est = matrix(NA, 
                            nrow=max_horizon+1,
                            ncol=length(k_shock_variable)*length(k_shocked_variable))
  n_var_schocked = length(k_shocked_variable)
  for (ii_kk in 1:length(k_shock_variable)){
    selected_columns <-  ((  (ii_kk-1)*n_var_schocked+1):( ii_kk*n_var_schocked))
    IRFs_X_point_est[,c(selected_columns)] = IRFs_X_point_est_all[[ii_kk]][, k_shocked_variable]
  }
  # end
  #===================================
  
  tictoc::tic()
  # IRFs_X_boot = get_IRFs_bootstrap(mboot = 200, n, m=m, Sigma_hat, Sigma_eps_hat, A_comp,
  #                                f_k_inverse, adj_mat,
  #                                p_lag, K, max_horizon, k_shock_variable = 1,
  #                                lambda_cv_adj, tau_cv_adj, seed = 0)
  
  flag_shock_Z = FALSE # compute irf on Z or on X
  IRFs_X_boot = get_IRFs_bootstrap_selected_shocks(mboot = 200, n, m=m, Sigma_hat, Sigma_eps_hat, A_comp,
                                                   f_k_inverse, adj_mat,
                                                   p_lag, K, max_horizon,
                                                   k_shock_variable = k_shock_variable,
                                                   lambda_cv_adj, tau_cv_adj, seed_num = 0,
                                                   k_shocked_variable = k_shocked_variable,
                                                   flag_block = FALSE,'None','None',
                                                   flag_shock_Z = flag_shock_Z,
                                                   flag_cum = flag_cumsum,
                                                   delta_shock=delta_shock)
  tictoc::toc()
  
  
  
  #  IRFs_X_boot_CI = get_CI_boot(IRFs_X_boot, alpha = alpha_irf_ci)
  
  IRFs_X_boot_CI = get_CI_bootCentered(IRFs_X_boot,  alpha_irf_ci, IRFs_X_point_est)
  
  
  IRF_lwr_raw = t(IRFs_X_boot_CI$IRFs_X_boot_lwr)
  IRF_upr_raw = t(IRFs_X_boot_CI$IRFs_X_boot_upr)
  IRF_avg_raw = t(IRFs_X_boot_CI$IRFs_X_boot_mean)
  IRF_avg_raw = t(IRFs_X_point_est) # select point estimate
  
  
  #scale
  IRF_lwr = IRF_lwr_raw*(10/IRF_avg_raw[1,1])
  IRF_upr = IRF_upr_raw*(10/IRF_avg_raw[1,1])
  IRF_avg = IRF_avg_raw*(10/IRF_avg_raw[1,1])
  
  #####
  #IRF_lwr1 = t(IRF_lwr_raw*(10/IRF_avg_raw[1,1]))
  #IRF_upr1 = t(IRF_upr_raw*(10/IRF_avg_raw[1,1]))
  #IRF_avg1 = t(IRF_avg_raw*(10/IRF_avg_raw[1,1]))
  
  
  par(mar = c(3, 3, 2, 2))
  kk = 0
  
  
  
  for (ikshock in 1:length(k_shock_variable)){
    for (ikschoked in 1:length(k_shocked_variable)){
      kk = kk+1#nonsense: should use indexes in loop 
      
      filename_plot = sprintf("%s -> %s",col_names_X[k_shock_variable[ikshock]],
                              col_names_X[k_shocked_variable[ikschoked]])
      if (flag_shock_Z==TRUE){
        
        filename_plot = paste(filename_plot,'_Z',sep="")
        
      }
      #if using kk, what is the point of using k_shocked_variable[ikschoked]?
      png(paste("irf_macro/",filename_plot,sep=""), width = 800, height = 600)
      plot(IRF_avg[kk,],type="l", cex.axis=2,cex.main=2,
           # NOME MAIN con ->
           main =sprintf("%s -> %s",col_names_X[k_shock_variable[ikshock]],
                         col_names_X[k_shocked_variable[ikschoked]]),xlab="",ylab="",
           col='black',lty =1, lwd=3,
           ylim = c(min(IRF_lwr[kk,]), max(IRF_upr[kk,])))
      lines(IRF_lwr[kk,], type="l",col='red', lty=3)
      lines(IRF_upr[kk,], type="l",col='red', lty=3)
      #lines(irf_ols[col_names_X[k_shocked_variable[ikschoked]],], col=4)
      abline(h=0)
      dev.off()
      
    }
  }
  
  
  
}






# ####test 
# 
# IRFs_X_boot0 = get_IRFs_bootstrap(mboot = 200, n, m, Sigma_hat, Sigma_eps_hat, A_comp,
#                                f_k_inverse, adj_mat,
#                                p_lag, K, max_horizon, k_shock_variable,
#                                lambda_cv_adj, tau_cv_adj, seed_num = 0, delta_shock=1,
#                                flag_shock_Z=FALSE, flag_cum = TRUE)
#   
# 
# IRFs_X_boot = get_IRFs_bootstrap_selected_shocks(mboot = 200, n, m=m, Sigma_hat, Sigma_eps_hat, A_comp,
#                                                    f_k_inverse, adj_mat,
#                                                    p_lag, K, max_horizon,
#                                                    k_shock_variable = k_shock_variable,
#                                                    lambda_cv_adj, tau_cv_adj, seed_num = 0,
#                                                    k_shocked_variable = k_shocked_variable,
#                                                    flag_block = FALSE,'None','None',
#                                                    flag_shock_Z = flag_shock_Z,
#                                                    flag_cum = flag_cumsum,
#                                                    delta_shock=delta_shock)
#   
#   
#   


##################################################################
#############################################replicate kanzig code
####################################################################

get_IRF_from_A_single<-function(A,K, S1,max_horizon){
  #A is in companion for if not a VAR(1)  
  p_lagK = dim(A)[1] 
  irf_X   = matrix(0,max_horizon+1,K )
  Ai        = diag(p_lagK)
  
  for (i in 1:(max_horizon+1)){
    irf_X[i,] = Ai[1:K,1:K]%*%S1
    Ai = A %*% Ai
    
  }
  return(irf_X)
}

get_A_companion <-function(A){
  K = dim(A)[1]
  p_lagK = dim(A)[2]
  p_lag  = p_lagK/K
  if (p_lag!=1){
    
    A_companion = matrix(0,K*p_lag,K*p_lag)
    A_companion[1:K,] = A 
    n_blocks <- p_lag-1
    A_companion[(K+1):(K*p_lag),1:(K*(p_lag-1))] = diag(K*(p_lag-1))
    A = A_companion}  
  return(A)
  
}


############irfs from kanzig code converted from matlab to R

get_irf_kanzig= function (X,p_lag,proxy_raw,max_horizon=50){
  
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  Xlag  = cbind(rep(1,dim(W)[1]),W[,(K+1):K_W])
  
  var_A  = solve(t(Xlag)%*%Xlag)%*%(t(Xlag)%*%W[,1:K])
  resid  =  W[,1:K] - Xlag%*%var_A 
  proxy  = proxy_raw[(p_lag+1):length(proxy_raw)]
  
  s11 = cov(cbind(resid[,1],proxy))
  s11 = s11[1,2]
  s1  = cov(cbind(proxy,resid))
  s1  = s1[2:dim(s1)[1],1]
  s1  = s1/s11
  
  var_A_companion = get_A_companion(t(var_A[2:dim(var_A)[1],]))
  
  
  IRF_X01 = get_IRF_from_A_single(var_A_companion,K,s1,max_horizon)
  
  IRF_X01 = apply(IRF_X01, 2, cumsum)
  return(IRF_X01)
}

proxy = data[,2]
IRF_kanzig =  get_irf_kanzig(X,p_lag,proxy)
colnames(IRF_kanzig) = col_names_X



##############irfs from adj matrix and var residuals 



get_irf_var_adj_mat = function(X,p_lag,adj_mat,max_horizon=50){
  
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  Xlag  = cbind(rep(1,dim(W)[1]),W[,(K+1):K_W])
  
  var_A  = solve(t(Xlag)%*%Xlag)%*%(t(Xlag)%*%W[,1:K])
  var_A_companion = get_A_companion(t(var_A[2:dim(var_A)[1],]))
  
  
  resid  =  W[,1:K] - Xlag%*%var_A 
  Sigma_eps_hat0 = cov(resid)#t(resid)%*%resid/dim(resid)[1]
  svar_parms0    = get_SVAR_params(Sigma_eps_hat0, adj_mat)
  s1_adj         = t(svar_parms0$PI)%*%svar_parms0$H%*%svar_parms0$PI
  s1_adj         = s1_adj[,1]
  IRF_X01_adj    = get_IRF_from_A_single(var_A_companion,K,s1_adj,max_horizon)
  
  IRF_X01_adj = apply(IRF_X01_adj, 2, cumsum)
  return(IRF_X01_adj)
}

IRF_var_adj_mat = get_irf_var_adj_mat(X,p_lag,adj_mat)
colnames(IRF_var_adj_mat) = col_names_X


###################using our functions no transform#############################

get_irf_no_transform = function(X,p_lag,adj_mat,max_horizon=50){
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  Sigma_hat0     = cov(W)
  Theta_hat0     = solve(Sigma_hat0)
  
  var_params     = get_VAR_params(Theta_hat0,K)
  Sigma_eps_hat0 = var_params$Sigma_eps
  A_hat0         = var_params$A
  
  var_A_companion = get_A_companion(A_hat0)
  #find the SVAR parameters
  svar_params0 = get_SVAR_params(Sigma_eps_hat0, adj_mat)
  PI0          = svar_params0$PI
  D0           = svar_params0$D
  H0           = svar_params0$H
  
  s1_adj         = t(PI0)%*%H0%*%PI0
  s1_adj         = s1_adj[,1]
  IRF_X01_adj    = get_IRF_from_A_single(var_A_companion,K,s1_adj,max_horizon)
  
  IRF_X01_adj = apply(IRF_X01_adj, 2, cumsum)
  return(IRF_X01_adj)
}


IRF_no_transform =get_irf_no_transform(X,p_lag,adj_mat) 
colnames(IRF_no_transform) = col_names_X


###################using our functions for latent#############################

get_irf_latent = function(X,p_lag,adj_mat,max_horizon=50){
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  Sigma_hat0l = huge.npn(W, npn.func="skeptic",verbose = FALSE)
  Theta_hat0l     = solve(Sigma_hat0l)
  
  var_params     = get_VAR_params(Theta_hat0l,K)
  Sigma_eps_hat0l = var_params$Sigma_eps
  A_hat0l         = var_params$A
  
  var_A_companion_l = get_A_companion(A_hat0l)
  #find the SVAR parameters
  svar_params0 = get_SVAR_params(Sigma_eps_hat0l, adj_mat)
  PI0l          = svar_params0$PI
  D0l           = svar_params0$D
  H0l           = svar_params0$H
  
  s1_adjl         = t(PI0l)%*%H0l%*%PI0l
  s1_adjl         = s1_adjl[,1]
  IRF_X01_adj    = get_IRF_from_A_single(var_A_companion_l,K,s1_adjl,max_horizon)
  
  IRF_X01_adj = apply(IRF_X01_adj, 2, cumsum)
  return(IRF_X01_adj)
}


IRF_latent =get_irf_latent(X,p_lag,adj_mat) 
colnames(IRF_latent) = col_names_X

###################using our functions for latent, sparse#############################

get_irf_latent_sparse = function(X,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj,max_horizon=50){
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  Sigma_hat_l = huge.npn(W, npn.func="skeptic",verbose = FALSE)
  mb_res      = huge(Sigma_hat_l , lambda =lambda_cv_adj ,method = "mb",verbose = FALSE)
  Omega_hat_l = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
  Theta_hat_l = get_Theta(Sigma_hat_l, Omega_hat_l)
  Sigma_eps_hat_l  = get_Sigma_eps(Theta_hat_l,K)
  
  var_params     = get_VAR_params(Theta_hat_l,K)
  Sigma_eps_hat_l = var_params$Sigma_eps
  A_hat_l         = var_params$A
  
  var_A_companion_l = get_A_companion(A_hat_l)
  #find the SVAR parameters
  svar_params = get_SVAR_params(Sigma_eps_hat_l, adj_mat)
  PI_l          = svar_params$PI
  D_l           = svar_params$D
  H_l           = svar_params$H
  
  s1_adj_l         = t(PI_l)%*%H_l%*%PI_l
  s1_adj_l         = s1_adj_l[,1]
  IRF_X01_adj    = get_IRF_from_A_single(var_A_companion_l,K,s1_adj_l,max_horizon)
  
  IRF_X01_adj = apply(IRF_X01_adj, 2, cumsum)
  return(IRF_X01_adj)
}


IRF_latent_sparse =get_irf_latent_sparse(X,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj)
colnames(IRF_latent_sparse) = col_names_X



get_irf_latent_sparse2 = function(X,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj,max_horizon=50){
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  Sigma_hat  = huge.npn(W, npn.func="skeptic",verbose = FALSE)
  mb_res      = huge(Sigma_hat , lambda =lambda_cv_adj ,method = "mb",verbose = FALSE)
  Omega_hat   = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
  Theta_hat   = get_Theta(Sigma_hat, Omega_hat)
  Sigma_eps_hat   = get_Sigma_eps(Theta_hat,K)
  
  # ------------------------------------------------------------ get the IRFs
  svar_params = get_SVAR_params(Sigma_eps_hat, adj_mat)
  PI          = svar_params$PI
  D           = svar_params$D
  H           = svar_params$H
  
  svar_params_comp = get_svar_params_companion(Sigma_eps_hat, Theta_hat, H, PI,  K)
  # extract the parameters
  PI_comp = svar_params_comp$PI
  A_comp = svar_params_comp$A
  H_comp = svar_params_comp$H
  # IRFs on Z
  UPsi_comp = get_IRF_Z(H_comp,PI_comp, A_comp, max_horizon, Lag=1)
  
  
  e_shock = matrix(0, nrow = K, ncol = 1)
  e_shock[1] = delta_shock#/H[k_shock_variable,k_shock_variable] 
  IRF_X = matrix(0,nrow=max_horizon+1,ncol=K)
  for (s in 1:(max_horizon+1)){
    IRF_X[s,]   = UPsi_comp[[s]][1:K,1:K]%*%e_shock
  }
  
  if (flag_cum == TRUE){
    IRF_X = apply(IRF_X, 2, cumsum)
  }
  return(IRF_X)
}


IRF_latent_sparse2 =get_irf_latent_sparse(X,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj)
colnames(IRF_latent_sparse2) = col_names_X

###################using our functions sparse#############################

get_irf_sparse = function(X,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj,max_horizon=50, seed_num =0){
  W         = get_W(X,p_lag)
  K         =dim(X)[2]   
  K_W       = dim(W)[2]
  
  f_k_inverse = get_f_k_inverse(X)
  
  
  Sigma_hat = huge.npn(W, npn.func="skeptic",verbose = FALSE)
  mb_res      = huge(Sigma_hat , lambda =lambda_cv_adj ,method = "mb",verbose = FALSE)
  Omega_hat = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
  Theta_hat = get_Theta(Sigma_hat, Omega_hat)
  Sigma_eps_hat  = get_Sigma_eps(Theta_hat,K)
  
  IRF_X = get_IRF_X(m = 500, adj_mat, Sigma_eps_hat,
                    Theta_hat, f_k_inverse, max_horizon = 50,
                    k_shock_variable = 1, delta_shock = 1, seed_num=seed_num, 
                    flag_shock_Z=FALSE, flag_cum = TRUE)
  
  scale = 1/IRF_X[1,1]
  
  return(IRF_X*scale)
}



IRF_sparse = get_irf_sparse(X,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj)  
colnames(IRF_sparse) = col_names_X



#####plots 4 checks
kk = 6
print(col_names_X[kk])
ylim_min = min(rbind(IRF_var_adj_mat[,kk],IRF_kanzig[,kk],IRF_latent[,kk],IRF_latent_sparse[,kk],IRF_sparse[,kk]))
ylim_max = max(rbind(IRF_var_adj_mat[,kk],IRF_kanzig[,kk],IRF_latent[,kk],IRF_latent_sparse[,kk],IRF_sparse[,kk]))

ylim = c(ylim_min,ylim_max)
plot(IRF_kanzig[,kk],type="l", cex.axis=2,cex.main=2,xlab="",ylab="", ylim =ylim, col = 2)
#lines(IRF_var_adj_mat[,kk], type="l",col=2, lty=3)
#lines(IRF_latent[,kk], type="l",col=2, lty=1)
#lines(IRF_latent_sparse[,kk], type="l",col=2, lty=4)
#lines(IRF_latent_sparse2[,kk], type="l",col=2, lty=4)
lines(IRF_sparse[,kk], type="l",col=2, lty=5)


######plot results



shock_size  = 10#%
IRF_kanzig2plot = IRF_kanzig*shock_size
IRF_sparse2plot = IRF_sparse*shock_size 

par(mar = c(3, 3, 2, 2))

for (ikshock in 1:length(k_shock_variable)){
  for (ikshocked in 1:length(k_shocked_variable)){
    indshock   = k_shock_variable[ikshock]
    indshocked = k_shocked_variable[ikshocked]
    filename_plot = sprintf("%s -> %s",col_names_X[indshock],
                            col_names_X[indshocked])
    
    filename_plot = paste('irf_macro/kanzig_vs_gcvar_',filename_plot,sep="")
    # if (flag_shock_Z==TRUE){
    #   
    #   filename_plot = paste(filename_plot,'_Z',sep="")
    #   
    # }
    png(filename_plot, width = 800, height = 600)
    ylim_min = min(rbind(IRF_kanzig2plot[,indshocked],IRF_sparse2plot[,indshocked]))
    ylim_max = max(rbind(IRF_kanzig2plot[,indshocked],IRF_sparse2plot[,indshocked]))
    
    
    plot(IRF_kanzig2plot[,indshocked],type="l", cex.axis=2,cex.main=2,
         # NOME MAIN con ->
         main =sprintf("%s -> %s",col_names_X[indshock],
                       col_names_X[indshocked]),xlab="",ylab="", col = 2, lty=3,
         ylim = c(ylim_min, ylim_max))
    lines(IRF_sparse2plot[,indshocked], col=4)
    abline(h=0)
    dev.off()
    
  }
}

########################PLOT AGAIN, BUT WITH BOOTSTRAP INTERVALS###############
############============bootstrap 4 confidence intervals#######################

get_IRFs_bootstrap_simple_cum <- function(mboot = 200, n, m, Sigma, Sigma_eps, A_comp,
                                          f_k_inverse, adj_mat,
                                          p_lag, K, max_horizon, k_shock_variable,
                                          lambda_cv_adj, tau_cv_adj, seed_num = 0, delta_shock=1){
  # this function will compute IRFs for a selected k_shock_variable
  # on all the other variables
  #set.seed(seed)
  
  # parallalized 
  num_cores <- detectCores()
  registerDoParallel(cores = num_cores)
  IRFs_X_boot <- foreach(kboot = 1:mboot, .combine = "cbind") %dopar% {
    X_samples = get_X_samples_bootstrap (n, K, p_lag, Sigma, Sigma_eps, A_comp, f_k_inverse,seed_num=(seed_num+2)*kboot )
    
    IRFs_X = get_irf_sparse(X_samples,p_lag,adj_mat,lambda_cv_adj,tau_cv_adj,
                            max_horizon=50, seed_num =(seed_num+1)*kboot)  
    #colnames(IRF_sparse) = col_names_X
    
    
    # Return IRFs_X for current iteration
    IRFs_X
  }
  # Convert the result to the desired array format
  IRFs_X_boot <- array(IRFs_X_boot, dim = c(max_horizon + 1, K, mboot))
  return(IRFs_X_boot)
  
}


tictoc::tic()

flag_shock_Z = FALSE # compute irf on Z or on X
f_k_inverse = get_f_k_inverse(X)

IRFs_X_boot = get_IRFs_bootstrap_simple_cum(mboot = 200, n, m, Sigma_hat, Sigma_eps_hat, A_comp,
                                            f_k_inverse, adj_mat,
                                            p_lag, K, max_horizon, k_shock_variable,
                                            lambda_cv_adj, tau_cv_adj, seed = 0, delta_shock=1)

tictoc::toc()



IRFs_X_boot_CI = get_CI_bootCentered(IRFs_X_boot, alpha = alpha_irf_ci,IRF_sparse )

IRF_lwr_raw = (IRFs_X_boot_CI$IRFs_X_boot_lwr)
IRF_upr_raw = (IRFs_X_boot_CI$IRFs_X_boot_upr)
#IRF_avg_raw = (IRFs_X_boot_CI$IRFs_X_boot_mean)
IRF_avg_raw = IRF_sparse


#scale
IRF_sparse_lwr = IRF_lwr_raw#*(1/IRF_avg_raw[1,1])
IRF_sparse_upr = IRF_upr_raw#*(1/IRF_avg_raw[1,1])
IRF_sparse_avg = IRF_avg_raw#*(1/IRF_avg_raw[1,1])


#################PLOTS: THIS SHOULD JUST BE A FUNCTION- COULD BE IMPROVED

shock_size  = 10#%
IRF_kanzig2plot = IRF_kanzig*shock_size
IRF_sparse2plot = IRF_sparse*shock_size 

IRF_sparse_lwr2plot = IRF_sparse_lwr*shock_size
IRF_sparse_upr2plot = IRF_sparse_upr*shock_size
IRF_sparse_avg2plot = IRF_sparse_avg*shock_size

par(mar = c(3, 3, 2, 2))

for (ikshock in 1:length(k_shock_variable)){
  for (ikshocked in 1:length(k_shocked_variable)){
    indshock   = k_shock_variable[ikshock]
    indshocked = k_shocked_variable[ikshocked]
    filename_plot = sprintf("%s -> %s",col_names_X[indshock],
                            col_names_X[indshocked])
    
    filename_plot = paste('irf_macro/kanzig_vs_gcvar_boot_',filename_plot,sep="")
    # if (flag_shock_Z==TRUE){
    #   
    #   filename_plot = paste(filename_plot,'_Z',sep="")
    #   
    # }
    png(filename_plot, width = 800, height = 600)
    ylim_min = min(rbind(IRF_kanzig2plot[,indshocked],IRF_sparse2plot[,indshocked],IRF_sparse_lwr2plot[,indshocked]))
    ylim_max = max(rbind(IRF_kanzig2plot[,indshocked],IRF_sparse2plot[,indshocked],IRF_sparse_upr2plot[,indshocked]))
    
    
    plot(IRF_kanzig2plot[,indshocked],type="l", cex.axis=2,cex.main=2,
         # NOME MAIN con ->
         main =sprintf("%s -> %s",col_names_X[indshock],
                       col_names_X[indshocked]),xlab="",ylab="", col = 'blue', lty=5,lwd=3,
         ylim = c(ylim_min, ylim_max))
    lines(IRF_sparse2plot[,indshocked], col='black',lty=1,lwd=3)
    
    lines(IRF_sparse_lwr2plot[,indshocked], type="l",col='red', lty=3)
    lines(IRF_sparse_upr2plot[,indshocked], type="l",col='red', lty=3)
    
    # polygon(c(0:max_horizon, rev(0:max_horizon)), 
    #         c(IRF_sparse_lwr2plot[,indshocked],rev(IRF_sparse_upr2plot[,indshocked])),
    #         col = "#6BD7AF")
    # 
    abline(h=0)
    dev.off()
    
  }
}





