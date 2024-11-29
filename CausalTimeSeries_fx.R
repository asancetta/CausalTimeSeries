####################
# Cordoni Francesco and Alessio Sancetta 19 May 2023
# script for implementation of the procedure of the paper 
# "Consistent Causal Inference for High Dimensional Time Series" Cordoni Francesco and Alessio Sancetta
#####################


#-------------------------------------------- Main Algorithm and Function

#---------------------Main requirements



# execute only the first time

skip_install_requirments = TRUE#only set true if previously installed

if (!skip_install_requirments)
{

requiredPackages = c("igraph","pcalg", "glasso","R.matlab", "huge",
                       "tictoc","vars","flare","Matrix","MASS","SID","ks",
                       "assert", "foreach", "doParallel","gRbase", "expm", "tsDyn")
# 
install.packages("BiocManager")
BiocManager::install(c("graph","RBGL", "Rgraphviz"))
#                      
#install.packages(requiredPackages)

install.packages("igraph")
install.packages("pcalg")
install.packages("glasso")
#install.packages("R.matlab")
install.packages("huge")
install.packages("tictoc")
install.packages("vars")
install.packages("flare")
install.packages("Matrix")
install.packages("MASS")
install.packages("SID")
install.packages("ks")
install.packages("assert")
install.packages("foreach")
install.packages("doParallel")
install.packages("gRbase")
install.packages("expm")
install.packages("tsDyn")


} 


require(igraph)
require(pcalg)
require(glasso)
#require(R.matlab)
require(huge)
require(tictoc)
require(vars)
require(flare)
library(Matrix)
require(SID)
require(ks)
require(assert)
library(foreach)
library(doParallel)
require(pracma)
require(gRbase)
require(expm)
require(tsDyn)


#-------------functions to estimate the latent VAR model

get_W<-function(Xall,p_lag){
  K  = dim(Xall)[2]
  W  = Xall[(p_lag+1):nrow(Xall),] 
  for( p in 1:p_lag){
    Xall_lag_p = Xall[(p_lag+1-p):(nrow(Xall)-p),]
    W =matrix(c(W,Xall_lag_p),nrow(Xall_lag_p),(p+1)*K)
  }
  return(W)
}

get_Omega_thresholded <- function (Omega, tau){
  
  OUT = NULL
  aux_Omega = Omega
  OUT$old_Omega = Omega
  diag(aux_Omega)<-0
  ix  = aux_Omega!=0
  ix_thresholded = ( abs(aux_Omega) >= tau )
  ix_thresholded = ix_thresholded *1
  OUT$Omega = aux_Omega * ix_thresholded
  diag(  OUT$Omega) <- diag( OUT$old_Omega)
  OUT$Omega =  ( OUT$Omega + t( OUT$Omega ) )*0.5
  return(OUT)
  
}

get_Omega<- function(beta_hat, tau){
  Omega_hat      = get_Omega_thresholded(beta_hat, tau = tau)    
  Omega_hat      = Omega_hat$Omega
  return(Omega_hat)
}

get_Omega_from_mb_old<-function(mb_res, lambda, tau){
  
  ind4mb = which(mb_res$lambda== lambda)
  #ensure that mb was computed for the given lambda
  assert(length(ind4mb)==1,msg='lambda is not in mb_res$lambda in get_Omega_from_mb')
  beta_hat  = as.matrix(mb_res$beta[[ind4mb]])
  Omega     = get_Omega(beta_hat, tau)
  return(Omega)
}

get_Omega_from_mb<-function(mb_res, lambda, tau){
  mb_lambdas         = mb_res$lambda
  #if there are many lambdas choose the scale for comparison
  if(length(mb_lambdas)>1){
    min_lambda_decay   =  abs(diff(mb_lambdas))/(mb_lambdas[1:(length(mb_res$lambda)-1)]+0.00001)
    min_lambda_decay   = min(min_lambda_decay)  
  }else{min_lambda_decay = mb_lambdas}#if there is only one lambda in mb, use the first for scale in comparison
  
  cond4mb  =  abs(mb_res$lambda-lambda)/lambda<0.01*min_lambda_decay
  ind4mb = which(cond4mb)
  #ensure that mb was computed for the given lambda
  assert(length(ind4mb)==1,msg='lambda is not in mb_res$lambda in get_Omega_from_mb')
  beta_hat  = as.matrix(mb_res$beta[[ind4mb]])
  Omega     = get_Omega(beta_hat, tau)
  return(Omega)
}



get_Theta <- function(Sigma, Omega){
  OUT = NULL
  OUT$refit = 1*(Omega !=0)
  
  B_non_zero  = OUT$refit
  
  N = dim(Sigma)[1]
  
  Eye = diag(N)
  Theta = matrix(0,N,N)
  #wi1 be the vector of non-zero elements of wi
  #, and then define (Bi)p×s0
  #as a known matrix with either 0’s
  #or 1’s such that Bi * wi1 = wi
  
  for (ii in 1:N){
    
    aux_path = as.matrix( B_non_zero ) +diag(N) # obsolete -  this was necessary in a previous version of the code 
    # to ensure that in the diagonal there is a non zero element
    auxB = aux_path[,ii]
    Bi = diag(N); Bi = Bi[,which(auxB != 0)]
    
    Theta[,ii] = Bi%*% solve( t(Bi) %*% Sigma %*% Bi) %*% t(Bi) %*% Eye[,ii]
    
  }
  Theta = (t(Theta)+Theta)/2
  return(Theta)
}


get_Sigma_eps<-function(Theta,K){
  Theta_11  = Theta[1:K, 1:K]
  Sigma_eps = solve(Theta_11)
  return(Sigma_eps)
}


get_VAR_params<-function(Theta,K, Sigma_eps =NULL){
  if( is.null(Sigma_eps))
  {
    Theta_11  = Theta[1:K, 1:K]
    Sigma_eps = solve(Theta_11)
  }
  A = - Sigma_eps %*% Theta[1:K, (K+1):ncol(Theta)]
  out = list(A = A, Sigma_eps = Sigma_eps)
  return(out)
}

get_neg_likelihood <- function(Sigma, Omega) {
  ot = as.numeric(unlist(determinant(Omega)))
  #if (ot[2]<=0) warning("Precision matrix estimate is not positive definite!")
  tmp = (sum(diag(Sigma%*%Omega))  - ot[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    return(Inf)
  }
}


get_mb_aic<- function(Sigma_hat,mb_res,n, tau_rel_set =c(0))
{
  lambdas = mb_res$lambda
  nlambda = length(lambdas)
  K_W = nrow(mb_res$df)#number of variables/equations
  aic    = matrix(c(0),nlambda,length(tau_rel_set)) 
  for (l in 1:nlambda){
    for (k in 1:length(tau_rel_set)){
      tau        = tau_rel_set[k]*lambdas[[l]]
      beta_hat   = as.matrix(mb_res$beta[[l]])
      Omega_hat  = get_Omega(beta_hat, tau)
      #number of non zero parameters accounting for symmetry
      df         = K_W+((sum(Omega_hat!=0)-K_W))/2
      Theta_hat  = get_Theta(Sigma_hat, Omega_hat)
      aic[l,k]   = get_neg_likelihood(Sigma_hat, Theta_hat)+2*df/n
    }
  }
  return(aic)
}


#------------functions for choosing lambda and tau via cross-validation

part.cv <- function(n, fold) {
  
  ntest = floor(n/fold)
  ntrain = n-ntest
  
  ind = sample(n)
  
  trainMat = matrix(NA, nrow=ntrain, ncol=fold)
  testMat = matrix(NA, nrow=ntest, ncol=fold)
  
  nn = 1:n
  
  for (j in 1:fold) {
    sel = ((j-1)*ntest+1):(j*ntest)
    testMat[,j] = ind[sel ]
    sel2 =nn[ !(nn %in% sel) ]
    trainMat[,j] = ind[sel2]
  }
  
  return(list(trainMat=trainMat, testMat=testMat))
}


get_lambda_tau_set<-function(lambdas4cv, tau_rel_set){
  # lambdas4cv     = lambdas[(ind_lambda_aic-n_lambdas_cv_interval)
  #                          :(ind_lambda_aic+n_lambdas_cv_interval)]
  n_lambdas4cv   = length(lambdas4cv)
  n_tau4cv       = length(tau_rel_set)
  lambda_tau_set = matrix(c(rep(lambdas4cv,n_tau4cv),kronecker(tau_rel_set,lambdas4cv)),
                          n_lambdas4cv*n_tau4cv,n_tau4cv)
  return(lambda_tau_set)
}


get_obj_cv_util<-function(W_est,W_test,lambda, tau){
  
  K_W            = ncol(W_test)
  Sigma_hat_test = huge.npn(W_test, npn.func="skeptic",verbose = FALSE)
  Sigma_hat_est  = huge.npn(W_est, npn.func="skeptic",verbose = FALSE)
  mb_res_est     = huge(Sigma_hat_est,lambda = lambda,method = "mb",verbose = FALSE)
  beta_hat_est   = mb_res_est$beta[[1]]
  Omega_hat_est  = get_Omega(beta_hat_est, tau)
  #number of non zero parameters accounting for symmetry
  df             = K_W+((sum(Omega_hat_est!=0)-K_W))/2
  Theta_hat_est  = get_Theta(Sigma_hat_est, Omega_hat_est)
  
  neg_loglik    = get_neg_likelihood(Sigma_hat_test, Theta_hat_est)
  out           = NULL
  out$obj       = neg_loglik
  out$df        = df
  return(out)
}

get_obj_cv<- function(seed_num, fold, W, lambda,tau){
  # evaluation of CV function using fold for 1 value of the grid
  
  n   = nrow(W)
  K_W = ncol(W) 
  
  set.seed(seed_num) 
  ind_cv = part.cv(n, fold)
  
  loss.re = matrix(0, nrow = fold, ncol = 1)
  number_non_zero_theta = matrix(0, nrow = fold, ncol = 1)
  for (i in 1:fold) {
    W_est  = W[ind_cv$trainMat[,i], ]
    W_test = W[ind_cv$testMat[,i], ]
    
    obj_cv = get_obj_cv_util(W_est,W_test,lambda, tau)
    
    loss.re[i] = loss.re[i]  + obj_cv$obj
    number_non_zero_theta[i] = obj_cv$df
  }
  
  loss.mean = mean(loss.re)
  loss.se = sd(loss.re)/sqrt(fold)
  number_non_zero_theta = mean(number_non_zero_theta)
  return(c(loss.mean,loss.se,number_non_zero_theta ))
}

get_lambda_tau_cv<-function(seed_num_cv, fold_cv, W, lambda_tau_set){
  #cross validation to choose lambda and tau
  n_lambda_tau_set = nrow(lambda_tau_set)
  obj_cv           = matrix(c(0), n_lambda_tau_set)
  se_obj_cv        = matrix(c(0), n_lambda_tau_set)
  for(i in 1:n_lambda_tau_set){
    lambda    = lambda_tau_set[i,1]
    tau       = lambda_tau_set[i,2] 
    res_cv    = get_obj_cv(seed_num_cv, fold_cv, W, lambda,tau)
    obj_cv[i] = res_cv[1] 
    se_obj_cv[i] = res_cv[2] 
    
    message2print = paste('progress equal to ',
                          floor(100*i/n_lambda_tau_set), "%")
    print(message2print)
  }
  
  ind_cv    = which.min(obj_cv) 
  lambda_cv = lambda_tau_set[ind_cv,1]
  tau_cv    = lambda_tau_set[ind_cv,2] 
  out = list(lambda= lambda_cv, tau=tau_cv, obj = obj_cv)
  return(out)
}


get_ind4lambda_tau_sets<-function(lambdas, tau_rel_set, lambda_cv,tau_cv,
                                  flag_minimise_false_neg = FALSE){
  # choose a smaller penalty to avoid false negatives: important in PC
  
  ind_lambda_tmp   = which(lambda_cv== lambdas)
  ind_tau_tmp      = which(abs((tau_cv/lambda_cv)-tau_rel_set)<=1e-4)  
  if (flag_minimise_false_neg== TRUE){
    
    ind_lambda_tmp = which(lambda_cv== lambdas)
    ind_tau_tmp = which(abs((tau_cv/lambda_cv)-tau_rel_set)<=1e-4)
    if( ind_lambda_tmp<length(lambdas))
    {
      lambda_cv = lambdas[ind_lambda_tmp+1]
      tau_cv = tau_rel_set[ind_tau_tmp]*lambda_cv
    }else{ print('WARNING: lambda_cv on the boundary')}  
  }
  out = list(lambda = lambda_cv,tau = tau_cv, 
             ind4lambdas =ind_lambda_tmp,  
             ind4tau_rel_set =ind_tau_tmp)
}


#----- functions to identify the causal graph

get_pc_alg <- function(Sigma_eps,n, alpha=0.01, zero_restr = NULL){
  # if not null, zero_restr must be symmetric and contrain TRUE where no edge is present 
  K = nrow(Sigma_eps)
  pc_res <-pc(suffStat = list(C =  cov2cor(Sigma_eps) , n = n),
              indepTest = gaussCItest, ## (partial correlation test)
              alpha = alpha, p = K,
              fixedGaps = zero_restr )
  #fixedGaps logical symmetric matrix of dimension p*p. If entry [i,j] is TRUE, the result is
  #guaranteed to have no edge between nodes i and j.
  
  adj_mat <- t(as(pc_res,"amat"))
  
  # we take the transpose since pcalg report the edges by row
  # recall that the PCalgo package produce adjacent matrix which are the transposed of the usual adj. matrix.
  # indeed (as done by igrpah pacakge, used to plot a graph) A[i,j] = 1 means   i-->j (parents(j) are in the j column) , while paclgo provide
  # as(*, "amat") amat[a,b] = 0  and  amat[b,a] = 1   implies a --> b, therefore in the internal function we have used the transpose
  # try this as example
  # g_directed_all
  # plot(g_directed_all)
  # g # this is what we use and what igrpah uses
  # as(g_directed_all,"amat")
  # graph_from_adjacency_matrix(g)
  # plot(graph_from_adjacency_matrix(g)) 
  return(list(results_pc =pc_res, adj_mat  = adj_mat ))
}


# --------- functions to find SVAR params

get_topological_order_from_graph<-function (adj_mat){
  #adj_mat is an adjacency matrix - $adj_mat output from get_pc_alg 
  # A topological ordering of a directed graph is a linear ordering of its vertices such that, for
  # every edge (u → v), u comes before v in the ordering. A topological ordering is possible if
  # and only if the graph has no directed cycles, that is, if it is a directed acyclic graph (DAG).
  # Any DAG has at least one topological ordering
  #  PI is the permutation matrix s.t.  the i th element of PI eps_t is not a parent of i-1
  require(gRbase)
  ix  = topo_sortMAT(adj_mat, index = TRUE)
  N = ncol(adj_mat)
  
  PI=  (diag(N))
  PI = PI[ix,]
  
  new_graph = adj_mat[ix,ix]
  
  return(list(DAG_oriented = new_graph, topological_order = ix,
              PI = PI))
}


get_PI<-function(adj_mat){
  res_topo = get_topological_order_from_graph(adj_mat)
  PI = res_topo$PI  
  return(PI)
}

get_SVAR_params<- function(Sigma_eps, adj_mat){
  
  PI       = get_PI(adj_mat)
  
  K = ncol(Sigma_eps)
  
  Delta = matrix(0,nrow = K,ncol = K)
  for(i in 1:K){
    V_i = which(adj_mat[,i] == 1)*1
    #B_i = union(V_i,c((K+1):(2*K)))
    if (length(V_i)>0){
      B_i = V_i
      Sigma_Bi_Bi = Sigma_eps[B_i,B_i]
      Sigma_Bi_i = Sigma_eps[B_i,i]
      d_i = solve(Sigma_Bi_Bi) %*%  Sigma_Bi_i
      
      
      
      Delta[i, V_i] = t(d_i[1:length(V_i)])
    }
  }
  
  D =  PI %*% Delta  %*% solve(PI)  
  #all(D[upper.tri(D)]==0)
  H1 = diag(K) - D 
  H = solve(H1)
  return(list(D = D, H = H, PI=PI))
}


#--------------------- functions for impulse response

# #this just returns normal quantiles: wrong!
# get_f_k_inverse_old <- function(X){
#   K           = ncol(X)
#   f_k_inverse = list()
#   for (k in 1:K){
#     # Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages, American Statistician 50, 361--365. 10.2307/2684934.
#     # recommended type 8. 
#     # use the default of type = 7
#     f_k_inverse[[k]] = function(x) quantile(X[,k],probs=pnorm(x), type=7)
#   }
#   return(f_k_inverse)
# }
# 

get_f_k_inverse <- function(X){
  K           = ncol(X)
  f_k_inverse = lapply(1:K, function(k){
    # Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages, American Statistician 50, 361--365. 10.2307/2684934.
    # recommended type 8. 
    # use the default of type = 7
    force(k)
    function(x){ quantile(X[,k],probs=pnorm(x), type=7)}
  })
  return(f_k_inverse)
}


get_Z2X<-function(Z,f_k_inverse ){
  #f_k_inverse is the output of get_f_k_inverse(X)
  K = ncol(Z)
  out = matrix(0,nrow(Z),ncol(Z))
  for (k in 1:K){
    out[,k]   = f_k_inverse[[k]](Z[,k])
  }
  return(out)
}


get_svar_params_companion <-function(Sigma_eps, Theta, H, PI, K){
  # Parameter for SVAR in companion form 
  #theta is the precision matrix of W
  K_W    = nrow(Theta)
  p_lag  = K_W/K-1 
  A      = - Sigma_eps %*% Theta[1:K,(K+1):K_W]
  if (p_lag!=1){
    
    A_companion = matrix(0,K*p_lag,K*p_lag)
    A_companion[1:K,] = A 
    n_blocks <- p_lag-1
    A_companion[(K+1):(K*p_lag),1:(K*(p_lag-1))] = diag(K*(p_lag-1))
    
    A_original = A
    A = A_companion
    # augment H for the companion we can consider only the first row and columns
    H_original = H
    PI_original = PI
    H_augmented = diag(K*p_lag)
    H_augmented[1:K, 1:K] = H
    PI_augmented  = diag(K*p_lag)
    PI_augmented[1:K, 1:K] = PI
    
    PI = PI_augmented
    H = H_augmented
  }else{
    PI_original = PI
    H_original = H
  }
  
  return(list(H = H , A = A, PI = PI,
              H_original = H_original, PI_original = PI_original
  ))
}


# get_IRF_Z_old<-function(H,PI,A, max_horizon, Lag=1){
#   
#   UPsi = list()
#   UPsi[[1]] = solve(PI)  %*% H #this is Psi_0
#   
#   N = ncol(A)
#   
#   for (i in 2:(max_horizon+1)){
#     aux = matrix(0,nrow = N, ncol = N )
#     j = 2;  
#     #-->  
#     # Psi_i = sum_{j=1}^i A_j Psi_{i-j},  i = 1,2,.., max_horizon,   where A_j = 0 , when j>p
#     # but since we start to count from 1 and not from 0 we have 
#     # (The list of Psi start from 1 )
#     # Psi_i = sum_{j=2}^i A_{j-1} Psi_{i-j} , i = 2,3,.., max_horizon+1,  where A_j = 0 , when (j-1)>p
#     # 
#     while ( (j-1)<= Lag ){
#       aux  = aux + A %*% UPsi[[ i-j+1 ]]
#       j = j + 1;
#     }
#     UPsi[[ i ]] = aux
#   }
#   
#   
#   return(UPsi)
# }


get_IRF_Z<-function(H,PI,A, max_horizon, Lag=1){
  
  UPsi = list()
  UPsi[[1]] = t(PI)%*% H%*%PI #this is Psi_0
  
  N = ncol(A)
  
  for (i in 2:(max_horizon+1)){
    aux = matrix(0,nrow = N, ncol = N )
    j = 2;  
    #-->  
    # Psi_i = sum_{j=1}^i A_j Psi_{i-j},  i = 1,2,.., max_horizon,   where A_j = 0 , when j>p
    # but since we start to count from 1 and not from 0 we have 
    # (The list of Psi start from 1 )
    # Psi_i = sum_{j=2}^i A_{j-1} Psi_{i-j} , i = 2,3,.., max_horizon+1,  where A_j = 0 , when (j-1)>p
    # 
    while ( (j-1)<= Lag ){
      aux  = aux + A %*% UPsi[[ i-j+1 ]]
      j = j + 1;
    }
    UPsi[[ i ]] = aux
  }
  
  
  return(UPsi)
}

#########################################

get_IRF_Z_simple<-function(H,PI,A, max_horizon, Lag=1){
  #to be used wiht VAR1 only (or VAR1 companion form of a VARp)
  UPsi = list()
  UPsi[[1]] = t(PI)%*% H%*%PI #this is Psi_0
  N = ncol(A)
  for (i in 2:(max_horizon+1)){
    UPsi[[ i ]] = A %*% UPsi[[ i-1 ]]
  }
  return(UPsi)
}

############################################

get_conditional_VAR_from_summands<-function(AsZ_samples,UPsiXi,UPsiXiShocked,horizon){
  
  first_term   = AsZ_samples[[horizon+1]]
  last_term    = UPsiXiShocked[[horizon+1]]
  middle_terms = 0*last_term
  if (horizon>0){
    for(s in 1:(horizon+1)){
      middle_terms =middle_terms+UPsiXi[[s]]    
    }
  }
  Z_shocked = first_term+middle_terms+last_term
  return(Z_shocked)
}


# 
# get_IRF_X_old <-function(m = 500, adj_mat, Sigma_eps,
#                      Theta, f_k_inverse, max_horizon = 50,
#                      k_shock_variable = 1, delta_shock = 1, seed_num=1, flag_shock_Z=FALSE){
#   
#   #k_shock_variable - index of variable shock
#   #delta_shock - size of shock
#   
#   K           = ncol(Sigma_eps)
#   K_W         = nrow(Theta) 
#   Sigma       = solve(Theta)
#   
#   svar_params = get_SVAR_params(Sigma_eps, adj_mat)
#   PI          = svar_params$PI
#   D           = svar_params$D
#   H           = svar_params$H
#   
#   # get the parameter of the VAR in companion form
#   svar_params_comp = get_svar_params_companion(Sigma_eps, Theta, H, PI,  K)
#   # extract the parameters
#   PI_comp = svar_params_comp$PI
#   A_comp = svar_params_comp$A
#   H_comp = svar_params_comp$H
#   
#   # IRFs on Z
#   UPsi_comp = get_IRF_Z(H_comp,PI_comp, A_comp, max_horizon, Lag=1)
#   
#   #Covariance matrix of shocks
#   # solve(H) %*% PI %*% Sigma_eps %*% t(solve(H) %*% PI)
#   Sigma_xi_pi = solve(H) %*% PI %*% Sigma_eps %*% t(solve(H) %*% PI) 
#   # ?? this matrix is PI Sigma_xi t(PI), will be used only to sample  the shock on row 583
#   # to avoid numerical rounding error we impose that the out of diagonal elements of Sigma_xi are equal to zero
#   
#   Sigma_xi_pi = diag(diag(Sigma_xi_pi))
#   
#   # Simulate latent vector Z for MC integration
#   if (!seed_num == "NA"){
#     set.seed((m+1)*seed_num) 
#   }
#   Sigma_Z_samples = Sigma[(K+1):K_W, (K+1):K_W]
#   eS <- eigen(Sigma_Z_samples, symmetric = TRUE)
#   ev <- eS$values
#   if ( !all(ev >0 )){
#     # possible numeric error, therefore we fix all the eigenvalues which are  negative to the smallest positive egeinvalues
#     iev_neg <- which(ev<=0)
#     ev[iev_neg] = min(ev[which(ev>0)])
#     Sigma_Z_samples = eS$vectors %*% diag(ev) %*%solve(eS$vectors) 
#   } 
#   Z_samples = MASS::mvrnorm(n = m, mu =  rep(0, K*p_lag), 
#                             Sigma=Sigma_Z_samples)
#   
#   
#   #A^s*Z
#   AsZ_samples  = list()
#   #recursively construct A_comp^s*Z_comp, avoiding unnecessary calculations
#   s   = 1
#   As1 = A_comp
#   AsZ_samples[[1]]   = Z_samples%*%t(As1[1:K,])
#   As0 = As1
#   for(s in 2:(max_horizon+1)){
#     As1 = As0%*%A_comp
#     AsZ_samples[[s]]   = Z_samples%*%t(As1[1:K,])
#     As0 = As1
#   }
#   
#   #find the index corresponding to the desired variable that is shocked (variable k_shock_variable) after permutation
#   e_shock = matrix(0, nrow = K, ncol = 1)
#   # H[k,k] = 1 by definition, but may divide by std of shock
#   e_shock[k_shock_variable] = delta_shock/H[k_shock_variable,k_shock_variable] 
#   e_shock = PI %*% e_shock# this will be used to identify the correct shock ?? in the old  code we had already used the permutation
#   
#   
#   #Upsilon_s*xi
#   UPsiXi      = list()
#   #Upsilon_s*xiDelta
#   UPsiXiDelta = list()
#   #Upsilon_s*xi0
#   UPsiXi0    = list()
#   for (s in 1:(max_horizon+1)){
#     #simulate the structural innovations for MC integration
#     if (!seed_num == "NA"){
#       set.seed(s*seed_num)
#     }
#     
#     xi_samples = MASS::mvrnorm(n = m , mu =  rep(0, K), Sigma=Sigma_xi_pi)
#     
#     #set the shocked innovation to delta
#     xi_samples_delta = xi_samples
#     xi_samples_delta[,which(e_shock!=0)] = delta_shock 
#     
#     #set the shocked innovation to zero
#     xi_samples_0 = xi_samples
#     xi_samples_0[,which(e_shock!=0)] = 0 
#     
#     #compute the s^th summand in the MA representation of VAR; only need to first K times K entries
#     UPsiXi[[s]]      = xi_samples%*%t(UPsi_comp[[s]][1:K,1:K])  
#     UPsiXiDelta[[s]] = xi_samples_delta%*%t(UPsi_comp[[s]][1:K,1:K])  
#     UPsiXi0[[s]]     = xi_samples_0%*%t(UPsi_comp[[s]][1:K,1:K])  
#   }
#   
#   IRF_X = matrix(0,nrow=max_horizon+1,ncol=K)
#   for (s in 1:(max_horizon+1)){
#     horizon     = s-1
#     Z_delta    = get_conditional_VAR_from_summands(AsZ_samples,UPsiXi,UPsiXiDelta,horizon)
#     Z_0         = get_conditional_VAR_from_summands(AsZ_samples,UPsiXi,UPsiXi0,horizon)
#     X_delta     = get_Z2X(Z_delta, f_k_inverse)
#     X_0         = get_Z2X(Z_0, f_k_inverse)
#     
#     if (flag_shock_Z){
#       IRF_X[s,]   = colMeans(Z_delta)-colMeans(Z_0) # ?? IRF on Z 
#     }else{
#       IRF_X[s,]   = colMeans(X_delta)-colMeans(X_0) 
#     }
# 
#   }  
#   return(IRF_X)
# }


get_IRF_X <-function(m = 500, adj_mat, Sigma_eps,
                     Theta, f_k_inverse, max_horizon = 50,
                     k_shock_variable = 1, delta_shock = 1, seed_num=0, 
                     flag_shock_Z=FALSE, flag_cum = FALSE){
  
  #k_shock_variable - index of variable shock
  #delta_shock - size of shock
  #flag_cum TRUE to compute cumulative irfs
  
  K           = ncol(Sigma_eps)
  K_W         = nrow(Theta) 
  Sigma       = solve(Theta)
  
  svar_params = get_SVAR_params(Sigma_eps, adj_mat)
  PI          = svar_params$PI
  D           = svar_params$D
  H           = svar_params$H
  
  # get the parameter of the VAR in companion form
  svar_params_comp = get_svar_params_companion(Sigma_eps, Theta, H, PI,  K)
  # extract the parameters
  PI_comp = svar_params_comp$PI
  A_comp = svar_params_comp$A
  H_comp = svar_params_comp$H
  # IRFs on Z
  UPsi_comp = get_IRF_Z(H_comp,PI_comp, A_comp, max_horizon, Lag=1)
  
  #Covariance matrix of shocks
  # solve(H) %*% PI %*% Sigma_eps %*% t(solve(H) %*% PI)
  
  invH     = solve(H)
  PIinvHPI = t(PI)%*%invH %*% PI
  Sigma_xi = PIinvHPI%*% Sigma_eps %*% t(PIinvHPI) 
  Sigma_xi = diag(diag(Sigma_xi))# may not be diag due to numerical rounding
  
  # Simulate latent vector Z for MC integration
  if (!seed_num == "NA"){
    set.seed((m+1)*seed_num) 
  }
  Sigma_Z_samples = Sigma[(K+1):K_W, (K+1):K_W]
  eS <- eigen(Sigma_Z_samples, symmetric = TRUE)
  ev <- eS$values
  if ( !all(ev >0 )){
    # possible numeric error, therefore we fix all the eigenvalues which are  negative to the smallest positive egeinvalues
    iev_neg <- which(ev<=0)
    ev[iev_neg] = min(ev[which(ev>0)])
    Sigma_Z_samples = eS$vectors %*% diag(ev) %*%solve(eS$vectors) 
  } 
  Z_samples = MASS::mvrnorm(n = m, mu =  rep(0, K*p_lag), 
                            Sigma=Sigma_Z_samples)
  
  
  #A^s*Z
  AsZ_samples  = list()
  #recursively construct A_comp^s*Z_comp, avoiding unnecessary calculations
  s   = 1
  As1 = A_comp
  AsZ_samples[[1]]   = Z_samples%*%t(As1[1:K,])
  As0 = As1
  for(s in 2:(max_horizon+1)){
    As1 = As0%*%A_comp
    AsZ_samples[[s]]   = Z_samples%*%t(As1[1:K,])
    As0 = As1
  }
  
  
  #find the index corresponding to the desired variable that is shocked (variable k_shock_variable) after permutation
  e_shock = matrix(0, nrow = K, ncol = 1)
  e_shock[k_shock_variable] = delta_shock
  
  #Upsilon_s*xi
  UPsiXi      = list()
  #Upsilon_s*xiDelta
  UPsiXiDelta = list()
  #Upsilon_s*xi0
  UPsiXi0    = list()
  for (s in 1:(max_horizon+1)){
    #simulate the structural innovations for MC integration
    if (!seed_num == "NA"){
      set.seed(s*seed_num)
    }
    
    xi_samples = MASS::mvrnorm(n = m , mu =  rep(0, K), Sigma=Sigma_xi)
    
    #set the shocked innovation to delta
    xi_samples_delta = xi_samples
    xi_samples_delta[,which(e_shock!=0)] = delta_shock 
    
    #set non-shocked innovations to zero
    xi_samples_0 = xi_samples
    xi_samples_0[,which(e_shock!=0)] = 0 
    
    
    #compute the s^th summand in the MA representation of VAR; only need to first K times K entries
    UPsiXi[[s]]      = xi_samples%*%t(UPsi_comp[[s]][1:K,1:K])  
    UPsiXiDelta[[s]] = xi_samples_delta%*%t(UPsi_comp[[s]][1:K,1:K])  
    UPsiXi0[[s]]     = xi_samples_0%*%t(UPsi_comp[[s]][1:K,1:K])  
  }
  
  IRF_X = matrix(0,nrow=max_horizon+1,ncol=K)
  for (s in 1:(max_horizon+1)){
    horizon     = s-1
    Z_delta    = get_conditional_VAR_from_summands(AsZ_samples,UPsiXi,UPsiXiDelta,horizon)
    Z_0         = get_conditional_VAR_from_summands(AsZ_samples,UPsiXi,UPsiXi0,horizon)
    #Z_0         = get_conditional_VAR_from_summands(AsZ_samples,UPsiXi,UPsiXi,horizon)
    X_delta     = get_Z2X(Z_delta, f_k_inverse)
    X_0         = get_Z2X(Z_0, f_k_inverse)
    
    if (flag_shock_Z){
      IRF_X[s,]   = colMeans(Z_delta)-colMeans(Z_0) # ?? IRF on Z
      
    }else{
      IRF_X[s,]   = colMeans(X_delta)-colMeans(X_0) 
    }
    
    
  }
  if (flag_cum == TRUE){
    IRF_X = apply(IRF_X, 2, cumsum)
  }
  
  return(IRF_X)
}


get_Z_samples_bootstrap <- function(n, K, p_lag, Sigma, Sigma_eps, A_comp,seed_num=0){
  
  # tic() 
  K_W         = nrow(Sigma)
  Z_samples = matrix(0,nrow = n+1, ncol = K*p_lag )
  
  Sigma_Z_samples = Sigma[(K+1):K_W, (K+1):K_W]
  eS <- eigen(Sigma_Z_samples, symmetric = TRUE)
  ev <- eS$values
  if ( !all(ev >0 )){
    # possible numeric error, therefore we fix all the eigenvalues which are  negative to the smallest positive egeinvalues
    iev_neg <- which(ev<=0)
    ev[iev_neg] = min(ev[which(ev>0)])
    Sigma_Z_samples = eS$vectors %*% diag(ev) %*%solve(eS$vectors) 
  } 
  set.seed(seed_num)
  Z_samples[1,] = MASS::mvrnorm(n = 1, mu =  rep(0, K*p_lag),
                                Sigma=Sigma_Z_samples)
  
  set.seed(seed_num+1)
  eps = mvtnorm::rmvnorm(n+1, mean = rep(0,K), sigma = Sigma_eps)
  eps_companion = matrix(0, nrow = n+1, ncol = K*p_lag)
  eps_companion[,1:K] = eps
  
  for (tt in 2:(n+1)){
    Z_samples[tt,] = A_comp %*% Z_samples[tt-1,] +  eps_companion[tt,]
  }
  
  return(Z_samples)
}



get_X_samples_bootstrap <- function(n, K, p_lag, Sigma, Sigma_eps, A_comp, f_k_inverse,seed_num=0){
  
  # tic() 
  K_W         = nrow(Sigma)
  Z_samples = matrix(0,nrow = n+1, ncol = K*p_lag )
  
  Sigma_Z_samples = Sigma[(K+1):K_W, (K+1):K_W]
  eS <- eigen(Sigma_Z_samples, symmetric = TRUE)
  ev <- eS$values
  if ( !all(ev >0 )){
    # possible numeric error, therefore we fix all the eigenvalues which are  negative to the smallest positive egeinvalues
    iev_neg <- which(ev<=0)
    ev[iev_neg] = min(ev[which(ev>0)])
    Sigma_Z_samples = eS$vectors %*% diag(ev) %*%solve(eS$vectors) 
  } 
  set.seed(seed_num)
  Z_samples[1,] = MASS::mvrnorm(n = 1, mu =  rep(0, K*p_lag),
                                Sigma=Sigma_Z_samples)
  
  set.seed(seed_num+1)
  eps = mvtnorm::rmvnorm(n+1, mean = rep(0,K), sigma = Sigma_eps)
  eps_companion = matrix(0, nrow = n+1, ncol = K*p_lag)
  eps_companion[,1:K] = eps
  
  for (tt in 2:(n+1)){
    Z_samples[tt,] = A_comp %*% Z_samples[tt-1,] +  eps_companion[tt,]
  }
  # toc()
  # 0.149 sec elapse
  #---------------
  # recover the X
  # tic()
  # Z_samples = VAR.sim(B = A_comp[1:K, ], lag = 12, n = n, include = "none", varcov = Sigma_eps)
  # toc()
  #  0.388 second my code is faster thanks to the companion form
  X_samples     = get_Z2X(Z_samples[,1:K], f_k_inverse)
  
  return(X_samples)
}


get_X_samples_block_bootstrap <- function(X, day_number,seed_num =0){
  # block bootstrap over day
  # tic() 
  day = unique(day_number)
  T_day = length(day)
  
  
  #sample daily  
  set.seed(seed_num)
  iboot = sample(1:T_day,   replace = TRUE)
  selected_day = day[iboot]
  ix_day = NULL
  for (jj in 1: T_day ){
    ix_day = c(ix_day, which(day_number == selected_day[jj]) )
  }
  X_samples =X[ ix_day,]
  
  return(X_samples)
}


get_IRFs_bootstrap <- function(mboot = 200, n, m, Sigma, Sigma_eps, A_comp,
                               f_k_inverse, adj_mat,
                               p_lag, K, max_horizon, k_shock_variable,
                               lambda_cv_adj, tau_cv_adj, seed_num = 0, delta_shock=1,
                               flag_shock_Z=FALSE, flag_cum = FALSE){
  # this function will compute IRFs for a selected k_shock_variable
  # on all the other variables
  #set.seed(seed)
  # IRFs_X_boot = array(0,c(max_horizon+1, K, mboot) )
  
  
  
  # parallalized 
  num_cores <- detectCores()
  registerDoParallel(cores = num_cores)
  IRFs_X_boot <- foreach(kboot = 1:mboot, .combine = "cbind") %dopar% {
    X_samples = get_X_samples_bootstrap (n, K, p_lag, Sigma, Sigma_eps, A_comp, f_k_inverse, (seed_num+2)*kboot)
    
    # Algo 1
    W         = get_W(X_samples,p_lag)
    K_W       = ncol(W)
    Sigma_hat_boot = huge.npn(W, npn.func="skeptic",verbose = FALSE)
    
    # Algorithm 2 and 3
    # we do not perform cv again in the interest of time 
    
    mb_res  = huge(Sigma_hat_boot , lambda =lambda_cv_adj ,method = "mb",verbose = FALSE)
    Omega_hat_boot      = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
    Theta_hat_boot      = get_Theta(Sigma_hat_boot, Omega_hat_boot)
    Sigma_eps_hat_boot  = get_Sigma_eps(Theta_hat_boot,K)
    
    # ------------------------------------------------------------ get the IRFs
    
    f_k_inverse_boot = get_f_k_inverse(X_samples)
    
    IRFs_X <- get_IRF_X(m = m, adj_mat, Sigma_eps_hat_boot,
                        Theta_hat_boot, f_k_inverse_boot, max_horizon = max_horizon,
                        k_shock_variable = k_shock_variable, delta_shock = delta_shock, seed_num=(seed_num+1)*kboot,
                        flag_shock_Z=flag_shock_Z, flag_cum = flag_cum)
    
    
    # Return IRFs_X for current iteration
    IRFs_X
  }
  # Convert the result to the desired array format
  IRFs_X_boot <- array(IRFs_X_boot, dim = c(max_horizon + 1, K, mboot))
  return(IRFs_X_boot)
  
  
}





get_IRFs_bootstrap_selected_shocks <- function(mboot = 200, n, m, Sigma, Sigma_eps, A_comp,
                                               f_k_inverse, adj_mat,
                                               p_lag, K, max_horizon, k_shock_variable,
                                               lambda_cv_adj, tau_cv_adj, seed_num = 0, k_shocked_variable,
                                               flag_block = FALSE, X, day_number,
                                               flag_shock_Z = FALSE, flag_cum  =FALSE,
                                               delta_shock = 1){
  
  # Please README
  # this function will compute IRFs for a vector of 'k_shock_variable'
  # on selected_variables 'k_shocked_variable' This will help saving time
  # for each bootstrap cycle the matrix 
  # IRFs_X will be a (horizon+1) x (k_shocked_variable *k_shock_variable)
  # suppose we have 
  # k_shock_variable=[2 , 3 , 5]
  # k_shocked_variable =[1, 10],
  # then
  # IRFs_X = [ IRFs_X of the 1st variable to be shocked, i.e., 1 to a shock on the 1st variable on k_shock_variable, i.e., 2,
  #  then we have IRFs_X of the 2st variable to be shocked, i.e., 10 to a shock on the 1st variable on k_shock_variable, i.e., 2,
  # IRFs_X of the 1st variable to be shocked, i.e., 1 to a shock on the 2st variable on k_shock_variable, i.e., 3,
  #  then we have IRFs_X of the 2st variable to be shocked, i.e., 10 to a shock on the 2st variable on k_shock_variable, i.e., 3, ....]
  #set.seed(seed)#UNCLEAR WHAT THIS DOES WHEN YOU LOOP
  # IRFs_X_boot = array(0,c(max_horizon+1, K, mboot) )
  
  
  
  # parallalized 
  num_cores <- detectCores()
  registerDoParallel(cores = num_cores)
  IRFs_X_boot <- foreach(kboot = 1:mboot, .combine = "cbind") %dopar% {
    
    if (flag_block){
      X_samples = get_X_samples_block_bootstrap (X, day_number,(seed_num+2)*kboot)
    }else{
      X_samples = get_X_samples_bootstrap (n, K, p_lag, Sigma, Sigma_eps, A_comp, f_k_inverse,(seed_num+2)*kboot)
    }
    
    # Algo 1
    W         = get_W(X_samples,p_lag)
    K_W       = ncol(W)
    Sigma_hat_boot = huge.npn(W, npn.func="skeptic",verbose = FALSE)
    
    # Algorithm 2 and 3
    # we do not perform again cv for time reasoning
    
    mb_res  = huge(Sigma_hat_boot , lambda =lambda_cv_adj ,method = "mb",verbose = FALSE)
    Omega_hat_boot      = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
    Theta_hat_boot      = get_Theta(Sigma_hat_boot, Omega_hat_boot)
    Sigma_eps_hat_boot  = get_Sigma_eps(Theta_hat_boot,K)
    
    # ------------------------------------------------------------ get the IRFs
    
    f_k_inverse_boot = get_f_k_inverse(X_samples)
    
    # copmute the IRF_X for all the k_shock_variable and save only those column of interest 
    # k_shocked_variable
    # Use lapply to iterate over k_shock_variable and extract the desired column
    IRFs_X_list <- lapply(k_shock_variable, function(i_k_shock_variable) {
      IRFs_X_i <- get_IRF_X(
        m = m, adj_mat, Sigma_eps_hat_boot, Theta_hat_boot, f_k_inverse_boot,
        max_horizon = max_horizon, k_shock_variable = i_k_shock_variable,
        delta_shock = delta_shock, seed_num = (seed_num+1)*kboot, flag_shock_Z=flag_shock_Z,
        flag_cum)
      
      # Extract the desired column from IRFs_X_i
      IRFs_X_i[, k_shocked_variable]
    })
    
    # Combine the results into a matrix using do.call and cbind
    IRFs_X_par <- do.call(cbind, IRFs_X_list)
    
    
    # IRFs_X_par = NULL
    # for (i_k_shock_variable in k_shock_variable){
    #   IRFs_X_i <- get_IRF_X(m = m, adj_mat, Sigma_eps_hat_boot,
    #                         Theta_hat_boot, f_k_inverse_boot, max_horizon = max_horizon,
    #                         k_shock_variable = i_k_shock_variable, delta_shock = 1, seed_num=1*kboot)
    #   IRFs_X_par = cbind(IRFs_X_par, IRFs_X_i[,c(k_shocked_variable)])
    # }
    
    
    
    
    # Return IRFs_X for current iteration
    IRFs_X_par 
  }
  # Convert the result to the desired array format
  IRFs_X_boot <- array(IRFs_X_boot, dim = c(max_horizon + 1,length(k_shock_variable)* length(k_shocked_variable), mboot))
  return(IRFs_X_boot)
  
  
}



# get_CI_boot_old <-function(IRFs_X_boot, alpha){
#   
#   IRFs_X_boot_mean = apply(IRFs_X_boot, c(1,2), mean)
#   mboot = dim(IRFs_X_boot)[3]
#   q_alpha = qnorm(1-alpha/2)
#   IRFs_X_boot_lwr = IRFs_X_boot_mean - q_alpha * apply(IRFs_X_boot, c(1,2), sd)/sqrt(mboot)
#   IRFs_X_boot_upr = IRFs_X_boot_mean + q_alpha * apply(IRFs_X_boot, c(1,2), sd)/sqrt(mboot)
#   
#   return(list(IRFs_X_boot_lwr = IRFs_X_boot_lwr, 
#               IRFs_X_boot_mean = IRFs_X_boot_mean, IRFs_X_boot_upr = IRFs_X_boot_upr))
# }

get_CI_boot <-function(IRFs_X_boot, alpha){
  
  IRFs_X_boot_mean = apply(IRFs_X_boot, c(1,2), mean)
  mboot = dim(IRFs_X_boot)[3]
  IRFs_X_boot_lwr = apply(IRFs_X_boot, c(1,2), quantile, probs = alpha/2)
  IRFs_X_boot_upr = apply(IRFs_X_boot, c(1,2), quantile, probs = 1-alpha/2)
  IRFs_X_boot_med = apply(IRFs_X_boot, c(1,2), quantile, probs = 1/2)
  
  return(list(IRFs_X_boot_lwr = IRFs_X_boot_lwr, 
              IRFs_X_boot_mean = IRFs_X_boot_mean, IRFs_X_boot_upr = IRFs_X_boot_upr,
              IRFs_X_boot_med = IRFs_X_boot_med))
}


get_CI_bootCentered <-function(IRFs_X_boot, alpha, IRFs_X_point_est){
  
  boot_quantiles   = get_CI_boot(IRFs_X_boot, alpha)
  shift4center     = IRFs_X_point_est-boot_quantiles$IRFs_X_boot_med
  IRFs_X_boot_lwr  = boot_quantiles$IRFs_X_boot_lwr+shift4center 
  IRFs_X_boot_upr  = boot_quantiles$IRFs_X_boot_upr+shift4center
  IRFs_X_boot_mean = boot_quantiles$IRFs_X_boot_mean
  
  return(list(IRFs_X_boot_lwr = IRFs_X_boot_lwr, 
              IRFs_X_boot_mean = IRFs_X_boot_mean,
              IRFs_X_boot_upr = IRFs_X_boot_upr))
}



############################end of tests


