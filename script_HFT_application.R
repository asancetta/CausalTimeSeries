####################
# Cordoni Francesco and Alessio Sancetta 19 May 2023
# script for implementation of the procedure of the paper 
# "Consistent Causal Inference for High Dimensional Time Series" Cordoni Francesco and Alessio Sancetta
#####################

# In this script you can find an R implementation of all the algorithm discussed in the paper
# and a step-by-step guide for using the procedure in an empirical application.

source("CausalTimeSeries_fx.R")

#################################################### 
#----------------------------------------------------------------------------------#
#                                                                   
# Load DATA and setting parameters                         
#----------------------------------------------------------------------------------#
namefile = 'Xall_volume_time.csv';#use full path if needed
 
# name of the file with intermediate results, use full path if needed
LocalSavedResults = "HFT_application_PC_result.RData"


flag_first_run = TRUE # set this to TRUE if it is the first run, 
#otherwise you can load the data from a local folder to save time for further analysis

# flag_instrument = FALSE
DATA = NULL

loading_from_local_save_results = FALSE # instead of running all the time you can save the data in local 
# and loading  that file, BUT THEN YOU HAVE TO RELOAD THE 
# FUNCTIONS IN THE HEADER OF THIS FILE IF YOU WANT TO RUN THE IRF


#M-B lasso graph calculation via Lasso or CLIME
flag_MB_CLIME = "MB"# Select CLIME for CLIME and MB for LASSO
#  (BE AWARE that estimation with CLIME in this application is quite slow,
# indeed you will  see little differences with LASSO graph and some edges will not be oriented)
# we suggest to use LASSO i.e. MB


if(!loading_from_local_save_results){
  aux_data = read.csv(namefile, header=TRUE)#read.csv(paste(InputDir,namefile,sep=""),header=TRUE)
   
  Xall = aux_data[,-1 ] # the first column is timestamp object
  
  DATA$Xall = as.matrix(Xall)
  DATA$names_nodes = as.character(colnames(DATA$Xall))
  #get the names
  
  DATA$rowX = aux_data[,1]
   
  
  #-------------------------------------------------------------
  #this should be the input
  X           = DATA$Xall
  col_names_X = DATA$names_nodes 
   
  
   
  #----------------------------settings
  p_lag = 1# number of lags in latent VAR
  # number of smoothing parameters to construct an initial grid of lambdas in Algorithm 2
  nlambda = 20
  # ratio between smallest over largest lamda parameter
  lambda_min_ratio = 0.01
  # set of values for tau divided by lambda, e.g. 2 means tau=2*lambda, in the paper 
  tau_rel_set = c(1.2,2)# a grid size of tau relative to lambda (tau = tau_rel_set*lambda) 
  fold_cv     = 5# number of folds for CV to choose smoothing parameters in Algorithm 2
  seed_num_cv = 1# seed for exact replication
  
  # number of lambdas to consider below and above lambda_aic
  #(lambda_aic is an estimated smoothing parameter around which we build a grid of lambdas for cross validation )
  n_lambdas_cv_interval = 3 
  flag_minimise_false_neg = TRUE# selects a smaller penalty to reduce false negatives in Algorithm 2
  alpha_pc    =0.05 # the alpha parameter in Algorithm 4
  thresh_stability =0.25# a parameter to select ???
  
  
  K = ncol(X) 
  n = nrow(X)
  
  
  
  #--------Algorithm 1
  
  W         = get_W(X,p_lag)
  K_W       = ncol(W)
  Sigma_hat = huge.npn(W, npn.func="skeptic",verbose = FALSE)
  
  #--------Algorithm 2 and 3 using cross-validation
  
  #M-B lasso graph calculation via Lasso
  
  if (flag_MB_CLIME == "MB"){
	  mb_res  = huge(Sigma_hat,nlambda = nlambda,lambda.min.ratio =lambda_min_ratio,method = "mb",verbose = FALSE)
   }else{
	  # for CLIME we can use the sugm function which has the same output of the huge, no need to change any other functions
    # it is slow compared to LASSO
	  mb_res  = sugm(Sigma_hat, nlambda = nlambda,lambda.min.ratio =lambda_min_ratio,method = "clime",verbose = FALSE)
   }
  
  lambdas = mb_res$lambda# the set of lambdas generated within huge for M-B procedure
  res_aic = get_mb_aic(Sigma_hat,mb_res,n, c(0))#get_mb_aic(mb_res, c(0))# thresh=0 for initial value of lambda and a set of lambdas
  #the index in lambdas that minimises AIC
  ind_lambda_aic = which.min(res_aic)
  #a subset of lambdas around the AIC minimiser to be used as grid in cross validation
   
  
  # select at least 7 elements for CV
  if ((ind_lambda_aic+n_lambdas_cv_interval)>length(lambdas)){
    lambdas4cv     = tail(lambdas, n=7)
  }else{
    if ((ind_lambda_aic-n_lambdas_cv_interval)<0){
      lambdas4cv     = lambdas[1:7]
    }else{
      lambdas4cv     = lambdas[(ind_lambda_aic-n_lambdas_cv_interval)
                               :(ind_lambda_aic+n_lambdas_cv_interval)]
    }
  }
  #the set of all allowed combinations of lambda and tau for cross validation
  lambda_tau_set = get_lambda_tau_set(lambdas4cv, tau_rel_set)
   
  if (flag_first_run){
    lambda_tau_cv = get_lambda_tau_cv(seed_num_cv, fold_cv, W, lambda_tau_set)
    save(lambda_tau_cv,file = "lambda_tau_cv_hft_empirical.RData")
  }else{
    load("lambda_tau_cv_hft_empirical.RData")
  }
  
  lambda_cv    = lambda_tau_cv$lambda
  tau_cv       = lambda_tau_cv$tau
  obj_cv       = lambda_tau_cv$obj
  #this plot shows as many curves as length(tau_rel_set) one on the right of the other
  plot(obj_cv)
  
  # lambda and tau, possibly adjusted to reduce false positives
  # and their index in lambdas and tau_rel_set above (the latter is not used)
  # if not adjusted, it just returns the input plus the inted in lambdas and tau_rel_set above
  lambda_tau_res = get_ind4lambda_tau_sets(lambdas, tau_rel_set, lambda_cv,tau_cv,
                                           flag_minimise_false_neg = flag_minimise_false_neg)
  
  lambda_cv_adj  = lambda_tau_res$lambda
  tau_cv_adj     = lambda_tau_res$tau
  #reuse the results from mb_res rather than recomputing
  Omega_hat      = get_Omega_from_mb(mb_res, lambda_cv_adj, tau_cv_adj)
  Theta_hat      = get_Theta(Sigma_hat, Omega_hat)
  Sigma_eps_hat  = get_Sigma_eps(Theta_hat,K)
  
  
  
  #----------------- Algorithm 5
  
  #let PC skip edges that are zero in Theta11, set to NULL otherwise
  zero_rest      = Theta_hat[1:K,1:K]==0
  pc_res         = get_pc_alg(Sigma_eps_hat,n, alpha=alpha_pc, 
                              zero_restr = zero_rest)
  
  iplotPC(pc_res$results_pc, labels = col_names_X)
  
  
  
  g = graph_from_adjacency_matrix(pc_res$adj_mat)
  #pc_res$adj_mat[which("SPYBookImb10"==DATA$names_nodes),which("CSCOBookImb1"==DATA$names_nodes)] =0 
  
  
  is.DAG(g)
  save.image(file=LocalSavedResults)
}else{
  # load the environment to save time for graphical analysis
  load(LocalSavedResults)
  print("PLEASE RELOAD THE FUNCTIONS IN THE HEADER OF THIS FILE")
  source("CausalTimeSeries_fx.R")
  }

# plot the graph




# We further explored the causal relationships by employing a grid search with thinner penalty 
# parameters for cross-validation. Interestingly, the resulting causal graph remained largely unchanged, 
# with only a few variations observed (5 edges out of 169). However, to ensure a more conservative approach 
# and mitigate any potential impact of noise on the accurate identification of the causal graph, 
# we opted to use a more conservative penalty. 
# This decision helps minimize the risk of noise interfering with the correct identification of causal relations.

# nlambda = 20
# lambda_min_ratio = 0.01
# boundary warning
# DAG
# 169

# nlambda = 20*2
# lambda_min_ratio = 0.01/2
#  SPYBookImb10 -> CSCOBookImb1 which generates many cycles
# [1] "DISBookImb4" "DISBookImb6" "DISBookImb5" "DISBookImb4"
# [1] "KOBookImb1" "KOBookImb2" "KOBookImb3" "KOBookImb1"
# [1] "SPYBookImb1" "SPYBookImb3" "SPYBookImb2" "SPYBookImb1"
# 179

# nlambda = 20*10 
# lambda_min_ratio = 0.01/10
#[1] "DISBookImb4" "DISBookImb6" "DISBookImb5" "DISBookImb4"
# [1] "KOBookImb5"  "KOBookImb9"  "KOBookImb10" "KOBookImb5" 
# [1] "SPYBookImb1" "SPYBookImb3" "SPYBookImb2" "SPYBookImb1"
# 176 edges

aux = pc_res$adj_mat
row.names(aux) = DATA$names_nodes
colnames(aux) = DATA$names_nodes 

# plot the graphs. A lot of Graphical options
{
  
  
title_main = sprintf("")

name_file = sprintf("Lasso-PC_%s_ratio_1.pdf","HFT")
pdf(name_file, width = 60, height = 40, pointsize = 30) 
{
# old graphical option
# { g1 <- graph_from_adjacency_matrix( aux)
#   
# 
#   
#     y<-c( 
#       rep(350+c(150,250),K/10),#first level
#       rep(350+c(150,250),K/10),
#       rep(200+c(150,250),K/10),
#       rep(200+c(150,250),K/10),
#       rep(50+c(150,250),K/10))
#     y[seq(1,12,by = 4)] = 175+350
#     y[seq(4,12,by = 4)] = 225+350 
#     
#     y[seq(1,12,by = 4)+12] = 175+350
#     y[seq(4,12,by = 4)+12] = 225+350 
#     
#     y[seq(1,12,by = 4)+12*2] = 175+200
#     y[seq(4,12,by = 4)+12*2] = 225+200 
#     
#     y[seq(1,12,by = 4)+12*3] = 175+200
#     y[seq(4,12,by = 4)+12*3] = 225+200
#     # -1*c(rep(c(100,175),11),100))
#     x<- c(linspace(0, 1000, n = 2*K/5),
#           linspace(0, 1000, n = 2*K/5),
#           linspace(0, 1000, n = K/5))
#     l<-data.frame(x,y)
#     l<-as.matrix(l)
#  
#   
#   E(g1)$color <- "royalblue" 
#   
#   plot(g1,layout=l, edge.arrow.size=0.8, vertex.label.color="black", 
#        vertex.size=800,
#        vertex.label.dist=225,vertex.label.cex=1.6,edge.curved=0.1,
#        ylim=c(100,700),xlim=c(0,1150),axes=FALSE,asp=0, rescale = F,
#        main = "")}
}
{ g1 <- graph_from_adjacency_matrix( aux)
  
  
  
  y<-c( 
    rep(350+c(150,250),K/10),#first level
    rep(350+c(150,250),K/10),
    rep(200+c(150,250),K/10),
    rep(200+c(150,250),K/10),
    rep(50+c(150,250),K/10))
  
  y[1] = 1000
  y[c(2,3)] = 900
  y[c(4,5)] = 800
  y[c(6,7)] = 700
  y[c(8,9)] = 600
  y[c(10,11)] = 500
  y[c(12)] = 400
  
  y[1+12] = 1000
  y[c(2,3)+12] = 900
  y[c(4,5)+12] = 800
  y[c(6,7)+12] = 700
  y[c(8,9)+12] = 600
  y[c(10,11)+12] = 500
  y[c(12)+12] = 400
  
  y[1+12*2] = 1000
  y[c(2,3)+12*2] = 900
  y[c(4,5)+12*2] = 800
  y[c(6,7)+12*2] = 700
  y[c(8,9)+12*2] = 600
  y[c(10,11)+12*2] = 500
  y[c(12)+12*2] = 400
  
  y[1+12*3] = 300
  y[c(2,3)+12*3] = 200
  y[c(4,5)+12*3] = 100
  y[c(6,7)+12*3] = 0
  y[c(8,9)+12*3] = -100
  y[c(10,11)+12*3] = -200
  y[c(12)+12*3] = -300
  
  y[1+12*4] = 300
  y[c(2,3)+12*4] = 200
  y[c(4,5)+12*4] = 100
  y[c(6,7)+12*4] = 0
  y[c(8,9)+12*4] = -100
  y[c(10,11)+12*4] = -200
  y[c(12)+12*4] = -300
  
 
  # -1*c(rep(c(100,175),11),100))
  x<- c(linspace(0, 1000, n = 2*K/5),
        linspace(0, 1000, n = 2*K/5),
        linspace(0, 1000, n = K/5))
  
  x[seq(2,11,by=2)] = 0
  x[seq(3,11,by=2)] = 200
  x[c(1,12)] = 100
  
  x[seq(2,11,by=2)+12] = 0+300
  x[seq(3,11,by=2)+12] = 200+300
  x[c(1,12)+12] = 100+300
  
  x[seq(2,11,by=2)+12*2] = 0+300*2
  x[seq(3,11,by=2)+12*2] = 200+300*2
  x[c(1,12)+12*2] = 100+300*2
  
  x[seq(2,11,by=2)+12*4] = 0+200-50
  x[seq(3,11,by=2)+12*4] = 200+200-50
  x[c(1,12)+12*4] = 100+200-50
  
  x[seq(2,11,by=2)+12*3] = 0+200+300+50
  x[seq(3,11,by=2)+12*3] = 200+200+300+50
  x[c(1,12)+12*3] = 100+200+300+50

  
  l<-data.frame(x,y)
  l<-as.matrix(l)
  
  
  E(g1)$color <- "royalblue" 
}  
  plot(g1,layout=l, edge.arrow.size=0.8, vertex.label.color="black", 
       vertex.size=800,
       vertex.label.dist=500,vertex.label.cex=1.6,edge.curved=0.2,
       ylim=c(-350,1050),xlim=c(0,850),axes=FALSE,asp=0, rescale = F,
       main = "")

dev.off() 


# plot only the book of 1 instrument CSCO
 

# select CSCO
sel_var = "CSCO"
iCSCO = which(grepl(sel_var,DATA$names_nodes, )== TRUE)
aux_csco = aux[iCSCO, iCSCO]

title_main = sprintf("")

name_file = sprintf("Lasso-PC_%s_ratio_1.pdf","CSCO")
pdf(name_file, width = 60, height = 40, pointsize = 30) 
 
{ g_csco <- graph_from_adjacency_matrix( aux_csco)
  
  
  
  y<-c( 
    rep(350+c(150,250),K/10) )
  
  y[1] = 1000
  y[c(2,3)] = 900
  y[c(4,5)] = 800
  y[c(6,7)] = 700
  y[c(8,9)] = 600
  y[c(10,11)] = 500
  y[c(12)] = 400
  
   
  # -1*c(rep(c(100,175),11),100))
  x<- c(1:12)
  
  x[seq(2,11,by=2)] = 0
  x[seq(3,11,by=2)] = 200
  x[c(1,12)] = 100
  
   
  l<-data.frame(x,y)
  l<-as.matrix(l)
  
  
  E(g_csco)$color <- "royalblue" 
}  

plot(g_csco,layout=l, edge.arrow.size=0.8, vertex.label.color="black", 
     vertex.size=400,
     vertex.label.dist=250,vertex.label.cex=1.6,edge.curved=0.1,
     ylim=c(350,1050),xlim=c(-50,250),axes=FALSE,asp=0, rescale = F,
     main = "")

dev.off()

# select AMZN CSCO SPY
# we plot the structure of AMZN and the child of AMZN, i.e., where the information from a change in the shape 
# of the book of AMZN will flow in which other instrument
sel_var = "AMZN"
iAMZN= which(grepl(sel_var,DATA$names_nodes, )== TRUE)
aux_AMZN = aux[iAMZN, iAMZN]

# select only those nodes which are child of AMZN BOOK-RET-TRADE IMB
sum(aux[iAMZN,-iAMZN])
aux_CHILD_AMZN = aux[iAMZN,-iAMZN]
i_to_select = NULL
for (ii in 1:ncol(aux_CHILD_AMZN)){
  if( any(aux_CHILD_AMZN[,ii]==1) ){
    i_to_select = c(i_to_select, ii)
  }
}
iAMZN_child = c(DATA$names_nodes[iAMZN], colnames(aux_CHILD_AMZN)[i_to_select])
aux_amzn_child = aux[iAMZN_child,iAMZN_child]
# for graphical reason we do not plot the relation among the children of AMZN
aux_amzn_child[ colnames(aux_CHILD_AMZN)[i_to_select], colnames(aux_CHILD_AMZN)[i_to_select]] = 0

title_main = sprintf("")

name_file = sprintf("Lasso-PC_%s_ratio_1.pdf","AMZN_and_CHILDREN")
pdf(name_file, width = 60, height = 40, pointsize = 30) 

{ g_amzn_child <- graph_from_adjacency_matrix( aux_amzn_child)
  
  
  
  y<-c( 
    rep(350+c(150,250),K/10),
    1:8)
  
  y[1] = 1000
  y[c(2,3)] = 900
  y[c(4,5)] = 800
  y[c(6,7)] = 700
  y[c(8,9)] = 600
  y[c(10,11)] = 500
  y[c(12)] = 400
  
  y[c(13:16)] = 200
  
  y[c(17:20)] = 100
  
  # -1*c(rep(c(100,175),11),100))
  x<- c(1:12)
  
  x[seq(2,11,by=2)] = 50
  x[seq(3,11,by=2)] = 150
  x[c(1,12)] = 100
  
  x[c(13,17)] = 25
  x[c(14,18)] = 75

  x[c(16,20)] = 125
  
  x[c(14,19)] = 175
  x[c(15)] = 75
  
  x[c(17:20)]  =   x[c(17:20)] +25
   
  l<-data.frame(x,y)
  l<-as.matrix(l)
  
  
  E(g_amzn_child)$color <- "royalblue" 
}  

plot(g_amzn_child,layout=l, edge.arrow.size=1.2, vertex.label.color="black", 
     vertex.size=300,
     vertex.label.dist=250,vertex.label.cex=1.6,edge.curved=0*c(rep(0.1,21), rep(0,8)),
     ylim=c(50,1050),xlim=c(0,225),axes=FALSE,asp=0, rescale = F,
     main = "")

dev.off()



# select CSCO with parents and Children
#  
sel_var = "CSCO"
iCSCO= which(grepl(sel_var,DATA$names_nodes, )== TRUE)
aux_CSCO = aux[iCSCO, iCSCO]

# select only those nodes which are child of AMZN BOOK-RET-TRADE IMB
sum(aux[iCSCO,-iCSCO])
aux_CHILD_CSCO = aux[iCSCO,-iCSCO]
i_to_select = NULL
for (ii in 1:ncol(aux_CHILD_CSCO)){
  if( any(aux_CHILD_CSCO[,ii]==1) ){
    i_to_select = c(i_to_select, ii)
  }
}
iCSCO_child = c(DATA$names_nodes[iCSCO], colnames(aux_CHILD_CSCO)[i_to_select])
aux_csco_child = aux[iCSCO_child,iCSCO_child]
# for graphical reason we do not plot the relation among the children of AMZN
aux_csco_child[ colnames(aux_CHILD_CSCO)[i_to_select], colnames(aux_CHILD_CSCO)[i_to_select]] = 0


# select the parents of CSCO- there are only two parents
sum(aux[-iCSCO,iCSCO])
aux_PA_CSCO = aux[-iCSCO,iCSCO]
i_to_select_pa = NULL
for (ii in 1:nrow(aux_PA_CSCO)){
  if( any(aux_PA_CSCO[ii,]==1) ){
    i_to_select_pa = c(i_to_select_pa, ii)
  }
}
row.names(aux_PA_CSCO)[i_to_select_pa] #AMZN Ret and AMZN TradeImb

iCSCO_child_pa = c(row.names(aux_PA_CSCO)[i_to_select_pa], DATA$names_nodes[iCSCO], colnames(aux_CHILD_CSCO)[i_to_select])
aux_csco_child_pa = aux[iCSCO_child_pa,iCSCO_child_pa]
aux_csco_child_pa[ colnames(aux_CHILD_CSCO)[i_to_select], colnames(aux_CHILD_CSCO)[i_to_select]] = 0
aux_csco_child_pa[row.names(aux_PA_CSCO)[i_to_select_pa],colnames(aux_CHILD_CSCO)[i_to_select]]  = 0

g_csco_child_pa <- graph_from_adjacency_matrix( aux_csco_child_pa)

title_main = sprintf("")

name_file = sprintf("Lasso-PC_%s_ratio_1.pdf","CSCO_PARENTS_and_CHILDREN")
pdf(name_file, width = 60, height = 40, pointsize = 30) 
# tiff("test.tiff", units="cm", width=60, height=40, res=300)


{ g_csco_child_pa <- graph_from_adjacency_matrix( aux_csco_child_pa)
  
  
  
  y<-c(1:2, 
    rep(350+c(150,250),K/10),
    1:9)
  
  y[1:2] = 1100
  y[3] = 1000
  y[2+c(2,3)] = 900
  y[2+c(4,5)] = 800
  y[2+c(6,7)] = 700
  y[2+c(8,9)] = 600
  y[2+c(10,11)] = 500
  y[2+c(12)] = 400
  y[13] = 400
  
  y[2+c(13:16)] = 200
  
  y[2+c(17:21)] = 100
  
  y[16] = 100
  y[23] = 200
  y[19] = 200

  # -1*c(rep(c(100,175),11),100))
  x<- c(1:12)
  x[1] = 125
  x[2] = 75
  x[seq(4,12,by=2)] = 50
  x[seq(5,13,by=2)] = 150
  x[c(3,12+2)] = 100
  x[13] = 200
  
  x[15] = 25
  x[17] = 75
  x[18] = 125
  x[20] = 50
  x[21] = 100
  
  x[16] = 150
  x[23] = 175
  
  x[19] = 250
  x[22] = 225
  l<-data.frame(x,y)
  l<-as.matrix(l)
  
  
  E(g_csco_child_pa)$color <- "royalblue" 
}  

plot(g_csco_child_pa,layout=l, edge.arrow.size=1.2, vertex.label.color="black", 
     vertex.size=300,
     vertex.label.dist=250,vertex.label.cex=1.6,edge.curved=0,
     ylim=c(50,1150),xlim=c(0,275),axes=FALSE,asp=0, rescale = F,
     main = "")

dev.off()


}

# interesting nodes. Plot the parents of the selected node
sel_var = "SPYRet"
ix_var = which(DATA$names_nodes == sel_var)
DATA$names_nodes[which(g[,ix_var]==1)]

sel_var = "CSCORet"
ix_var = which(DATA$names_nodes == sel_var)
DATA$names_nodes[which(g[,ix_var]==1)]

sel_var = "CSCOBookImb2"
ix_var = which(DATA$names_nodes == sel_var)
DATA$names_nodes[which(g[,ix_var]==1)]

sel_var = "CSCOBookImb1"
ix_var = which(DATA$names_nodes == sel_var)
DATA$names_nodes[which(g[,ix_var]==1)]


# sel_var = "SPYBookImb1"
# ix_var = which(DATA$names_nodes == sel_var)
# DATA$names_nodes[which(g[,ix_var]==1)]
# 
# 
# sel_var = "SPYTradeImb"
# ix_var = which(DATA$names_nodes == sel_var)
# DATA$names_nodes[which(g[,ix_var]==1)]
 
#----------------------  

#----------------------------------------------------------------------------------#
#                                                                   
# GET IRFS                              
#----------------------------------------------------------------------------------#
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
  
  find_directe_cycles = function(g) {
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
  
  
  cycles = find_directe_cycles(g)
  
  
  if (length(cycles)>=1){
    for (i in cycles){
      print(col_names_X[i])
    }
  }
}else{
  max_horizon = 20
  # number of sampling to compute the IRFs
  m = 500
  delta_shock = 1
  # SELECT THE VARIABLE TO SHOCK
  alpha_irf_ci = 0.1# confidence interval for irfs of size 1-alpha_irf_ci, 10% as in Kanzig 
  
  flag_cumsum = FALSE#compute the cumulative irfs
  k_shock_variable = which(col_names_X %in%c("CSCOBookImb1","CSCOBookImb2","CSCORet","SPYRet"))   #
  # where to look the response
  k_shocked_variable = which(col_names_X %in%c("CSCORet","SPYRet"))   #
  
  # convert datetime object to days
  date <- as.Date(DATA$rowX, format = "%Y-%m-%d")
  day_number <- as.numeric(date)
  # this object will be used in the block bootstrap, i.e., we will bootstrap on the day,
  # preserving the intraday pattern
  
  load_from_local = TRUE # Set FALSE if you want to run (it is quite slow)
  LocalFileIrfs ="irf_hft/IRFs_HFT.RData" 
  # YOU CAN CHOOSE THE GRAPH YOU WANT BY SELECTING THE RELATED Adjacency matrix
  adj_mat     = pc_res$adj_mat
  #find the SVAR parameters
  svar_params = get_SVAR_params(Sigma_eps_hat, adj_mat)
  PI          = svar_params$PI
  D           = svar_params$D
  H           = svar_params$H
  # ------------------------------------------------------------ get the IRFs 
  # point estimate
  f_k_inverse = get_f_k_inverse(X)
  
  
  
  # IRFs_X = get_IRF_X(m = m, adj_mat, Sigma_eps_hat,
  #                  Theta_hat, f_k_inverse, max_horizon = max_horizon,
  #                  k_shock_variable = k_shock_variable, delta_shock = delta_shock, seed_num=1)
  
  # bootstrap for confidence interval 
  svar_params_comp = get_svar_params_companion(Sigma_eps_hat, Theta_hat, H, PI,  K)
  A_comp = svar_params_comp$A
  
  #===================================
  # do not do the cumulative
  IRFs_X_point_est_all = NULL
  for (ii_kk in 1:length(k_shock_variable)){
    
    i_k_shock_variable = k_shock_variable[ii_kk]
    IRFs_X_point_est_all[[ii_kk]] <- get_IRF_X(m = m, adj_mat, Sigma_eps_hat,
                                                      Theta_hat, f_k_inverse, max_horizon = max_horizon,
                                                      k_shock_variable = i_k_shock_variable, 
                                                      delta_shock = delta_shock, seed_num=0,
                                                      flag_shock_Z=FALSE, flag_cum = FALSE)
  }
  # reformat point est for get_CI_bootCentered
  # restrict the IRF to only the k_shocked_variable to be consistent with IRFs_X_boot
  IRFs_X_point_est = matrix(NA, 
                            nrow=max_horizon+1,
                            ncol=length(k_shock_variable)*length(k_shocked_variable))
  n_var_schocked = length(k_shocked_variable)
  for (ii_kk in 1:length(k_shock_variable)){
    selected_columns <-  ((  (ii_kk-1)*n_var_schocked+1):( ii_kk*n_var_schocked))
    IRFs_X_point_est[,c(selected_columns)] = IRFs_X_point_est_all[[ii_kk]][, k_shocked_variable]
  }
  
  
  #===================================
  
  
  if (!load_from_local){
    tictoc::tic()
    # IRFs_X_boot = get_IRFs_bootstrap_selected_shocks(mboot = 200, n, m=m, Sigma_hat, Sigma_eps_hat, A_comp,
    #                                  f_k_inverse, adj_mat,
    #                                  p_lag, K, max_horizon,
    #                                  k_shock_variable = k_shock_variable,
    #                                  lambda_cv_adj, tau_cv_adj, seed = 0,
    #                                  k_shocked_variable = k_shocked_variable,
    #                                  flag_block = FALSE,
    #                                  X, day_number )
    # tictoc::toc()
    flag_shock_Z = FALSE # compute irf on Z or on X
    IRFs_X_boot = get_IRFs_bootstrap_selected_shocks(mboot = 200, n, m=m, 
                                                     Sigma_hat, Sigma_eps_hat, A_comp,
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
    save(IRFs_X_boot,file=LocalFileIrfs)
  }else{
    load(LocalFileIrfs)
  }
  
  
  IRFs_X_boot_CI = get_CI_bootCentered(IRFs_X_boot,  alpha_irf_ci, IRFs_X_point_est)
  
  IRF_lwr_raw = t(IRFs_X_boot_CI$IRFs_X_boot_lwr)
  IRF_upr_raw = t(IRFs_X_boot_CI$IRFs_X_boot_upr)
  IRF_avg_raw = t(IRFs_X_boot_CI$IRFs_X_boot_mean)
  IRF_avg_raw = t(IRFs_X_point_est)
  
  
  # #scale
  IRF_lwr = IRF_lwr_raw 
  IRF_upr = IRF_upr_raw 
  IRF_avg = IRF_avg_raw 

  # IRFs_X_boot_CI = get_CI_boot(IRFs_X_boot, alpha = 0.01)
  
  # takes the cumulative IRFs since in thi example we have diff the variables
  # K_IRFs = length(k_shocked_variable)*length(k_shock_variable)
  # IRFs_X = matrix(0,nrow = max_horizon+1,ncol = K_IRFs)
  # plot not cumulative
  plot_ci = FALSE
  par(mar = c(2, 2, 2, 2))
  kk = 0
  for (ikshock in 1:length(k_shock_variable)){
    for (ikschoked in 1:length(k_shocked_variable)){
      kk = kk+1
      if (plot_ci){
      png(sprintf("irf_hft/IRFs_%s_to_%s.png",DATA$names_nodes[k_shock_variable[ikshock]],
                  DATA$names_nodes[k_shocked_variable[ikschoked]]), 
          width = 800, height = 400)
      }else{
        png(sprintf("irf_hft/point_estimate_IRFs_%s_to_%s.png",DATA$names_nodes[k_shock_variable[ikshock]],
                    DATA$names_nodes[k_shocked_variable[ikschoked]]), 
            width = 800, height = 400)
      }
      plot(IRF_avg[kk,],type="l", cex.axis=2,cex.main=2,
           # NOME MAIN con ->
           main =sprintf("%s -> %s",DATA$names_nodes[k_shock_variable[ikshock]],
                         DATA$names_nodes[k_shocked_variable[ikschoked]]),
           xlab="",ylab="", col = "black",lty=1,lwd=3,
           ylim = c(min(IRF_lwr[kk,]), max(IRF_upr[kk,])))
      # lines(IRFs_X_boot_CI$IRFs_X_boot_lwr[,kk], type="l",col="red", lty=3)
      # lines(IRFs_X_boot_CI$IRFs_X_boot_upr[,kk], type="l",col="red", lty=3)
      if (plot_ci){
      lines(IRF_lwr[kk,], type="l",col='red', lty=3)
      lines(IRF_upr[kk,], type="l",col='red', lty=3)
      }
      abline(h=0)
      dev.off()
    }
  }
  
  
}
