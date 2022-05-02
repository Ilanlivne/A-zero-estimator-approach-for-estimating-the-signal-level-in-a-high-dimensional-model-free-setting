


  rm(list = ls())
  packages <- c("future", "lubridate","furrr", "ggplot2","tidyr","dplyr","latex2exp","gridExtra","Rcpp","scales","stargazer","cccp")
  
  start_time <- Sys.time()
  lapply(packages, require, character.only = TRUE)
  source("run_sim.R")
  Rcpp::sourceCpp("zero_terms_function.cpp")
  #options(future.globals.maxSize= 2672864000)
  
  
  ########################################
  #### initial parameters settings
  ######################################
  
  n_cores =72  
  n_simulation = 20
  n_pop<- 1*10^6
  n = p = 300
  c = 6  #number of large betas. sparse/dense parameter.
  tau2 = 2 
  sig_eps = 1  
  for (prop_large in c( 0.1, 0.3, 0.5, 0.7, 0.9)) {
  
  
  
################################################################
###### semi-supervised data 
##########################################################
  
  X_pop<-matrix(
     rexp(n_pop *p,1)-1
     ,nrow=n_pop,ncol=p)
  
  lamda_m2 = 1.42071^2   #[1+E(Xsin(X))]^2
  tau2_B=prop_large*tau2
  large = sqrt(  tau2_B/(lamda_m2*c) )
  small = sqrt( (tau2-tau2_B)/( (p-c)*lamda_m2 ) )  
  delta<-c(rep(large, c), rep(small,p-c))
  
  
  ############################################
  ### covariate selection algorithm 
  ##############################################
  selection_algoritm <- function(X,Y){
    # step 1: calculate \hat\beta_j^2 for j=1,...,p:  
    W_3 = as.data.frame(X)*Y
    mean_squared_W<- colMeans(W_3^2)  #calculate first element of beta_square_hat
    var_W<-sapply(W_3,var) # #calculate the second element of beta
    beta_square_hat<- mean_squared_W- var_W # calculate vector of beta2_hat
    # step2:  calculate the differences lamda_j for j=2,...,p:  
    dt <- data.frame(
      delta = delta, delta_type = c(rep("big",c),rep("small",p-c)),
      j = 1:p,  beta_square_hat   ) %>% arrange(beta_square_hat) %>%
      mutate( index = 1:n(),  lag_1 = lag(beta_square_hat), 
              diff = beta_square_hat - lag_1
      )%>%dplyr::select(-lag_1) %>% filter(diff != "NA")  
    #step 3: Select the covariates of S_gamma.
    #calculate  j_star:
    index_star = dt %>%  mutate(max_diff = max(diff, na.rm = T))%>% 
      filter(diff == max_diff) %>% dplyr::select(index) %>% unlist()
    if (p-index_star<2) {
      index_star = index_star  %>% unlist()-1
    }
    dt<- dt %>% mutate(pred =if_else(index >= index_star,"big","small") )
    return(dt)
  }  
  
  
  Y_pop =  X_pop%*%delta+  sin(X_pop)%*%delta + rexp(n_pop,1)-1
    
  ##############################################
  #Best linear predictor (BLP) parameters
  BLP<-(t(X_pop)%*%Y_pop)/n_pop
  (alpha <- mean(Y_pop))
  (sigma2<- var(as.vector(Y_pop))-tau2)   # check for sigma2 = 1 
  
  (mehane <- var(double_dist_sum(X_pop)))
  sim_fn <- function(setting_param_list){
  Rcpp::sourceCpp("zero_terms_function.cpp")
  seed <-  setting_param_list$seed
  set.seed(seed)
  p = n

##################################################################
    ##  sample data  ###
 X<-matrix(
  rexp(n *p,1)-1,nrow=n,ncol=p)
  Y = X%*%delta + sin(X)%*%delta +  rnorm(n,0,sig_eps)    
  
##############################################################################################333 
###                 naive estimator
##################################################################################
  W = as.data.frame(X)*Y
  mean_squared_W<- colMeans(W^2)  #calculate first element of beta_square_hat
  var_W<-sapply(W,var) # #calculate the second element of beta
  beta_square_hat<- mean_squared_W - var_W # calculate vector of beta2_hat
  tau2_hat <- sum(beta_square_hat)
  naive<-if_else(tau2_hat<0,0,tau2_hat)
  

  dt<- selection_algoritm(X,Y)   
  estimated_indexes <- filter(dt,pred == "big") %>% dplyr::select(j)
  mehane_selection <- var(double_dist_sum(X_pop[,estimated_indexes$j]))
    
  
  #####################################################
  #####       Estimator T_{\hat g}    ##########      
  X<- as.matrix(X)
  g <- mean(double_dist_sum(X))
  mone_s<-func22(X,Y)*4/(n*(n-1))
  (a_star_hat <- mone_s /mehane)
  single_coeff_estimator <-  naive - a_star_hat*g
  ####################################################################
  #####         Estimator T_{\hat h}
  #######################################################################
  
  X_selection <- as.matrix(X[,1:c])
  g_selection <- mean(double_dist_sum(X_selection))
  mone_selection<-func22(X_selection,Y)*4/(n*(n-1))
  a_star_hat_selection <- mone_selection /mehane_selection
  single_coeff_estimator_selection <-  naive - a_star_hat_selection*g_selection
  ####################################################################3
  
  #  create final table
  Naive <- naive
  Single <- if_else(single_coeff_estimator<0,0,single_coeff_estimator)
  Single_selection <- if_else(single_coeff_estimator_selection<0,0,single_coeff_estimator_selection)
  
  
  return(tibble( Naive = naive
                 ,T_1 = Single
                 ,T_2 = Single_selection
                 
  ))
}


  #######################################################################################3
  
  ss=451  #seed
  setting_param_list <-  list(seed = ss:(n_simulation+ss))
  result <-run_sim(sim_fn, setting_param_list, nworkers = n_cores)
  raw_storage<-gather(result, method, value, 3:ncol(result)) 
  end_time <- Sys.time()
  total_time  <- round(end_time - start_time)
  
  
  
  
  ###############################################
  #    Summary Table
  ##############################################
  
  
  stat<- raw_storage %>% group_by(method) %>%
    summarise(Mean =round(mean(value),2)
              ,SD = round(sd(value),3)
              ,MSE=mean((value-tau2)^2) 
              ,SD_RMSE = round(sqrt(var((value-tau2)^2)/(4*MSE*n_simulation)),3)
    ) %>% mutate(RMSE = round(sqrt(MSE),3), bias = round(tau2-Mean,2)
                 ,tau2=round(tau2,1), n=n , prop_large = prop_large   )
  
  RMSE_Naive<- stat %>% filter(method == "Naive") %>%  dplyr::select(RMSE)%>%pull()
  
  
  stat <- stat %>%   dplyr::select(prop_large,Estimator = method,bias
                                   ,SD
                                   ,RMSE
                                   ,SD_RMSE
                                   
  )%>%  mutate(tau2 = round(tau2,1), sigma2 = round(sigma2,2),
               n_sim = n_simulation
               , "# large betas" =c
               , prop_large =  label_percent()(prop_large)
               , n = n, n_pop = n_pop , total_time=total_time
               , pct_imp = round((RMSE  /RMSE_Naive - 1) * 100,2) 
               ,ss
  )
  
  mytime <- format(Sys.time(), "%dday_%Hhour_%Mminute_%Sseconds")
  
  
  
  
  #####################################################
  ###  outputs
  #####################################################
  
  latex_table <- stat%>%  mutate(p=n)%>%
    dplyr::select(prop_large,  tau2, n=n, Estimator, Bias=bias, SE=SD, RMSE
                  ,"% Change" = pct_imp
                  , SE.RMSE = SD_RMSE
    ) %>%  
    arrange(prop_large, factor(Estimator, levels = Estimator[c(1,2,3)]))
  
  sink(file = paste0(mytime,".tex"))
  stargazer(latex_table,summary = FALSE,rownames = FALSE, digits = 2)
  sink(file = NULL)
  
  
}


stat