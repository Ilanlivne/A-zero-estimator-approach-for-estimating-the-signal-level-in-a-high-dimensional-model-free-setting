rm(list = ls())
packages <- c("ggplot2","tidyr","dplyr","latex2exp","gridExtra" ,"Rcpp","scales","stargazer","cccp")
lapply(packages, require, character.only = TRUE)
source("run_sim_multicore.R")
Rcpp::sourceCpp("zero_terms_functions.cpp")
start_time <- Sys.time()


########################################
#### initial parameters settings
######################################

n_cores =72  
n_simulation = 100
n_pop<- 1*10^6
n = p = 300
c =6  #number of large betas. sparse/dense parameter.
tau2 = 1
sig_eps = 1
for (prop_large in c( 0.3, 0.5, 0.7, 0.9)) {
B=100


################################################################
###### semi-supervised data 
##########################################################
lamda_m2 = 1.42071^2   #[1+E(Xsin(X))]^2


tau2_B=prop_large*tau2
large = sqrt(  tau2_B/(lamda_m2*c) )
small = sqrt( (tau2-tau2_B)/( (p-c)*lamda_m2 ) )  
delta<-c(rep(large, c), rep(small,p-c))

X_pop<-matrix(rexp(n_pop *p,1)-1 ,nrow=n_pop,ncol=p)
Y_pop =  X_pop%*%delta+  sin(X_pop)%*%delta +  rnorm(n_pop,0,sig_eps)
  
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
            #         diff = beta_square_hat - lag_1
            diff = beta_square_hat-lag_1
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



##################
#BLP parameters
BLP<-(t(X_pop)%*%Y_pop)/n_pop
(alpha <- mean(Y_pop))
(sigma2<- var(as.vector(Y_pop))-tau2)   # check for sigma2 = 1 



(mehane <- var(double_dist_sum(X_pop))/n)  # var of g_n
#(mehane_selection <- var(double_dist_sum(X_pop[,1:c]))/n)




sim_fn <- function(setting_param_list){
  Rcpp::sourceCpp("zero_terms_functions.cpp")
  
  EigenPrism <- function(y,X,invsqrtSig=NULL,alpha=0.05,target='beta2',zero.ind=c(),diagnostics=T){
    # Author: Lucas Janson (statweb.stanford.edu/~ljanson)
    # Runs EigenPrism procedure for estimating and generating confidence
    #  intervals for variance components in high-dimensional linear model:
    #       y = X%*%beta + e,   rows of X iid~ N(0,Sig),   e iid~ N(0,sigma^2)
    #  Requires cccp package for solving second order cone optimization.
    #  Note confidence interval endpoints may lie outside parameter domain, so it may be appropriate
    #   to clip them after the fact.
    # 
    # Inputs:
    #  y: response vector of length n (will automatically be centered)
    #  X: n by p design matrix; columns will automatically be centered and scaled to variance 1;
    #      should not contain intercept column, since both y and X will be centered
    #  invsqrtSig: if columns of X not independent, p by p positive definite matrix which is the square-root
    #               of the inverse of Sig, where Sig is the *correlation* matrix of the X (default is identity)
    #  alpha: significance level for confidence interval (default = 0.05)
    #  target: target of estimation/inference
    #		  'beta2' (default) is the squared 2-norm of the coefficient vector: sum(beta^2)
    #           'sigma2' is the noise variance sigma^2
    #           'heritability' is the fraction of variance of y explained by X%*%beta: t(beta)%*%Sig%*%beta/var(y)
    #  zero.ind: vector of which indices of the weight vector w to constrain to zero (default is none)
    #  diagnostics: boolean (default = T) for whether to generate diagnostic plots for the V_i, lambda_i, and w_i
    #  
    # Outputs:
    #  estimate: unbiased estimate of the target (for heritability, only approximately unbiased)
    #  CI: 100*(1-alpha)% confidence interval for target
    
    # Get dimensionality of problem
    n = nrow(X)
    p = ncol(X)
    
    # Transform y and X to proper form
    y = y-mean(y)
    X = scale(X,T,T)*n/(n-1)
    if(!is.null(invsqrtSig)) X = X%*%invsqrtSig
    
    # Take singular value decomposition and rescale singular values
    svd = svd(X)
    lambda = svd$d^2/p
    
    # Defined cone-constrained linear problem to optimize weights; [v; w] is vector of optimization variables
    q = c(1,rep(0,n)) #coefficient vector in objective function
    A = rbind(c(0,rep(1,n)),c(0,lambda)) #matrix for linear constraints
    b = c(0,1) #vector for linear constraints
    if(target=='sigma2') b = c(1,0) #switch constraints if target is sigma^2
    # Constrain some weights to be zero if desired
    if(!is.null(zero.ind)){
      A = rbind(A,cbind(rep(0,length(zero.ind)),diag(rep(1,n))[zero.ind,]))
      b = c(b,rep(0,length(zero.ind)))
    }
    # Define second-order cone constraints
    soc1 = socc(diag(c(1/4,rep(1,n))),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
    soc2 = socc(diag(c(1/4,lambda)),c(-1/2,rep(0,n)),c(1/4,rep(0,n)),1/2)
    prob = dlp(as.vector(q),A,as.vector(b),list(soc1,soc2))
    
    # Solve optimization problem and extract variables
    opt = cps(prob,ctrl(trace=F))
    v = getx(opt)[1]
    w = getx(opt)[-1]
    
    # Compute estimate and y's variance
    est = sum(w*(t(svd$u)%*%y)^2)
    yvar = sum(y^2)/n
    
    # Compute confidence interval
    CI = est + yvar*sqrt(v)*qnorm(1-alpha/2)*c(-1,1)
    if(target=='heritability'){
      est = est/yvar
      CI = CI/yvar
    }
    
    # Generate list with results
    result=list()
    result$estimate = est
    result$CI = CI
    
    # Generate diagnostic plots
    if(diagnostics){
      par(mfrow=c(1,3))
      
      # Check that eigenvectors are approximately Gaussian
      nV = floor(log10(n))
      srtV = svd$v[,10^(0:nV)]
      labs = c()
      for(i in 1:(nV+1)){
        srtV[,i] = sort(srtV[,i])
        ind = 10^(i-1)
        labs = c(labs,bquote(V[.(ind)]))
      }
      matplot(qnorm((1:p)/(p+1)),srtV,type="l",lwd=2,
              ylab="Quantiles of Eigenvectors",xlab="Gaussian Quantiles",
              main=expression(paste("Check Gaussianity of Eigenvectors ",V[i])))
      legend("topleft",as.expression(labs),col=1:(nV+1),lty=1:(nV+1),lwd=2)
      
      # Check that there are no outliers in the eigenvalues
      hist(lambda,main=expression(paste("Histogram of Normalized Eigenvalues ",lambda[i])),
           xlab=expression(lambda[i]))
      
      # Check that the weights are not dominated by just a few values
      srtw = sort(abs(w),T)
      plot(1:n,cumsum(srtw)/sum(srtw),type="l",lwd=2,
           main=expression(paste("Fraction of Total Weight in Largest k ",w[i])),
           xlab="k",ylab="Fraction of Total Weight")
    }
    
    return(result)
  }
  
  
  seed <-  setting_param_list$seed
  set.seed(seed)
  n <- setting_param_list$n
  p = n

  ##  sample data  ###
  X <- matrix(rexp(n * p,1)-1, nrow=n, ncol=p)
  Y <- X%*%delta + sin(X)%*%delta + rnorm(n,0,sig_eps)    
  
  
  Eigen_tau2 <- EigenPrism(y=Y,X=X,target='beta2',diagnostics=F)$estimate  
  
  
  #####        Single coefficient mixed with selection estimator     
  
  dt<- selection_algoritm(X,Y)   
  estimated_indexes <- filter(dt,pred == "big") %>% dplyr::select(j)
  
  X_selection <- as.matrix(X[,estimated_indexes$j])
  g_selection <- mean(double_dist_sum(X_selection))
  
  mehane_selection <- var(double_dist_sum(X_pop[,estimated_indexes$j]))
  
  
  
  X<- as.matrix(X)
  g_single <- mean(double_dist_sum(X))
  
  ###   Bootstraping step 
  Eigen_Pri_b = g_single_b  = g_single_selection_b = NA  
  for (b in 1:B) {    
    # step 3.1: Resample with replacement
    indices <- sample(1:nrow(X),size = nrow(X) ,replace=TRUE) 
    X_b <- X[indices,]
    Y_b <- Y[indices]
    # step 3.2: calculate initials estimator:  EigenPrism
    Eigen_Pri_b[b] <- EigenPrism(y=Y_b,X=X_b,target='beta2',diagnostics=F)$estimate
    g_single_b[b] <- mean(double_dist_sum(X_b))
    g_single_selection_b[b] <- mean(double_dist_sum(X_b[,estimated_indexes$j]))
 
     }# end of bootstrapping loop
  
  
  ########single-selection##########
  cov_emp_single <- cov(Eigen_Pri_b, g_single_b)
  c_star_mean_single <- cov_emp_single/mehane 
  ########single-selection##########
  cov_emp_sin_sel <- cov(Eigen_Pri_b, g_single_selection_b)
  c_star_mean_sel <- cov_emp_sin_sel/mehane_selection 
  
  
  
  
  return(tibble( 
    Eigenprism = Eigen_tau2
    ,T_g_tilde = Eigen_tau2 - c_star_mean_single*g_single
    ,T_h_tilde = Eigen_tau2 - c_star_mean_sel*g_selection
    
  ))
}

######################################################################################

#########################
# quick check for sim_fn:
#############################
setting_param_list <- list(seed = 1,n=n)
sim_fn(setting_param_list)


###############################################################################
ss=451
setting_param_list <-  list(seed = ss:(n_simulation+ss), n = n)

# availableCores()
result <-run_sim(sim_fn, setting_param_list, nworkers = n_cores)
arrange(result, scenario)
raw_storage<-gather(result, method, value, 4:6) 
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

 RMSE_Eigen<- stat %>% filter(method == "Eigenprism") %>%  dplyr::select(RMSE)%>%pull()
 stat <- stat %>% dplyr::select(prop_large,Estimator = method, bias, SD, RMSE, SD_RMSE 
          )%>%  mutate(tau2 = round(tau2,1), sigma2 = round(sigma2,2),
                 n_sim = n_simulation
                 , "# large betas" =c
                 , prop_large =  label_percent()(prop_large)
                 , n = n, n_pop = n_pop , total_time=total_time
                 , pct_imp = round((RMSE  /RMSE_Eigen - 1) * 100,2))


 mytime <- format(Sys.time(), "%dday_%Hhour_%Mminute_%Sseconds")



 latex_table <- stat%>%  mutate(p=n)%>%
      dplyr::select(prop_large,  tau2, n=n, Estimator, Bias=bias, SE=SD, RMSE
                    ,"% Change" = pct_imp
                    , SE.RMSE = SD_RMSE
      ) %>%  
      arrange(prop_large, factor(Estimator, levels = Estimator[c(1,2,3)]))





}

stat
