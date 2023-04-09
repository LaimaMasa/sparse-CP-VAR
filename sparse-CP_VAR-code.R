jags_model <- '
model {

# 2.2 Section: Specification of shrinkage prior
for (i in 1:(K+1)) { 
  vec_deltaA[1:(N*k), i] ~ dmnorm.vcov(rep(0, N*k), Gamma[1:(N*k), 1:(N*k), i])    # def vec_deltaA per regime, i in 1 to K
}       

# when i=[1:K+1]
for (i in 1:(K+1)) {   
    for (j in 1:(N*k)) {  # for j in 1:N*k, for regime i in 1:K+1
        for (l in 1:(N*k)) {
            Gamma[j, l, i] <- ifelse(l==j, phi[i, j], 0)
        }                                                   
        phi[i, j] <- di[z[i,j]] * vari[i, j] 
        # d1=c0 and d2=c1 passed from R to JAGS
        vari[i, j] <- 1 / g[i, j]
        g[i, j] ~ dgamma(nu1/2, nu2/2) # nu1 and nu2 passed from R to JAGS
        zz[i, j] ~ dbern(1 - w[i, j])
            z[i, j] <- ifelse(zz[i, j]==0, 1, 2) 
            # JAGS does not handle index 0
        w[i, j] ~ dbeta(alpha, beta) # passed from R to JAGS
    }
}


# 2.3 Section: Sampling break dynamics
# states 
s[1] <- 1
for (t in 1:(capT)) {
    pi[t] ~ dbeta(K, (capT-K))
    ber[t] ~ dbern(pi[t])
    s[t+1] <- s[t] + ber[t]
}

# indicator function on states vector s, capT x (K+1)
for (t in 1:capT) {
 for (i in 1:(K+1)) {
     indic[t, i] <- ifelse(s[t]>=i, 1, 0)
  }
}


# 3.1 Section: Sampling the mean parameters
# manipulation step: make a vec out of vec_deltaA matrix
for (i in 1:(K+1)) {
  A_vec[((i-1)*N*k+1) : (i*N*k)] <- vec_deltaA[,i]
}

# create matrix A out of A_vec = A[1:N, 1:k_bar]
for (i in 1:k_bar) {
  A[1:N, i] <- A_vec[((i-1)*N+1) : (i*N)]
}

# create matrix Astar = N x k*(K+1), or N x k_bar
Astar[1:N, (0*k+1):(1*k)] <- A[, 1:k]
for (i in 1:K) {
  Astar[1:N, (i*k+1):((i+1)*k)] <- Astar[1:N, ((i-1)*k+1):(i*k)] + A[, (i*k+1):((i+1)*k)]
}

# create X_tilde = product between X and the indicator function, k_bar x T
tX <- t(X)
for (i in 1:(K+1)) {
  for (j in 1:k) {
    for (t in 1:capT) {
        foo[t, j, i]   <- tX[t, j] * indic[t, i]
    }
  }
  tX_tilde[1:capT, ((i-1)*k+1) : (i*k)] <- foo[1:capT, 1:k, i]
}

X_tilde <- t(tX_tilde)


# Likelihood

# Process model    A  N x k_bar    X_tilde  k_bar x capT
    # create matrix for the mean parameter
    AX_tilde <- A %*% X_tilde
    
    # create a precision matrix for Tau, invert it to obtain Sigma
    Tau[1:N, 1:N] ~ dwish(R_shape, k_df)
    Sigma[1:N, 1:N] <- inverse(Tau)

# Data model
    for (t in 1:capT) {
      Y[1:N, t] ~ dmnorm( AX_tilde[1:N, t], Tau[1:N, 1:N] )
    }

}
'

# load the required libraries libraries
library(runjags)
library(rockchalk) # mvrnorm
library(vars) # VAR
library(xlsx) # write.xlsx
library(HDInterval) # hdi
library(dplyr) # select
library(stats)

# for (l in 1:20) {
  
  folder_name <- "results/sim-study-C-120/"
  name_text <- paste0("sim-study-C-120_", 14) # l
  
  # function to simulate data from VAR(1)
  sim_VAR1 <- function(T, n, vecA, vecSigma) {
    A <- matrix(vecA, nrow = n, ncol = (n+1), byrow = FALSE) # by column
    Sigma <- matrix(vecSigma, nrow = n, ncol = n, byrow = FALSE) # by column
    y <- matrix(NA, T, n)
    y[1, ] <- A[, 1]
    set.seed(14) # seed per simulation, l
    for (t in 2:T) {
      y[t, ] <- mvrnorm(1, as.matrix(A[, 1]) + A[, 2:(n+1)] %*% y[t-1, ], Sigma)
    }
    return(y)
  }
  
  # Simulate the bivariate VAR(1) process according to DGP C parameters; see Dufays (2021)
  simT <- 120
  
  DGP_C <- rbind(
    sim_VAR1(T=floor(0.5*simT),            n=2, vecA=c(0, 0, .5, .1, .1, .5),   vecSigma = c(1, 0, 0, 1)),
    sim_VAR1(T=floor((2*simT/3-0.5*simT)), n=2, vecA=c(0, 0, -.5, .1, .1, -.5), vecSigma = c(1, 0, 0, 1)),
    sim_VAR1(T=floor((3*simT/4-2*simT/3)), n=2, vecA=c(0, 0, .5, .1, .1, .5),   vecSigma = c(1, 0, 0, 1)),
    sim_VAR1(T=floor(simT/4),              n=2, vecA=c(0, 0, -.5, .1, .1, -.5), vecSigma = c(1, 0, 0, 1)) )
  colnames(DGP_C) <- c("y1", "y2")
  
  # plot the simulated data and write the pdf file into the folder location
  pdf(file = paste0(folder_name, "sim_plot_", name_text, ".pdf"), width = 11.7, height = 8.3)
  plot(x=seq(1,length(DGP_C[,1])), y=DGP_C[,1], type = "l", col="blue",
       main = "Simulation", xlab="Time", ylab="y")
  lines(x=seq(1,length(DGP_C[,2])), y=DGP_C[,2], type = "l", col="black")
  dev.off()
  
  
  # assign simulated data to a variable 'data'
  data <- DGP_C
  
  
  # OLS VAR estimates on the full data set
  # select the VAR(p) lag according to the information criteria
  VARselect(data, lag.max = 8, type = "const")[["selection"]]
  
  # fit VAR(p) model, when p=1
  VAR_model <- VAR(data, p=1, type = "const")
  
  # obtain the VAR(1) results
  res_VAR <- summary(VAR_model)
  
  # obtain VAR(1) standard errors
  SE_y1 <- res_VAR$varresult$y1$coefficients[,2]
  SE_y2 <- res_VAR$varresult$y2$coefficients[,2]
  SE_VAR <- c(SE_y1, SE_y2) 
  
  # obtain VAR(1) parameter estimates
  E_y1 <- t(res_VAR$varresult$y1$coefficients[,1])
  E_y2 <- t(res_VAR$varresult$y2$coefficients[,1])
  E_VAR <- rbind(E_y1, E_y2)
  E_VAR <- cbind(E_VAR[, 3], E_VAR[, 1:2])
  VAR_est <- c(E_VAR)
  
  # prepare data and other parameters for JAGS
  burnin <- 20000
  sample <- 20000
  thin <- 2
  
  N <- ncol(data) # n variables in the data
  P <- 1          # p lags
  k <- (N*P) + 1  # k parameters to estimate per regime
  
  # formatting of variable Y
  Time <- nrow(data) 
  Y <- t(data[(P+1):Time,]) # drop dates and first p values
  capT <- ncol(Y)           # T, the lenght of the time series
  
  K <- round(log(capT), 0)  # K, implied number of breaks
  k_bar <- (K+1)*((N*P)+1)  # k_bar, total number of pars in CP-VAR(p)
  
  # create a matrix X out of lagged variables of Y
  X <- matrix(data = 0, nrow = (Time-P),  ncol = (N*P)) # N*P = 4*6 = 24
  for (p in 1:P) {
    X[ ,((p-1)*N+1):(p*N)] <- data[(P-(p-1)):(Time-p), ]
  }
  colnames(X) <-   rep(colnames(data), P)
  one_vec <- matrix(1, nrow = (Time-P), ncol = 1)
  X <- t(cbind(one_vec, X))
  
  
  # 2.2.1 Section: Calibration of the S&S hyper-parameters; see Dufays(2021)
  # define penalty function
  pen <- function(pi, capT) {
    val <- -( log(pi/(1-pi)) + log(capT) )
    return(val)
  }
  
  # simplify: take the mean of SEs from standard VAR model
  sig_hat <- mean(SE_VAR) 
  
  a <- (1/3) * sig_hat                
  b <- (a * (19*T*K+1)) / K
  
  pen99 <- pen(0.99, capT)
  pen75 <- pen(0.75, capT)
  
  w_99 <- a * (1-exp(pen99))
  w_75 <- a * (1-exp(pen75))
  
  # calibration of alpha and beta through the optimization function
  calib_fct <- function(th) {
    ( (pbeta(w_99, th[1], th[2]) - pbeta(w_75, th[1], th[2]) ) - 0.999)^2
  }
  
  optim_res <- optim(par=c(1,1), fn=calib_fct)
  a_cal <- optim_res$par[1]
  b_cal <- optim_res$par[2]
  
  # hyperparameters according to Malsinger-Walli & Wagner (2016)
  di <- c((1/10000), 1) # r = c0=1/10000 , c1=1
  nu <- c(10, 8)       # nu1=5 , nu2=4
  
  # R/df = expectation of cov matrix
  k_df <- N+1
  R_shape <- diag(((1^2)*k_df), N, N)
  Sig <-  R_shape/k_df # only for plotting
  
  
  # 3.1 Section: Sampling the mean parameters 
  # the code of JAGS model is defined in .txt file 
  
  vec_deltaA <- matrix(0, nrow = N*k, ncol = (K+1))
  for (i in 1:(K+1)) {
    vec_deltaA[, i] <- VAR_est
  }
  init_list <- list(vec_deltaA = vec_deltaA)
  
  # prepare the data for JAGS
  dat <- dump.format(list(Y=Y, X=X, capT=capT, N=N,
                          K=K, k=k, k_bar=k_bar,
                          alpha=a_cal, beta=b_cal,
                          di=di, nu1=nu[1], nu2=nu[2],
                          R_shape=R_shape, k_df=k_df))
  
  # initialize 3 chains with initial values and 
  # use different pseudo-random number generators
  inits1 <- dump.format(c(init_list, list( .RNG.name="base::Super-Duper", .RNG.seed=99999 )))
  inits2 <- dump.format(c(init_list, list( .RNG.name="base::Wichmann-Hill", .RNG.seed=1234 )))
  inits3 <- dump.format(c(init_list, list( .RNG.name="base::Mersenne-Twister", .RNG.seed=6666 )))
  
  # tell JAGS which latent variables to monitor
  monitor = c("A", "Astar", "Sigma", "s")
  
  # run the function that fits the models using JAGS
  results <- run.jags(jags_model, modules="glm",
                      monitor=monitor, data=dat, n.chains=3,
                      inits=c(inits1, inits2, inits3), 
                      plots = FALSE, method="parallel",
                      burnin=burnin, sample=sample, thin=thin)
  
  # read the summary of the results
  res0 <- add.summary(results)
  res <- res0$summaries
  
  # combine the MCMC chains
  chains <- rbind(results$mcmc[[1]], results$mcmc[[2]], results$mcmc[[3]])
  
  chains_A <- as.matrix(as.data.frame(chains) %>% dplyr::select(starts_with("A[")))
  chains_Astar <- as.matrix(as.data.frame(chains) %>% dplyr::select(starts_with("Astar[")))
  chains_s <- as.data.frame(chains) %>% dplyr::select(starts_with("s["))
  
  s <- t(rbind(apply(chains_s, 2, median), hdi(chains_s, 0.95))) # median
  colnames(s) <- c("median", "lower", "upper")
  
  tau <- rep(0, (nrow(s)))
  for (i in 1:(nrow(s)-1)) {
    tau[i+1] <- ifelse( s[i+1, 1] > s[i, 1], (i+1), 0)
  }
  cp <- as.data.frame(tau[tau>0])
  colnames(cp) <- "cp"
  
  # summarize inputs
  inputs <- data.frame(
    label = c("simT", "capT", "N" ,"P", "K", "k", "k_bar", "sig_hat", "a_cal", "b_cal", "nu1", "nu2", "c0", "c1", "burnin", "sample", "thin"),
    par = c(simT, capT, N, P, K, k, k_bar, sig_hat, a_cal, b_cal, nu[1], nu[2], di[1], di[2], burnin, sample, thin))
  
  
  # write all summary results in a workbook
  xlsx_name <- paste0(folder_name, "results_", name_text, ".xlsx")
  write.xlsx(as.data.frame(res), file = xlsx_name,
             sheetName = "res_table", append = FALSE)
  write.xlsx(cp, file = xlsx_name,
             sheetName="cp", append=TRUE)
  write.xlsx(inputs, file = xlsx_name,
             sheetName="inputs", append=TRUE)
  write.xlsx(data, file = xlsx_name,
             sheetName="data", append=TRUE)
  write.xlsx(Y, file = xlsx_name, 
             sheetName="Y", append=TRUE)
  write.xlsx(X, file = xlsx_name,
             sheetName="X", append=TRUE)
  
  
  # write all the cahins in excel (optional)
  # write.xlsx(chains, file = paste0(folder_name, "chains_", name_text, ".xlsx"), sheetName="chains", append=FALSE)
  
  
  # plot the results and write the pdf files into the folder location
  par_names <- colnames(chains)
  
  # chains
  pdf(file = paste0(folder_name, "chains_all_", name_text, ".pdf"), width = 8.3, height = 11.7)
  par(mfrow=c(12,2), mai=c(0.3, 0.6, 0.1, 0.3)) #c(0.02, 0.6, 0.02, 0.02))
  for (i in 1:length(par_names)) {
    plot(x=seq(1,nrow(chains)), y=chains[, i], type = "l", col="black",
         ylab=paste(par_names[i])) #xaxt='n'
    abline(v=sample, col="darkgrey")
    abline(v=(2*sample), col="darkgrey") # separate the three parallel chains
  }
  dev.off()
  
  # densities with 95% HDI 
  pdf(file = paste0(folder_name, "dens_all_", name_text, ".pdf"), width = 8.3, height = 11.7)
  par(mfrow=c(12,2), mai=c(0.3, 0.6, 0.1, 0.3))
  for (i in 1:length(par_names)) {
    plot(density(chains[, i]), ylab=paste(par_names[i]), main=" ", cex=0.75, lwd=1.5)
    abline(v=res[i, 2], col="black", lwd=1.5)
    abline(v=res[i, 1], col="red")
    abline(v=res[i, 3], col="red")
    lines(density(results[["mcmc"]][[1]][,i]), col="darkgrey")
    lines(density(results[["mcmc"]][[2]][,i]), col="darkgrey")
    lines(density(results[["mcmc"]][[3]][,i]), col="darkgrey")
  }
  dev.off()
  
  # ACF plot
  pdf(file = paste0(folder_name, "acf_A_", name_text, ".pdf"), width = 8.3, height = 11.7)
  par(mfrow=c(8,2), mai=c(0.3, 0.6, 0.1, 0.3))
  for (i in 1:length(colnames(chains_A))) {
    plot(acf(chains_A[, i], plot = FALSE), main="", ylab=paste(colnames(chains_A)[i]))
  }
  dev.off()
  
  pdf(file = paste0(folder_name, "acf_Astar_", name_text, ".pdf"), width = 8.3, height = 11.7)
  par(mfrow=c(8,2), mai=c(0.3, 0.6, 0.1, 0.3))
  for (i in 1:length(colnames(chains_Astar))) {
    plot(acf(chains_Astar[, i], plot = FALSE), main="", ylab=paste(colnames(chains_Astar)[i]))
  }
  dev.off()
  
  # summary results plot (chains, ECDF, density, ACF)
  pdf(file = paste0(folder_name, "sum_res_", name_text, ".pdf"), width = 8.3, height = 8.3)
  par(mfrow=c(10,2), mai = c(0.02, 0.6, 0.02, 0.02))
  plot(res0)
  dev.off()

# }