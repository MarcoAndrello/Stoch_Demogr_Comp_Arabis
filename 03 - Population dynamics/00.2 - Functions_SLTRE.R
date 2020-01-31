# Stochastic Life-Table Response Experiment (SLTRE) (Davison et al. 2013 Am Nat)
# considering F0, F1 and F2 and not F
# (The version considering F and not F0, F1 and F2 is used only as a basis for developing the SLTRE in presence of a seed bank
# and can be found in Old code\Other - Population dynamics\SLTRE_using_F.R)

# There are three version of the function:
 # - Main version
 # - Seed bank version: to perform the SLTRE for the Seed bank model
 # - Different reference population: to test how the SLTRE changes when a different reference populaiton is choses



# ..................................................Helper functions
# Helper functions
# - to reduce G and P matrices (50 to 4 stages)
# - to reduce F0, F1 and F2 vectors (50 to 4 stages)
# - to calculate log.lambda through Tuljapurkar's approximation
# These functions are the same for the three versions of the SLTRE function

# Reduce G and P matrices
reduce.mat <- function(mat.in) {
  P <- mat.in
  P1 <- cbind(rowMeans(P[,1:5]),
              rowMeans(P[,6:10]),
              rowMeans(P[,11:15]),
              rowMeans(P[,16:dim(P)[2]])
  )
  A <- rbind(colSums(P1[1:5,]),
             colSums(P1[6:10,]),
             colSums(P1[11:15,]),
             colSums(P1[16:dim(P)[1],])
  )
  A
}

# reduce F0, F1 and F2 vectors
# For F0 and F1, we take the mean over the five stages composing the same bin;
reduce.F <- function(F.in) {
  c(mean(F.in[1:5]),
    mean(F.in[6:10]),
    mean(F.in[11:15]),
    mean(F.in[16:length(F.in)])
  )
}
# For F2 (nsd), we take the sum because it is a frequency distribution
reduce.F2 <- function(F.in) {
  c(sum(F.in[1:5]),
    sum(F.in[6:10]),
    sum(F.in[11:15]),
    sum(F.in[16:length(F.in)])
  )
}

# Calculate log stochastic growth rate
# Calculates the log stochastic growth rate by Tuljapukar approximation.
# Modified from popbio::stoch.growth.rate by commenting the simulation part
stoch.growth.rate.fast <- function (matrices, prob = NULL, maxt = 50000, verbose = TRUE) 
{
  if (is.list(matrices)) {
    matrices <- matrix(unlist(matrices), ncol = length(matrices))
  }
  s <- sqrt(dim(matrices)[1])
  n <- dim(matrices)[2]
  if (is.null(prob)) {
    prob <- rep(1/n, n)
  }
  Abar <- numeric(s^2)
  Exy <- numeric(s^4)
  for (i in 1:n) {
    A <- matrices[, i]
    Exy <- Exy + prob[i] * kronecker(A, A)
    Abar <- Abar + prob[i] * A
  }
  C <- (Exy - kronecker(Abar, Abar)) * n/(n - 1)
  C <- matrix(C, nrow = s^2)
  Abar <- matrix(Abar, nrow = s)
  ev <- eigen(Abar)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  lambda <- Re(ev$values[lmax])
  W <- ev$vectors
  w <- abs(Re(W[, lmax]))
  V <- Conj(solve(W))
  v <- abs(Re(V[lmax, ]))
  S <- v %o% w
  r <- numeric(maxt)
  n0 <- w
  # for (t in 1:maxt) {
  #     if (verbose) {
  #         if (t == 1 || t%%10000 == 0) {
  #             print(paste("Calculating stochastic growth at time", 
  #                         t), quote = FALSE)
  #         }
  #     }
  #     col <- sample(1:n, 1, prob = prob)
  #     A <- matrix(matrices[, col], nrow = s)
  #     n0 <- A %*% n0
  #     N <- sum(n0)
  #     r[t] <- log(N)
  #     n0 <- n0/N
  # }
  # loglsim <- mean(r)
  # dse <- 1.96 * sqrt(var(r)/maxt)
  # CI <- c(loglsim - dse, loglsim + dse)
  Svec <- matrix(S, ncol = 1)
  tau2 <- t(Svec) %*% C %*% Svec
  loglams <- log(lambda) - tau2/(2 * lambda^2)
  # stoch <- list(approx = as.numeric(loglams), sim = loglsim, 
  #               sim.CI = CI)
  # stoch
  return(as.numeric(loglams))
}


# ......................................................................................Main function

calc.SLTRE <- function(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd, reduce.matr=TRUE) {
  
  num.sites <- length(matrices.G)
  num.years <- length(matrices.G[[1]])
  
  if (reduce.matr == T) {
    k <- 4 # Size of the reduced matrices
    
    # Transform lists into arrays and reduce matrices to k stages
    G <- P <- array(NA,c(k, k, num.sites, num.years))
    for (i.site in 1 : num.sites){
      for(i.year in 1 : num.years) {
        G[,,i.site,i.year] <- reduce.mat(matrices.G[[i.site]][[i.year]])
        P[,,i.site,i.year] <- reduce.mat(matrices.P[[i.site]][[i.year]])
      }
    }
    rm(matrices.G, matrices.P)
    
    # Reduce F0, F1 and F2
    F0 <- F1 <- F2 <- array(NA,c(k, num.sites, num.years))
    for (i.site in 1 : num.sites){
      for(i.year in 1 : num.years) {
        F0[,i.site,i.year] <- reduce.F(pfruit.all[[i.site]][,i.year])
        F1[,i.site,i.year] <- reduce.F(log.nfruits.all[[i.site]][,i.year])
        F2[,i.site,i.year] <- reduce.F2(nsd[[i.site]][,i.year])
      }
    }
  } else {
    k <- 50
    G <- P <- array(NA,c(k, k, num.sites, num.years))
    F0 <- F1 <- F2 <- array(NA,c(k, num.sites, num.years))
    # Transform lists into arrays
    for (i.site in 1 : num.sites){
      for(i.year in 1 : num.years) {
        G[,,i.site,i.year] <- matrices.G[[i.site]][[i.year]]
        P[,,i.site,i.year] <- matrices.P[[i.site]][[i.year]]
        F0[,i.site,i.year] <- pfruit.all[[i.site]][,i.year]
        F1[,i.site,i.year] <- log.nfruits.all[[i.site]][,i.year]
        F2[,i.site,i.year] <- nsd[[i.site]][,i.year]
      }
    }
  }
  
  # Calculate parameters
  surv <- array(NA,c(k,num.sites,num.years)) # k, year, site
  g <- array(NA,c((k*k),num.sites,num.years))
  for (i.site in 1 : num.sites){
    for(i.year in 1 : num.years) {
      surv[,i.site,i.year] <- colSums(P[,,i.site,i.year])
      g[,i.site,i.year] <- as.vector(G[,,i.site,i.year]) # g11, g21, g31, g41, g12, ... g44
    }
  }
  num.par <- k + (k*k) + k + k + k # Stage-specific surv(j), gamma(i,j), F0(j), F1(j), F2(i)
  ### erase matrices (A, G, P, F) here?
  
  # Arrange parameters in the same array
  park <- abind(surv, g, F0, F1, F2, along=1)
  
  # Calculate parameters of the reference site (mean site)
  survR <- apply(surv, c(1,3), mean)
  gR <- apply(g, c(1,3), mean)
  F0R <- apply(F0, c(1,3), mean)
  F1R <- apply(F1, c(1,3), mean)
  F2R <- apply(F2, c(1,3), mean)
  parkR <- array(abind(survR, gR, F0R, F1R, F2R, along=1), c(num.par, 1, num.years))
  
  # Add it to the array
  park <- abind(park, parkR, along=2)
  dimnames(park) <- list(par=c(1:num.par),site=c(1:(num.sites+1)),year=c(1:num.years))
  #rm(surv,g,F0,F1,F2,survR,gR,F0R,F1R,F2R,parkR)
  
  # Calculate means over years
  parkmean <- apply(park,c(1,2),mean)
  dimnames(parkmean) <- list(par=c(1:num.par),site=c(1:(num.sites+1)))
  
  # Calculate intrinsic growth rates (r) and elasticities in the mean population
  r <- array(NA,(num.sites+1))
  e <- array(NA,c(num.par,(num.sites+1)))
  
  for (i.site in 1 : (num.sites+1)) {
    res <- calc.elast.vital.rates(surv   = parkmean[1:k,i.site],
                                  growth = parkmean[(k+1) : (k+(k*k)), i.site],
                                  F0 = parkmean[(k+(k*k)+1) : (k+(k*k)+k), i.site],
                                  F1 = parkmean[(k+(k*k)+k+1) : (k+(k*k)+k+k), i.site],
                                  F2 = parkmean[(k+(k*k)+k+k+1) : (k+(k*k)+k+k+k), i.site])
    e[,i.site] <- res$elast
    r[i.site] <- res$r
  }
  dimnames(e) <- dimnames(parkmean)
  
  # Calculate coefficient of variation over years per par, per site 
  cv <- apply(park, c(1,2), function(x) { sd(x) / mean(x) } )
  
  # Calculate correlations rho between par, per site 
  rho <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    rho[,,i.site] <- cor(t(park[,i.site,]))
  }
  
  # Matrices of ee, per site 
  ee <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ek <- e[,i.site] * matrix(1,num.par,num.par)
    el <- t(ek)
    ee[,,i.site] <- ek*el
  }
  
  # Matrices of cc per site 
  cc <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ck <- cv[,i.site] * matrix(1,num.par,num.par)
    cl <- t(ck)
    cc[,,i.site] <- ck*cl
  }
  
  # Matrices of ccrho, per site 
  ccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ccrho[,,i.site] <- cc[,,i.site] * rho[,,i.site]
  }
  
  # Differences
  diffpark <- array(NA,c(num.par,(num.sites+1)))
  diffee <-
    meanee <- 
    diffcc <-
    meancc <-
    diffrho <-
    meanrho <-
    diffccrho <-
    meanccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
  
  for (i.site in 1 : num.sites) {
    #diffpark[,i.site]    <- parkmean[,i.site] - parkmean[,(num.sites+1)]         # Davison et al. 2013 eq. 7
    diffpark[,i.site]    <- log(parkmean[,i.site]) - log(parkmean[,(num.sites+1)]) # Davison et al. 2019 eq. 2
    diffee[,,i.site]     <-  ee[,,i.site] - ee[,,(num.sites+1)]
    meanee[,,i.site]     <- (ee[,,i.site] + ee[,,(num.sites+1)]) / 2
    diffcc[,,i.site]     <-  cc[,,i.site] - cc[,,(num.sites+1)]
    meancc[,,i.site]     <- (cc[,,i.site] + cc[,,(num.sites+1)]) / 2
    diffrho[,,i.site]    <-  rho[,,i.site] - rho[,,(num.sites+1)]
    meanrho[,,i.site]    <- (rho[,,i.site] + rho[,,(num.sites+1)]) / 2
    diffccrho[,,i.site]  <-  ccrho[,,i.site] - ccrho[,,(num.sites+1)]
    meanccrho[,,i.site]  <- (ccrho[,,i.site] + ccrho[,,(num.sites+1)]) / 2
  }
  
  # Calculating contributions (C_x) and net stochastic contributions (C_net_e)
  C_mu <- array(NA,c(num.par,num.sites))
  C_e     <-
    C_ccrho <-
    C_cc    <-
    C_rho   <- array(NA,c(num.par,num.par,num.sites))
  C_net_e     <- 
    C_net_ccrho <-
    C_net_cc    <-
    C_net_rho   <- array(NA,c(num.par,num.sites))
  
  for (i.site in 1 : num.sites) {
    # Contributions
    # C_mu[,i.site]     <- diffpark[,i.site] / parkmean[,(num.sites+1)] * e[,i.site]            # Davison et al. 2013, eq B4: not relevant here?
    C_mu[,i.site]     <- diffpark[,i.site] * e[,i.site]                                        # Davison et al. 2019, eq 2
    C_e[,,i.site]     <- -0.5 *  meanccrho[,,i.site] * diffee[,,i.site] 
    C_ccrho[,,i.site] <- -0.5 * meanee[,,i.site] * diffccrho[,,i.site]            
    C_cc[,,i.site]    <- -0.5 * meanee[,,i.site] * meanrho[,,i.site] * diffcc[,,i.site]
    C_rho[,,i.site]   <- -0.5 * meanee[,,i.site] * meancc[,,i.site] * diffrho[,,i.site]
    
    # Net stochastic contributions
    C_net_e[,i.site]     <- colSums(C_e[,,i.site])
    C_net_ccrho[,i.site] <- colSums(C_ccrho[,,i.site])
    C_net_cc[,i.site]    <- colSums(C_cc[,,i.site])
    C_net_rho[,i.site]   <- colSums(C_rho[,,i.site])
  }
  
  #Name dimensions
  dimnames(C_mu)        <-
    dimnames(C_net_e)     <-
    dimnames(C_net_ccrho) <-
    dimnames(C_net_cc)    <-
    dimnames(C_net_rho)  <- list(par=c(1:num.par), site=v.site)
  
  # Bind the four types of contributions into the same array
  C_net <- abind(C_mu, C_net_e, C_net_cc, C_net_rho, along=3)
  dimnames(C_net)  <- list(par=c(1:num.par), site=v.site, stat=c("mu","e","cc","rho"))
  cnet <- aperm(C_net,c(1,3,2))
  
  # # # #  This part works only with 4 stages # # # # # # # # # 
  # Reorder parameters: surv, retro, stasis, progr, fec
  cnet.or <- cnet
  new.ord <- c(1 : 4,                 # survival
               9, 13, 14, 17, 18, 19, # retrogression g12, g13, g23, g14, g24, g34
               5, 10, 15, 20,         # stasis g11, g22, g33, g44
               6, 7, 8, 11, 12, 16,   # progression g21, g31, g41, g32, g42, g43
               21 : 24,               # F0
               25 : 28,               # F1
               29 : 32)               # F2   
  cnet.or <- cnet.or[new.ord,,]
  
  # Aggregate parameters: surv, retro, stasis, progr, F0, F1, F2. Net contribution of each life-cycle component to each population through each statistic
  cnet <- array(NA,c(7,4,num.sites))
  for (i.site in 1 : num.sites){
    cnet[1,,i.site] <- colSums(cnet.or[1:4,,i.site])
    cnet[2,,i.site] <- colSums(cnet.or[6:10,,i.site])
    cnet[3,,i.site] <- colSums(cnet.or[11:14,,i.site])
    cnet[4,,i.site] <- colSums(cnet.or[15:20,,i.site])
    cnet[5,,i.site] <- colSums(cnet.or[21:24,,i.site])
    cnet[6,,i.site] <- colSums(cnet.or[25:28,,i.site])
    cnet[7,,i.site] <- colSums(cnet.or[29:32,,i.site])
  }
  dimnames(cnet) <- list(par=c("S","G-","G=","G+","F0","F1","F2"),
                         stat=c("mu","e","cc","rho"),
                         site = v.site)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Calculating log-lambda for the reference population
  
  # Defining matrices to construct the matrix model
  # Using the parameters contained in survR, gR, etc to construct matrices for each year for the reference population
  matrices.R <- list()
  for (i.year in 1 : 6) {
    survMat <- matrix(survR[,i.year],k,k,byrow=T)
    G <- matrix(gR[,i.year],k,k)
    F <- 0.02 * matrix(F2R[,i.year],k,1) %*% (F0R[,i.year]*10^F1R[,i.year])
    matrices.R[[i.year]] <- G*survMat + F
  }
  
  # Calculating log-lambda for the reference population
  log.lambda.R <- stoch.growth.rate.fast(matrices.R)
  
  
  return(list(cnet = cnet,
              log.lambda.R = log.lambda.R))
}





# ......................................................................................SEED BANK model

calc.SLTRE_SeedBank <- function(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd, reduce.matr=TRUE, germ) {
  
  num.sites <- length(matrices.G)
  num.years <- length(matrices.G[[1]])
  
  
  k <- 4 # Size of the reduced matrices. One extra stage for seedbank will be added later
  
  # Transform lists into arrays and reduce matrices to k stages
  # Input: the 50-stage matrices, not the 51-stage
  G <- P <- array(NA,c(k, k, num.sites, num.years))
  for (i.site in 1 : num.sites){
    for(i.year in 1 : num.years) {
      G[,,i.site,i.year] <- reduce.mat(matrices.G[[i.site]][[i.year]])
      P[,,i.site,i.year] <- reduce.mat(matrices.P[[i.site]][[i.year]])
    }
  }
  #rm(matrices.G, matrices.P)
  
  # Reduce F0, F1 and F2
  # Input: the 50-stage matrices, not the 51-stage
  F0 <- F1 <- F2 <-  array(NA,c(k, num.sites, num.years))
  for (i.site in 1 : num.sites){
    for(i.year in 1 : num.years) {
      F0[,i.site,i.year] <- reduce.F(pfruit.all[[i.site]][,i.year])
      F1[,i.site,i.year] <- reduce.F(log.nfruits.all[[i.site]][,i.year])
      F2[,i.site,i.year] <- reduce.F2(nsd[[i.site]][,i.year])
    }
  }
  
  # Calculate parameters
  surv <- array(NA,c(k,num.sites,num.years)) # k, year, site
  g <- array(NA,c((k*k),num.sites,num.years))
  for (i.site in 1 : num.sites){
    for(i.year in 1 : num.years) {
      surv[,i.site,i.year] <- colSums(P[,,i.site,i.year])
      g[,i.site,i.year] <- as.vector(G[,,i.site,i.year]) # g11, g21, g31, g41, g12, ... g44
    }
  }
  num.par <- k + (k*k) + k + k + k + 1 # Stage-specific surv(j), g(i,j), F0(j), F1(j), F2(i), gamma
  ### erase matrices (A, G, P, F) here?
  
  # Arrange parameters in the same array
  park <- abind(surv, g, F0, F1, F2, along=1)
  
  # SEED BANK MODEL: add germination rate
  # will have to be modified if varying germ per site
  park <- abind(park,array(germ,c(1,num.sites,num.years)),along=1)
  
  # Calculate parameters of the reference site (mean site)
  survR <- apply(surv, c(1,3), mean)
  gR <- apply(g, c(1,3), mean)
  F0R <- apply(F0, c(1,3), mean)
  F1R <- apply(F1, c(1,3), mean)
  F2R <- apply(F2, c(1,3), mean)
  germR <- rep(mean(germ),num.years)
  parkR <- array(abind(survR, gR, F0R, F1R, F2R, germR, along=1), c(num.par, 1, num.years))
  
  # Add it to the array
  park <- abind(park, parkR, along=2)
  dimnames(park) <- list(par=c(1:num.par),site=c(1:(num.sites+1)),year=c(1:num.years))
  #rm(surv,g,F0,F1,F2,survR,gR,F0R,F1R,F2R,parkR)
  
  # Calculate means over years
  parkmean <- apply(park,c(1,2),mean)
  dimnames(parkmean) <- list(par=c(1:num.par),site=c(1:(num.sites+1)))
  
  # Calculate intrinsic growth rates (r) and elasticities in the mean (over years) population for each site
  # SEED BANK MODEL: Using a modified calc.elast.vital.rates function
  r <- array(NA,(num.sites+1))
  e <- array(NA,c(num.par,(num.sites+1)))
  
  for (i.site in 1 : (num.sites+1)) {
    
    # ## MI SONO FERMATO QUA: COSTRUISCO IL CALCOLO ELASTICITA' PER LA SEED BANK
    # surv   = parkmean[1:k,i.site]
    # growth = parkmean[(k+1) : (k+(k*k)), i.site]
    # F0 = parkmean[(k+(k*k)+1) : (k+(k*k)+k), i.site]
    # F1 = parkmean[(k+(k*k)+k+1) : (k+(k*k)+k+k), i.site]
    # F2 = parkmean[(k+(k*k)+k+k+1) : (k+(k*k)+k+k+k), i.site]
    # germ = parkmean[num.par, i.site]
    # #################
    
    
    res <- calc.elast.vital.rates_SeedBank(surv   = parkmean[1:k,i.site],
                                           growth = parkmean[(k+1) : (k+(k*k)), i.site],
                                           F0 = parkmean[(k+(k*k)+1) : (k+(k*k)+k), i.site],
                                           F1 = parkmean[(k+(k*k)+k+1) : (k+(k*k)+k+k), i.site],
                                           F2 = parkmean[(k+(k*k)+k+k+1) : (k+(k*k)+k+k+k), i.site],
                                           germ = parkmean[num.par, i.site])
    e[,i.site] <- res$elast
    r[i.site] <- res$r
  }
  dimnames(e) <- dimnames(parkmean)
  
  # Calculate coefficient of variation over years per par, per site 
  cv <- apply(park, c(1,2), function(x) { sd(x) / mean(x) } )
  
  # Calculate correlations rho between par, per site 
  rho <- array(NA,c(num.par,num.par,(num.sites+1)))
  options(warn=-1)
  for (i.site in 1 : (num.sites+1)) {
    rho[,,i.site] <- cor(t(park[,i.site,]))
  }
  options(warn=0)
  # SEED BANK: when germ does not change maong sites, its correlations are NA
  rho[is.na(rho)] <- 0 # Set NA to zero
  
  # Matrices of ee, per site 
  ee <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ek <- e[,i.site] * matrix(1,num.par,num.par)
    el <- t(ek)
    ee[,,i.site] <- ek*el
  }
  
  # Matrices of cc per site 
  cc <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ck <- cv[,i.site] * matrix(1,num.par,num.par)
    cl <- t(ck)
    cc[,,i.site] <- ck*cl
  }
  
  # Matrices of ccrho, per site 
  ccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ccrho[,,i.site] <- cc[,,i.site] * rho[,,i.site]
  }
  
  # Differences
  diffpark <- array(NA,c(num.par,(num.sites+1)))
  diffee <-
    meanee <- 
    diffcc <-
    meancc <-
    diffrho <-
    meanrho <-
    diffccrho <-
    meanccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
  
  for (i.site in 1 : num.sites) {
    #diffpark[,i.site]    <- parkmean[,i.site] - parkmean[,(num.sites+1)]         # Davison et al. 2013 eq. 7
    diffpark[,i.site]    <- log(parkmean[,i.site]) - log(parkmean[,(num.sites+1)]) # Davison et al. 2019 eq. 2
    diffee[,,i.site]     <-  ee[,,i.site] - ee[,,(num.sites+1)]
    meanee[,,i.site]     <- (ee[,,i.site] + ee[,,(num.sites+1)]) / 2
    diffcc[,,i.site]     <-  cc[,,i.site] - cc[,,(num.sites+1)]
    meancc[,,i.site]     <- (cc[,,i.site] + cc[,,(num.sites+1)]) / 2
    diffrho[,,i.site]    <-  rho[,,i.site] - rho[,,(num.sites+1)]
    meanrho[,,i.site]    <- (rho[,,i.site] + rho[,,(num.sites+1)]) / 2
    diffccrho[,,i.site]  <-  ccrho[,,i.site] - ccrho[,,(num.sites+1)]
    meanccrho[,,i.site]  <- (ccrho[,,i.site] + ccrho[,,(num.sites+1)]) / 2
  }
  
  # Calculating contributions (C_x) and net stochastic contributions (C_net_e)
  C_mu <- array(NA,c(num.par,num.sites))
  C_e     <-
    C_ccrho <-
    C_cc    <-
    C_rho   <- array(NA,c(num.par,num.par,num.sites))
  C_net_e     <- 
    C_net_ccrho <-
    C_net_cc    <-
    C_net_rho   <- array(NA,c(num.par,num.sites))
  
  for (i.site in 1 : num.sites) {
    # Contributions
    # C_mu[,i.site]     <- diffpark[,i.site] / parkmean[,(num.sites+1)] * e[,i.site]            # Davison et al. 2013, eq B4: not relevant here?
    C_mu[,i.site]     <- diffpark[,i.site] * e[,i.site]                                        # Davison et al. 2019, eq 2
    C_e[,,i.site]     <- -0.5 *  meanccrho[,,i.site] * diffee[,,i.site] 
    C_ccrho[,,i.site] <- -0.5 * meanee[,,i.site] * diffccrho[,,i.site]            
    C_cc[,,i.site]    <- -0.5 * meanee[,,i.site] * meanrho[,,i.site] * diffcc[,,i.site]
    C_rho[,,i.site]   <- -0.5 * meanee[,,i.site] * meancc[,,i.site] * diffrho[,,i.site]
    
    # Net stochastic contributions
    C_net_e[,i.site]     <- colSums(C_e[,,i.site])
    C_net_ccrho[,i.site] <- colSums(C_ccrho[,,i.site])
    C_net_cc[,i.site]    <- colSums(C_cc[,,i.site])
    C_net_rho[,i.site]   <- colSums(C_rho[,,i.site])
  }
  
  #Name dimensions
  dimnames(C_mu)        <-
    dimnames(C_net_e)     <-
    dimnames(C_net_ccrho) <-
    dimnames(C_net_cc)    <-
    dimnames(C_net_rho)  <- list(par=c(1:num.par), site=v.site)
  
  # Bind the four types of contributions into the same array
  C_net <- abind(C_mu, C_net_e, C_net_cc, C_net_rho, along=3)
  dimnames(C_net)  <- list(par=c(1:num.par), site=v.site, stat=c("mu","e","cc","rho"))
  cnet <- aperm(C_net,c(1,3,2))
  
  # # # #  This part works only with 4 stages # # # # # # # # # 
  # Reorder parameters: surv, retro, stasis, progr, fec
  cnet.or <- cnet
  new.ord <- c(1 : 4,                 # survival
               9, 13, 14, 17, 18, 19, # retrogression g12, g13, g23, g14, g24, g34
               5, 10, 15, 20,         # stasis g11, g22, g33, g44
               6, 7, 8, 11, 12, 16,   # progression g21, g31, g41, g32, g42, g43
               21 : 24,               # F0
               25 : 28,               # F1
               29 : 32,               # F2
               33)                    # germ         
  cnet.or <- cnet.or[new.ord,,]
  
  # Aggregate parameters: surv, retro, stasis, progr, F0, F1, F2,germ.
  # These are the net contribution of each life-cycle component to each population through each statistic
  # Seed bank model : 8 rows instead of 7 in cnet
  cnet <- array(NA,c(8,4,num.sites))
  for (i.site in 1 : num.sites){
    cnet[1,,i.site] <- colSums(cnet.or[1:4,,i.site])
    cnet[2,,i.site] <- colSums(cnet.or[6:10,,i.site])
    cnet[3,,i.site] <- colSums(cnet.or[11:14,,i.site])
    cnet[4,,i.site] <- colSums(cnet.or[15:20,,i.site])
    cnet[5,,i.site] <- colSums(cnet.or[21:24,,i.site])
    cnet[6,,i.site] <- colSums(cnet.or[25:28,,i.site])
    cnet[7,,i.site] <- colSums(cnet.or[29:32,,i.site])
    cnet[8,,i.site] <- colSums(cnet.or[33,,i.site,drop=F])
  }
  dimnames(cnet) <- list(par=c("S","G-","G=","G+","F0","F1","F2","germ"),
                         stat=c("mu","e","cc","rho"),
                         site = v.site)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Calculating log-lambda for the reference population
  # Following the equations in the Appendix of the manuscript
  # Using the parameters contained in survR, gR, etc to construct matrices for each year for the reference population
  matrices.R <- list()
  n <- nrow(survR)          # Determine number of stages and matrix dimension
  num.years <- ncol(survR)  # Determine number of years
  dim.mat <- n + 1
  F1R <- 10^F1R             # Since the input F1 is the log.nfruits
  epsilon <- 0.02           # coefficient of conversion from fruits to recruits
  for (i.year in 1 : 6) {
    gMatR <- matrix(gR[,i.year],nrow=n) # Reshape growth into a matrix
    matrices.R[[i.year]] <- matrix(NA,nrow=dim.mat, ncol=dim.mat)
    matrices.R[[i.year]][1,1] <- (1-germR[i.year])
    for (k in 2 : dim.mat) {
      matrices.R[[i.year]][k,1] <- germR[i.year] * F2R[k-1,i.year]
    }
    for (l in 2 : dim.mat) {
      matrices.R[[i.year]][1,l] <- (1-germR[i.year]) * epsilon * F0R[l-1,i.year] * F1R[l-1,i.year]
    }
    for (k in 2 : dim.mat) {
      for (l in 2 : dim.mat) {
        matrices.R[[i.year]][k,l] <- survR[l-1,i.year]*gMatR[k-1,l-1] + germR[i.year] * epsilon * F0R[l-1,i.year] * F1R[l-1,i.year] * F2R[k-1,i.year]
      }
    }
  }

  # Calculating log-lambda for the reference population
  log.lambda.R <- stoch.growth.rate.fast(matrices.R)
  
  
  return(list(cnet = cnet,
              log.lambda.R = log.lambda.R))
}


# ......................................................................................Reference population analysis

calc.SLTRE_SensRef <- function(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd, reduce.matr=TRUE, scen) {
  
  num.sites <- length(matrices.G)
  num.years <- length(matrices.G[[1]])
  
  if (reduce.matr == T) {
    k <- 4 # Size of the reduced matrices
    
    # Transform lists into arrays and reduce matrices to k stages
    G <- P <- array(NA,c(k, k, num.sites, num.years))
    for (i.site in 1 : num.sites){
      for(i.year in 1 : num.years) {
        G[,,i.site,i.year] <- reduce.mat(matrices.G[[i.site]][[i.year]])
        P[,,i.site,i.year] <- reduce.mat(matrices.P[[i.site]][[i.year]])
      }
    }
    rm(matrices.G, matrices.P)
    
    # Reduce F0, F1 and F2
    F0 <- F1 <- F2 <- array(NA,c(k, num.sites, num.years))
    for (i.site in 1 : num.sites){
      for(i.year in 1 : num.years) {
        F0[,i.site,i.year] <- reduce.F(pfruit.all[[i.site]][,i.year])
        F1[,i.site,i.year] <- reduce.F(log.nfruits.all[[i.site]][,i.year])
        F2[,i.site,i.year] <- reduce.F2(nsd[[i.site]][,i.year])
      }
    }
  } else {
    k <- 50
    G <- P <- array(NA,c(k, k, num.sites, num.years))
    F0 <- F1 <- F2 <- array(NA,c(k, num.sites, num.years))
    # Transform lists into arrays
    for (i.site in 1 : num.sites){
      for(i.year in 1 : num.years) {
        G[,,i.site,i.year] <- matrices.G[[i.site]][[i.year]]
        P[,,i.site,i.year] <- matrices.P[[i.site]][[i.year]]
        F0[,i.site,i.year] <- pfruit.all[[i.site]][,i.year]
        F1[,i.site,i.year] <- log.nfruits.all[[i.site]][,i.year]
        F2[,i.site,i.year] <- nsd[[i.site]][,i.year]
      }
    }
  }
  
  # Calculate parameters
  surv <- array(NA,c(k,num.sites,num.years)) # k, year, site
  g <- array(NA,c((k*k),num.sites,num.years))
  for (i.site in 1 : num.sites){
    for(i.year in 1 : num.years) {
      surv[,i.site,i.year] <- colSums(P[,,i.site,i.year])
      g[,i.site,i.year] <- as.vector(G[,,i.site,i.year]) # g11, g21, g31, g41, g12, ... g44
    }
  }
  num.par <- k + (k*k) + k + k + k # Stage-specific surv(j), gamma(i,j), F0(j), F1(j), F2(i)
  ### erase matrices (A, G, P, F) here?
  
  # Arrange parameters in the same array
  park <- abind(surv, g, F0, F1, F2, along=1)
  
  # Calculate parameters of the reference site
  # Here is where the code differ form the function used in the Main Text 
  survR <- surv[,scen,]
  gR  <- g[,scen,]
  F0R <- F0[,scen,]
  F1R <- F1[,scen,]
  F2R <- F2[,scen,]
  
  parkR <- array(abind(survR, gR, F0R, F1R, F2R, along=1), c(num.par, 1, num.years))
  
  # Add it to the array
  park <- abind(park, parkR, along=2)
  dimnames(park) <- list(par=c(1:num.par),site=c(1:(num.sites+1)),year=c(1:num.years))
  #rm(surv,g,F0,F1,F2,survR,gR,F0R,F1R,F2R,parkR)
  
  # Calculate means over years
  parkmean <- apply(park,c(1,2),mean)
  dimnames(parkmean) <- list(par=c(1:num.par),site=c(1:(num.sites+1)))
  
  # Calculate intrinsic growth rates (r) and elasticities in the mean population
  r <- array(NA,(num.sites+1))
  e <- array(NA,c(num.par,(num.sites+1)))
  
  for (i.site in 1 : (num.sites+1)) {
    res <- calc.elast.vital.rates(surv   = parkmean[1:k,i.site],
                                  growth = parkmean[(k+1) : (k+(k*k)), i.site],
                                  F0 = parkmean[(k+(k*k)+1) : (k+(k*k)+k), i.site],
                                  F1 = parkmean[(k+(k*k)+k+1) : (k+(k*k)+k+k), i.site],
                                  F2 = parkmean[(k+(k*k)+k+k+1) : (k+(k*k)+k+k+k), i.site])
    e[,i.site] <- res$elast
    r[i.site] <- res$r
  }
  dimnames(e) <- dimnames(parkmean)
  
  # Calculate coefficient of variation over years per par, per site 
  cv <- apply(park, c(1,2), function(x) { sd(x) / mean(x) } )
  
  # Calculate correlations rho between par, per site 
  rho <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    rho[,,i.site] <- cor(t(park[,i.site,]))
  }
  
  # Matrices of ee, per site 
  ee <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ek <- e[,i.site] * matrix(1,num.par,num.par)
    el <- t(ek)
    ee[,,i.site] <- ek*el
  }
  
  # Matrices of cc per site 
  cc <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ck <- cv[,i.site] * matrix(1,num.par,num.par)
    cl <- t(ck)
    cc[,,i.site] <- ck*cl
  }
  
  # Matrices of ccrho, per site 
  ccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
  for (i.site in 1 : (num.sites+1)) {
    ccrho[,,i.site] <- cc[,,i.site] * rho[,,i.site]
  }
  
  # Differences
  diffpark <- array(NA,c(num.par,(num.sites+1)))
  diffee <-
    meanee <- 
    diffcc <-
    meancc <-
    diffrho <-
    meanrho <-
    diffccrho <-
    meanccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
  
  for (i.site in 1 : num.sites) {
    #diffpark[,i.site]    <- parkmean[,i.site] - parkmean[,(num.sites+1)]         # Davison et al. 2013 eq. 7
    diffpark[,i.site]    <- log(parkmean[,i.site]) - log(parkmean[,(num.sites+1)]) # Davison et al. 2019 eq. 2
    diffee[,,i.site]     <-  ee[,,i.site] - ee[,,(num.sites+1)]
    meanee[,,i.site]     <- (ee[,,i.site] + ee[,,(num.sites+1)]) / 2
    diffcc[,,i.site]     <-  cc[,,i.site] - cc[,,(num.sites+1)]
    meancc[,,i.site]     <- (cc[,,i.site] + cc[,,(num.sites+1)]) / 2
    diffrho[,,i.site]    <-  rho[,,i.site] - rho[,,(num.sites+1)]
    meanrho[,,i.site]    <- (rho[,,i.site] + rho[,,(num.sites+1)]) / 2
    diffccrho[,,i.site]  <-  ccrho[,,i.site] - ccrho[,,(num.sites+1)]
    meanccrho[,,i.site]  <- (ccrho[,,i.site] + ccrho[,,(num.sites+1)]) / 2
  }
  
  # Calculating contributions (C_x) and net stochastic contributions (C_net_e)
  C_mu <- array(NA,c(num.par,num.sites))
  C_e     <-
    C_ccrho <-
    C_cc    <-
    C_rho   <- array(NA,c(num.par,num.par,num.sites))
  C_net_e     <- 
    C_net_ccrho <-
    C_net_cc    <-
    C_net_rho   <- array(NA,c(num.par,num.sites))
  
  for (i.site in 1 : num.sites) {
    # Contributions
    # C_mu[,i.site]     <- diffpark[,i.site] / parkmean[,(num.sites+1)] * e[,i.site]            # Davison et al. 2013, eq B4: not relevant here?
    C_mu[,i.site]     <- diffpark[,i.site] * e[,i.site]                                        # Davison et al. 2019, eq 2
    C_e[,,i.site]     <- -0.5 *  meanccrho[,,i.site] * diffee[,,i.site] 
    C_ccrho[,,i.site] <- -0.5 * meanee[,,i.site] * diffccrho[,,i.site]            
    C_cc[,,i.site]    <- -0.5 * meanee[,,i.site] * meanrho[,,i.site] * diffcc[,,i.site]
    C_rho[,,i.site]   <- -0.5 * meanee[,,i.site] * meancc[,,i.site] * diffrho[,,i.site]
    
    # Net stochastic contributions
    C_net_e[,i.site]     <- colSums(C_e[,,i.site])
    C_net_ccrho[,i.site] <- colSums(C_ccrho[,,i.site])
    C_net_cc[,i.site]    <- colSums(C_cc[,,i.site])
    C_net_rho[,i.site]   <- colSums(C_rho[,,i.site])
  }
  
  #Name dimensions
  dimnames(C_mu)        <-
    dimnames(C_net_e)     <-
    dimnames(C_net_ccrho) <-
    dimnames(C_net_cc)    <-
    dimnames(C_net_rho)  <- list(par=c(1:num.par), site=v.site)
  
  # Bind the four types of contributions into the same array
  C_net <- abind(C_mu, C_net_e, C_net_cc, C_net_rho, along=3)
  dimnames(C_net)  <- list(par=c(1:num.par), site=v.site, stat=c("mu","e","cc","rho"))
  cnet <- aperm(C_net,c(1,3,2))
  
  # # # #  This part works only with 4 stages # # # # # # # # # 
  # Reorder parameters: surv, retro, stasis, progr, fec
  cnet.or <- cnet
  new.ord <- c(1 : 4,                 # survival
               9, 13, 14, 17, 18, 19, # retrogression g12, g13, g23, g14, g24, g34
               5, 10, 15, 20,         # stasis g11, g22, g33, g44
               6, 7, 8, 11, 12, 16,   # progression g21, g31, g41, g32, g42, g43
               21 : 24,               # F0
               25 : 28,               # F1
               29 : 32)               # F2   
  cnet.or <- cnet.or[new.ord,,]
  
  # Aggregate parameters: surv, retro, stasis, progr, F0, F1, F2. Net contribution of each life-cycle component to each population through each statistic
  cnet <- array(NA,c(7,4,num.sites))
  for (i.site in 1 : num.sites){
    cnet[1,,i.site] <- colSums(cnet.or[1:4,,i.site])
    cnet[2,,i.site] <- colSums(cnet.or[6:10,,i.site])
    cnet[3,,i.site] <- colSums(cnet.or[11:14,,i.site])
    cnet[4,,i.site] <- colSums(cnet.or[15:20,,i.site])
    cnet[5,,i.site] <- colSums(cnet.or[21:24,,i.site])
    cnet[6,,i.site] <- colSums(cnet.or[25:28,,i.site])
    cnet[7,,i.site] <- colSums(cnet.or[29:32,,i.site])
  }
  dimnames(cnet) <- list(par=c("S","G-","G=","G+","F0","F1","F2"),
                         stat=c("mu","e","cc","rho"),
                         site = v.site)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  # Calculating log-lambda for the reference population
  
  # Defining matrices to construct the matrix model
  # Using the parameters contained in survR, gR, etc to construct matrices for each year for the reference population
  matrices.R <- list()
  for (i.year in 1 : 6) {
    survMat <- matrix(survR[,i.year],k,k,byrow=T)
    G <- matrix(gR[,i.year],k,k)
    F <- 0.02 * matrix(F2R[,i.year],k,1) %*% (F0R[,i.year]*10^F1R[,i.year])
    matrices.R[[i.year]] <- G*survMat + F
  }
  
  # Calculating log-lambda for the reference population
  log.lambda.R <- stoch.growth.rate.fast(matrices.R)
  
  
  return(list(cnet = cnet,
              log.lambda.R = log.lambda.R))
}



