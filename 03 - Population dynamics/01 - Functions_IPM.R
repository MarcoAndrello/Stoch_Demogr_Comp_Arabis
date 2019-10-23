# Functions for construction and analysis of the  matrix population models 
# Marco Andrello
# 23/10/2019

# Function to append climatic data of the following yeat to year t
extend.env.dataset <- function(env0) {
    env.onlyclim <- ungroup(env0) %>% select(-c(Site, Year))
    env0 <-  mutate(env0, UniqueID = paste0(Site,"-",Year))  
    env_dat_follow <- env_dat_mean <- array(NA,dim(env.onlyclim))
    for (i in 1 : dim(env0)[1]) {
        current_UniqueID <- env0$UniqueID[i]
        current_year <- env0$Year[i]
        follow_year <- current_year + 1
        follow_UniqueID <- current_UniqueID
        substr(follow_UniqueID,5,8) <- as.character(follow_year)
        id <- which(env0$UniqueID == follow_UniqueID)
        env_dat_follow[i,] <- as.numeric(env.onlyclim[id,])
        env_dat_mean[i,] <- ( as.numeric(env.onlyclim[i,]) + as.numeric(env.onlyclim[id,]) ) / 2
    }
    colnames(env_dat_follow) <- paste0(colnames(env.onlyclim),"_f")
    colnames(env_dat_mean) <- paste0(colnames(env.onlyclim),"_m")
    env0 <- cbind(as.data.frame(env0),env_dat_follow,env_dat_mean)
    env0 <- as_tibble(env0) %>% select(-UniqueID)
    rm(current_UniqueID,follow_UniqueID,current_year,follow_year,env_dat_follow,env_dat_mean)
    return(env0)
}

# Function to sample models
sample_model <- function(object,weights) {
    n.models <- length(object)
    i.model <- sample(n.models,1,replace=T,prob=weights)
    object[i.model][[1]]
}

# Function to scale predictors
scale_predictors <- function(newdat,predictors.mean,predictors.sd) {
    predictors.mean %>% names %>% paste0(".or") %>% match(names(newdat)) -> id.var
    newdat %>% select(id.var) %>% scale(predictors.mean,predictors.sd) %>% as_tibble -> out
    names(out) <- names(predictors.mean)
    out
}

# Survival function
s.x <- function(newdat,site=NULL,year=NULL) {
    newdat %>%
        scale_predictors(predictors.surv.mean,predictors.surv.sd) %>%
        mutate(Site=site,Year=year,SiteYear=paste0(Site,Year)) ->
        newdat_scaled
    mod <- sample_model(d.surv,m.weight.surv)
    predict(mod, newdat_scaled, type="response", re.form= ~(1 | Site/SiteYear),allow.new.levels=T)
}

# Growth function
g.yx <- function(newdat,site=NULL,year=NULL) {
    newdat %>%
        scale_predictors(predictors.growth.mean,predictors.growth.sd) %>%
        mutate(Site=site,Year=year,SiteYear=paste0(Site,Year)) ->
        newdat_scaled
    mod <- sample_model(d.growth,m.weight.growth)
    ypred <- predict(mod, newdat_scaled, re.form= ~(1 | Site/SiteYear), allow.new.levels=T, type="response")
    theta <- getME(mod,"glmer.nb.theta")
    g <- matrix(NA,nrow=length(ypred),ncol=length(ypred))
    for (k in 1 : length(ypred)) {
        g[,k] <- dnbinom((y-1),mu=ypred[k],size=theta)
    }
    return(g)
}

# Recruit size function
newb.size.distr <- function(newdat,site=NULL,year=NULL) {
    newdat %>%
        slice(1) %>%
        mutate(T.mean.or = T.mean_f.or, T.mean.sq.or = T.mean_f.sq.or, T.range.or = T.range_f.or) %>%  # because we are predicting the nsd for the following year
        scale_predictors(predictors.recruit.mean,predictors.recruit.sd) %>%
        mutate(Site=site,Year=year,SiteYear=paste0(Site,Year)) ->
        newdat_scaled
    mod <- sample_model(d.recruit,m.weight.recruit)
    ypred <- predict(mod, newdat_scaled, type="response", re.form= ~(1 | Site/SiteYear),allow.new.levels=T)
    theta <- getME(mod,"glmer.nb.theta")
    dnbinom((y-1), mu=ypred, size=theta)
}

# Reproduction (probability of bearing fruits) function
p.fruit <- function(newdat,site=NULL,year=NULL) {
    newdat %>%
        scale_predictors(predictors.fruiting.mean,predictors.fruiting.sd) %>%
        mutate(Site=site,Year=year,SiteYear=paste0(Site,Year)) ->
        newdat_scaled
    mod <- sample_model(d.fruiting,m.weight.fruiting)
    predict(mod, newdat_scaled, type="response", re.form= ~(1 | Site/SiteYear),allow.new.levels=T)
}

# Fecundity (number of fruits)  function
n.fruits <- function(newdat,site=NULL,year=NULL) {
    newdat %>%
        scale_predictors(predictors.fruits.mean,predictors.fruits.sd) %>%
        mutate(Site=site,Year=year,SiteYear=paste0(Site,Year)) ->
        newdat_scaled
    mod <- sample_model(d.fruits,m.weight.fruits)
    ypred <- predict(mod, newdat_scaled, re.form= ~(1 | Site/SiteYear),allow.new.levels=T)
    pfruit <- p.fruit(newdat,site,year)
    out <- list(pfruit  = pfruit,
                log.nfruits = ypred,
                nfruits = pfruit * (10^ypred)
                )
    out
}

# Overall Fecundity function
f.yx <- function(newdat,site=NULL,year=NULL,nsd_year) {
    res.nfruits <- n.fruits(newdat,site,year)
    out <- list(pfruit = res.nfruits$pfruit,
                log.nfruits = res.nfruits$log.nfruits,
                F = 0.02 * matrix(nsd_year,dim(newdat)[1],1) %*% res.nfruits$nfruits # 0.02 is the median establishment coefficient
                )
    out
}

# Function to build the matrix model
build.IPM <- function(newdat,site=NULL,year=NULL,nsd_year) { 
    # Calculation of transition and fercundity matrices
    S <- s.x(newdat,site,year)
    G <- h * g.yx(newdat,site,year)
    P <- G 
    for(i in 1 : n) P[,i] = G[,i] * S[i] # growth/survival matrix
    res.f.yx    <- f.yx(newdat,site,year,nsd_year)
    pfruit      <- res.f.yx$pfruit
    log.nfruits <- res.f.yx$log.nfruits
    F           <- res.f.yx$F
    # Fix eviction in P
    # Since the growth function is bounded at 0, there is no eviction at the lower bound
    for(i in 1 : n) {
        G[n,i] <- G[n,i] + 1 - sum(G[,i])
        P[,i] <- G[,i] * S[i]
    }
    K <- P + F
    return(list(G = G,
                P = P,
                F = F,
                K = K,
                pfruit = pfruit,
                log.nfruits = log.nfruits))
}



# Elasticities of a = log lambda_s to changes in the mean and standard deviation of vital rates
# This function is adapted from the "stoch.sens" function of popbio (itself copied from the demogR package)
# according to equations given in Haridas & Tuljapurkar (2005 Am. Nat.) and Caswell (Aust. N. Z. J. Stat. 47(1), 2005, 75-85)
stoch.sens.Marco <- function (A, G, P, F, 
                              pfruit, log.nfruits, nsd,
                              tlimit = 100) {
    
    # These lines are for debugging:
    # tlimit=100 
    # matrices.K -> A
    # matrices.G -> G
    # matrices.P -> P
    # matrices.F -> F
    
    epsilon <- 0.02 # 0.02 is the median establishment coefficient
    
    k <- ncol(A[[1]])
    num.years <- length(A)
    
    # Calculate survival per year from P matrix
    surv <- array(NA,c(k,num.years))
    for (i.year in 1 : num.years){
        surv[,i.year] <- colSums(P[[i.year]])
    }
    
    # Calculate matrices of parameter values (surv, p.fruit, log.nfruits and recr.size)
    # Calculate derivatives of p(i,j) to surv and growth (equal to derivatives of A(i,j) to these vital rates)
    # Calculate derivatives of f(i,j) to p.fruit, nfruits and recr.size  (equal to derivatives of A(i,j) to these vital rates)
    S <- F0 <- F1 <- F2 <- 
        Sdev <- Gdev <- Fdev0 <- Fdev1 <- Fdev2 <- list()
    for (i.year in 1 : num.years){
        S[[i.year]] <- matrix(rep(surv[,i.year],k), k, k, byrow=T)
        F0[[i.year]] <- matrix(rep(pfruit[,i.year],k), k, k, byrow=T)
        F1[[i.year]] <- matrix(rep(log.nfruits[,i.year],k), k, k,byrow=T)
        F2[[i.year]] <- matrix(rep(nsd[,i.year],k), k, k, byrow=FALSE)
        Sdev[[i.year]] <- G[[i.year]] # Derivative of A to survival
        Gdev[[i.year]] <- S[[i.year]] # Derivative of A to growth
        Fdev0[[i.year]] <- epsilon * nsd[,i.year] %*% t(10^log.nfruits[,i.year]) # Derivative of A to F0
        # Fdev1[[i.year]] <- epsilon * log(10) * nsd[,i.year] %*% t( pfruit[,i.year] * 10^log.nfruits[,i.year] ) # derivative of A to the log F1; not used
        Fdev1[[i.year]] <- epsilon * nsd[,i.year] %*% t(pfruit[,i.year]) # Derivative of A to F1
        Fdev2[[i.year]] <- epsilon * rep(1,k) %*% t( pfruit[,i.year] * 10^log.nfruits[,i.year] ) # Derivative of A to F2
    }
    
    # Sample matrices and vectors
    sequence <- sample(num.years, tlimit, replace = TRUE)
    A <- A[sequence]
    S <- S[sequence]
    G <- G[sequence]
    F <- F[sequence]
    F0 <- F0[sequence]
    F1 <- F1[sequence]
    F2 <- F2[sequence]
    Sdev <- Sdev[sequence]
    Gdev <- Gdev[sequence]
    #   Fdev is equal to 1
    Fdev0 <- Fdev0[sequence]
    Fdev1 <- Fdev1[sequence]
    Fdev2 <- Fdev2[sequence]
    surv <- surv[,sequence]
    pfruit <- pfruit[,sequence]
    log.nfruits <- log.nfruits[,sequence]
    nsd <- nsd[,sequence]
    
    # Calculate right and left eigenvectors at each time step
    wvec <- rep(1/k, k)
    w <- cbind(wvec)
    r <- rep(0, tlimit)
    for (i in 1:tlimit) {
        a <- A[[i]]
        wvec <- a %*% wvec
        r[i] <- sum(wvec)
        wvec <- wvec/r[i]
        w <- cbind(w, wvec)
    }
    vvec <- rep(1/k, k)
    v <- cbind(vvec)
    for (i in rev(1:tlimit)) {
        a <- A[[i]]
        vvec <- vvec %*% a
        v <- cbind(t(vvec), v)
    }
    
    # Define sensitivity and elasticity matrices
    # (for survival, it is a vector)
    elasmat.G    <- 
        elasmat.G.mu <- 
        elasmat.G.si <-
        elasmat.F    <- 
        elasmat.F.mu <- 
        elasmat.F.si <- 
        matrix(0, nrow = k, ncol = k)
    elasmat.S    <- 
        elasmat.S.mu <- 
        elasmat.S.si <- 
        elasmat.F0    <- 
        elasmat.F0.mu <- 
        elasmat.F0.si <- 
        elasmat.F1    <- 
        elasmat.F1.mu <- 
        elasmat.F1.si <- 
        elasmat.F2    <- 
        elasmat.F2.mu <- 
        elasmat.F2.si <- array(0, k)
    
    # Calculate mean over years of vital rates
    G.mean <- apply(simplify2array(G), 1:2, mean)
    F.mean <- apply(simplify2array(F), 1:2, mean)
    surv.mean <- rowMeans(surv)
    pfruit.mean <- rowMeans(pfruit)
    log.nfruits.mean <- rowMeans(log.nfruits)
    nsd.mean <- rowMeans(nsd)
    
    for (t in 1:tlimit) {
        
        # Elasticity to survival
        for (j in 1 : k) {
            a <- surv.mean[j]
            elasmat.S.mu[j] <- elasmat.S.mu[j] +   ( ( sum(v[,t+1] * (Sdev[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
            a <- surv[j,t] - surv.mean[j]
            elasmat.S.si[j] <- elasmat.S.si[j] +   ( ( sum(v[,t+1] * (Sdev[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
        }
        
        # Elasticity to growth
        a <- G.mean * Gdev[[t]]
        elasmat.G.mu <- elasmat.G.mu + ((v[, t + 1] %*% t(w[, t]) * a)/as.numeric((r[t] * t(v[, t + 1]) %*% w[, t + 1])))
        a <- ( G[[t]] - G.mean )  * Gdev[[t]]
        elasmat.G.si <- elasmat.G.si + ((v[, t + 1] %*% t(w[, t]) * a)/as.numeric((r[t] * t(v[, t + 1]) %*% w[, t + 1])))
        
        # Elasticity to total fecundity
        a <- F.mean
        elasmat.F.mu <- elasmat.F.mu + ((v[, t + 1] %*% t(w[, t]) * a)/as.numeric((r[t] * t(v[, t + 1]) %*% w[, t + 1])))
        a <- F[[t]] - F.mean
        elasmat.F.si <- elasmat.F.si + ((v[, t + 1] %*% t(w[, t]) * a)/as.numeric((r[t] * t(v[, t + 1]) %*% w[, t + 1])))
        
        # Elasticity to reproduction (F0)
        for (j in 1 : k) {
            a <- pfruit.mean[j]
            elasmat.F0.mu[j] <- elasmat.F0.mu[j] +   ( ( sum(v[,t+1] * (Fdev0[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
            a <- pfruit[j,t] - pfruit.mean[j]
            elasmat.F0.si[j] <- elasmat.F0.si[j] +   ( ( sum(v[,t+1] * (Fdev0[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
        }
        
        # Elasticity to nfruits (F1)
        for (j in 1 : k) {
            # a <- log.nfruits.mean[j]
            a <- 10^log.nfruits.mean[j]
            elasmat.F1.mu[j] <- elasmat.F1.mu[j] +   ( ( sum(v[,t+1] * (Fdev1[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
            # a <- log.nfruits[j,t] - log.nfruits.mean[j]
            a <- 10^log.nfruits[j,t] - 10^log.nfruits.mean[j]
            elasmat.F1.si[j] <- elasmat.F1.si[j] +   ( ( sum(v[,t+1] * (Fdev1[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
        }
        
        # Elasticity to nsd (F2)(recruit size)
        for (i in 1 : k) {
            a <- nsd.mean[i]
            elasmat.F2.mu[i] <- elasmat.F2.mu[i] +   ( ( sum(v[i,t+1] * (Fdev2[[t]][i,] * w[,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
            a <- nsd[i,t] - nsd.mean[i]
            elasmat.F2.si[i] <- elasmat.F2.si[i] +   ( ( sum(v[i,t+1] * (Fdev2[[t]][i,] * w[,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
        }
        
    }
    
    elasmat.S.mu <- elasmat.S.mu/tlimit
    elasmat.S.si <- elasmat.S.si/tlimit
    elasmat.G.mu <- elasmat.G.mu/tlimit
    elasmat.G.si <- elasmat.G.si/tlimit
    elasmat.F0.mu <- elasmat.F0.mu/tlimit
    elasmat.F0.si <- elasmat.F0.si/tlimit
    elasmat.F1.mu <- elasmat.F1.mu/tlimit
    elasmat.F1.si <- elasmat.F1.si/tlimit
    elasmat.F2.mu <- elasmat.F2.mu/tlimit
    elasmat.F2.si <- elasmat.F2.si/tlimit
    elasmat.F.mu <- elasmat.F.mu/tlimit
    elasmat.F.si <- elasmat.F.si/tlimit
    
    # Sum elasticities by life-cycle component
    # Survival
    E.surv.mu <- sum(elasmat.S.mu)
    E.surv.si <- sum(elasmat.S.si)
    # Retrogressive growth
    temp <- elasmat.G.mu
    temp[lower.tri(temp,diag=T)] <- 0
    E.retro.growth.mu <- sum(temp)
    temp <- elasmat.G.si
    temp[lower.tri(temp,diag=T)] <- 0
    E.retro.growth.si <- sum(temp)
    # Stasis
    temp <- elasmat.G.mu
    E.stasis.mu <- sum(diag(temp))
    temp <- elasmat.G.si
    E.stasis.si <- sum(diag(temp))
    # Progressive growth
    temp <- elasmat.G.mu
    temp[upper.tri(temp,diag=T)] <- 0
    E.progr.growth.mu <- sum(temp)
    temp <- elasmat.G.si
    temp[upper.tri(temp,diag=T)] <- 0
    E.progr.growth.si <- sum(temp)
    # Total Fecundity
    E.fec.mu <- sum(elasmat.F.mu)
    E.fec.si <- sum(elasmat.F.si)
    # Reproduction
    E.F0.mu <- sum(elasmat.F0.mu)
    E.F0.si <- sum(elasmat.F0.si)
    # Fruit set
    E.F1.mu <- sum(elasmat.F1.mu)
    E.F1.si <- sum(elasmat.F1.si)
    # Recruit size
    E.F2.mu <- sum(elasmat.F2.mu)
    E.F2.si <- sum(elasmat.F2.si)
    
    elast <- c(E.surv.mu,
               E.surv.si,
               E.retro.growth.mu,
               E.retro.growth.si,
               E.stasis.mu,
               E.stasis.si, 
               E.progr.growth.mu,
               E.progr.growth.si,
               E.fec.mu,
               E.fec.si,
               E.F0.mu,
               E.F0.si,
               E.F1.mu,
               E.F1.si,
               E.F2.mu,
               E.F2.si)
    
    wmean <- rowMeans(w[,10:(tlimit+1)])
    list(elast=elast,
         w=wmean)
}



# Functions to calculate survival, retrogressive growth, stasis, progressive growth, total fecundity
# from matrices
calc.life.cycle.comp <- function(mat.K, mat.G, mat.F, mat.P, pfruit, log.nfruits) {
  
    # Calculate stable-stage distribution (w)
    eigen.an <- eigen(mat.K)
    lambda.det <- Re(eigen.an$values[1])
    w <- Re(eigen.an$vectors[,1]) / sum(Re(eigen.an$vectors[,1]))
    
    # Survival S
    mean.surv <- weighted.mean(colSums(mat.P),w)
    
    # Progressive growth G+
    lo.G.mat <- mat.G
    lo.G.mat[upper.tri(lo.G.mat,diag=T)] <- 0
    mean.progr <- weighted.mean(colSums(lo.G.mat),w) 

    # Retrogressive growth G-   
    up.G.mat <- mat.G
    up.G.mat[lower.tri(up.G.mat,diag=T)] <- 0
    mean.retro <- weighted.mean(colSums(up.G.mat),w) 
    
    # Stasis G=
    mean.stasis <- weighted.mean(diag(mat.G),w)
    
    # Total fecundity
    mean.fec <- weighted.mean(colSums(mat.F),w)
    
    # Reproduction F0
    mean.pfruit <- weighted.mean(pfruit,w)
    
    # Fecundity F1
    mean.nfruits     <- weighted.mean(10^log.nfruits,w) # Number of fruits conditional on fruiting
    mean.nfruits.net <- weighted.mean((pfruit*(10^log.nfruits)),w) # Number of fruits: not used
    
    return(list(lambda.det = lambda.det,
                w = w,
                mean.surv = mean.surv,
                mean.progr = mean.progr,
                mean.retro = mean.retro,
                mean.stasis = mean.stasis,
                mean.fec = mean.fec,
                mean.pfruit = mean.pfruit,
                mean.nfruits = mean.nfruits,
                mean.nfruits.net = mean.nfruits.net))
}

# Stochastic life expectancy following Tuljapurkar and Horvitz (2006)
# Equations for random i.i.d. environments
stoch.life.expectancy <- function(matrices.P) {
  # Mean
  n <- dim(matrices.P[[1]])[1]
  qbar <- matrix(0,nrow=n,ncol=n)
  for (i in 1 : length(matrices.P)){
    qbar <- qbar + matrices.P[[i]]
  }
  qbar <- qbar / length(matrices.P)   # Eq. 17 (equal probability pi_a of choosing a matrix)
  Nbar <- ginv(diag(n) - qbar)        # Eq. 19
  lexp.mean <- colSums(Nbar)
  
  # Variance
  eta_a <- eta_squared_a <- var_eta_a <- array(NA,dim=c(length(matrices.P),n))
  
  for (i in 1 : length(matrices.P)){
    Na <- diag(n) + Nbar %*% matrices.P[[i]]                        # Eq. 18: Nbar == (I-qbar)^-1
    eta_a[i,] <- colSums(Na)                                        # Eq. 20
    eta_squared_a[i,] <- colSums(Na + 2 * Nbar %*% (Na - diag(n)))  # Eq. 22
    var_eta_a[i,] <- eta_squared_a[i,] - eta_a[i,]^2                # Definition of variance
  }
  lexp.var <- colMeans(var_eta_a)
  
  stoch.life.expectancy <- data.frame(mean=lexp.mean, sd=sqrt(lexp.var))
  return(stoch.life.expectancy)
}




# Function to calculate the deterministic intrinsic growth rate and deterministic elasticities to lower-level vital rates 
# It is used to calculate the SLTRE
calc.elast.vital.rates <- function(surv, growth, F0, F1, F2){
    
    # Number of stages
    k <- length(surv)
    
    # Defining matrices to construct the matrix model
    survMat <- matrix(surv,k,k,byrow=T)
    G <- matrix(growth,k,k)
    F <- 0.02 * matrix(F2,k,1) %*% (F0*10^F1)
    A <- G*survMat + F
    
    # Calculating sensitivities, lambda and log-lambda
    sens <- sensitivity(A)
    lambda <- Re(eigen(A)$values[1])
    r = log(lambda)
    
    # Calculating elasticities to lower-level vital rates using the chain rule
    ES <- colSums(sens*G) * surv/lambda
    EG <- (sens*survMat) * G/lambda
    EF0 <- EF1 <- EF2 <- vector()
    for (j in 1 : k) {
        EF0[j] <- EF1[j] <- sum( F[,j] * sens[,j] ) / lambda # It can be shown analytically that EF0 and EF1 are equal (but the sensitivities are not)
    }
    for (i in 1 : k) {
        EF2[i] <- sum( F[i,] * sens[i,] ) / lambda
    }

    # Formatting results for output
    elast = c(ES, as.vector(EG), EF0, EF1, EF2)
    list(r = r, elast = elast)
}



# Stochastic Life-Table Response Experiment (SLTRE) (Davison et al. 2013 Am Nat)
# Two versions: 
# 1) considering F0, F1 and F2 and not F
# 2) considering F and not F0, F1 and F2; used only as a basis for developing the SLTRE in presence of a seed bank


# 1) considering F0, F1 and F2 and not F

# Helper functions to reduce G and P matrices and to aggregate F0, F1 and F2 over 50 stages into 4 stages
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

# For F0 and F1, we take the mean over the five stages composing the same bin;
reduce.F <- function(F.in) {
    c(mean(F.in[1:5]),
      mean(F.in[6:10]),
      mean(F.in[11:15]),
      mean(F.in[16:length(F.in)])
    )
}

# For nsd, we take the sum because it is a frequency distribution
reduce.F2 <- function(F.in) {
    c(sum(F.in[1:5]),
      sum(F.in[6:10]),
      sum(F.in[11:15]),
      sum(F.in[16:length(F.in)])
    )
}


calc.SLTRE <- function(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd) {
    
    # Size of the reduced matrices
    k <- 4
    
    num.sites <- length(matrices.G)
    num.years <- length(matrices.G[[1]])
    
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
    rm(surv,g,F0,F1,F2,survR,gR,F0R,F1R,F2R,parkR)
    
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
        diffpark[,i.site]    <- parkmean[,i.site] - parkmean[,(num.sites+1)]
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
        C_mu[,i.site] <- diffpark[,i.site] / parkmean[,(num.sites+1)] * e[,i.site]            # Davison et al. 2013, eq B4
        C_e[,,i.site]     <- -0.5 *  meanccrho[,,i.site] * diffee[,,i.site] 
        C_ccrho[,,i.site] <- -0.5 * meanee[,,i.site] * diffccrho[,,i.site]            
        C_cc[,,i.site]    <- -0.5 * meanee[,,i.site] * meanrho[,,i.site] * diffcc[,,i.site]
        C_rho[,,i.site]   <- -0.5 * meanee[,,i.site] * meancc[,,i.site] * diffrho[,,i.site]
        
        # Net stochastic contributions
        C_net_e[,i.site]     <- colSums(C_e[,,i.site])
        C_net_ccrho[,i.site] <- colSums(C_ccrho[,,i.site])
        C_net_cc[,i.site]    <-colSums(C_cc[,,i.site])
        C_net_rho[,i.site]   <-colSums(C_rho[,,i.site])
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
    
    return(cnet)
}



# # 2) Considering F and not F0, F1 and F2
# calc.SLTRE <- function(matrices.K, matrices.G, matrices.P, matrices.F){
#     num.sites <- length(matrices.K)
#     num.years <- length(matrices.K[[1]])
#     
#     # Transform lists into arrays and reduce matrices
#     A <- P <- G <- F <- array(NA,c(4, 4, num.sites, num.years))
#     for (i.site in 1 : num.sites){
#         for(i.year in 1 : num.years) {
#             A[,,i.site,i.year] <- reduce.mat(matrices.K[[i.site]][[i.year]])
#             P[,,i.site,i.year] <- reduce.mat(matrices.P[[i.site]][[i.year]])
#             G[,,i.site,i.year] <- reduce.mat(matrices.G[[i.site]][[i.year]])
#             F[,,i.site,i.year] <- reduce.mat(matrices.F[[i.site]][[i.year]])
#         }
#     }
#     rm(matrices.K, matrices.G, matrices.P, matrices.F)
#     k <- dim(A)[1]
#     
#     # Calculate parameters
#     num.par <- k + (k*k) + (k*k) # Stage-specific surv(j), gamma(i,j) and fec(i,j)
#     surv <- array(NA,c(k,num.sites,num.years)) # k, year, site
#     g <-
#         f <- array(NA,c((k*k),num.sites,num.years))
#     for (i.site in 1 : num.sites){
#         for(i.year in 1 : num.years) {
#             surv[,i.site,i.year] <- colSums(P[,,i.site,i.year])
#             g[,i.site,i.year] <- as.vector(G[,,i.site,i.year]) # g11, g21, g31, g41, g12, ... g44
#             f[,i.site,i.year] <- as.vector(F[,,i.site,i.year]) # same
#         }
#     }
#     ### erase matrices here?
#     
#     # Arrange parameters in the same array
#     park <- abind(surv, g, f, along=1)
#     
#     # Calculate parameters of the reference site (mean site)
#     survR <- apply(surv, c(1,3), mean)
#     gR <- apply(g, c(1,3), mean)
#     fR <- apply(f, c(1,3), mean)
#     parkR <- array(abind(survR, gR, fR, along=1), c(num.par, 1, num.years))
#     
#     # Add it to the array
#     park <- abind(park, parkR, along=2)
#     dimnames(park) <- list(par=c(1:num.par),site=c(1:(num.sites+1)),year=c(1:num.years))
#     rm(surv,g,f,survR,gR,fR,parkR)
#     
#     # Calculate means over years
#     parkmean <- apply(park,c(1,2),mean)
#     dimnames(parkmean) <- list(par=c(1:num.par),site=c(1:(num.sites+1)))
#     
#     # Calculate intrinsic growth rates (r) and elasticities in the mean population
#     r <- array(NA,(num.sites+1))
#     e <- array(NA,c(num.par,(num.sites+1)))
#     for (i.site in 1 : (num.sites+1)) {
#         res <- calc.elast.vital.rates(surv   = parkmean[1:k,i.site],
#                                       growth = parkmean[(k+1) : (k+(k*k)), i.site],
#                                       fec    = parkmean[((k+(k*k))+1) : num.par, i.site])
#         e[,i.site] <- res$elast
#         r[i.site] <- res$r
#     }
#     dimnames(e) <- dimnames(parkmean)
#     
#     # Calculate coefficient of variation over years per par, per site 
#     cv <- apply(park, c(1,2), function(x) { sd(x) / mean(x) } )
#     
#     # Calculate correlations rho between par, per site 
#     rho <- array(NA,c(num.par,num.par,(num.sites+1)))
#     for (i.site in 1 : (num.sites+1)) {
#         rho[,,i.site] <- cor(t(park[,i.site,]))
#     }
#     
#     # Matrices of ee, per site 
#     ee <- array(NA,c(num.par,num.par,(num.sites+1)))
#     for (i.site in 1 : (num.sites+1)) {
#         ek <- e[,i.site] * matrix(1,num.par,num.par)
#         el <- t(ek)
#         ee[,,i.site] <- ek*el
#     }
#     
#     # Matrices of cc per site 
#     cc <- array(NA,c(num.par,num.par,(num.sites+1)))
#     for (i.site in 1 : (num.sites+1)) {
#         ck <- cv[,i.site] * matrix(1,num.par,num.par)
#         cl <- t(ck)
#         cc[,,i.site] <- ck*cl
#     }
#     
#     # Matrices of ccrho, per site 
#     ccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
#     for (i.site in 1 : (num.sites+1)) {
#         ccrho[,,i.site] <- cc[,,i.site] * rho[,,i.site]
#     }
#     
#     
#     
#     # Differences
#     diffpark <- array(NA,c(num.par,(num.sites+1)))
#     diffee <-
#         meanee <- 
#         diffcc <-
#         meancc <-
#         diffrho <-
#         meanrho <-
#         diffccrho <-
#         meanccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
#     for (i.site in 1 : num.sites) {
#         diffpark[,i.site]    <- parkmean[,i.site] - parkmean[,(num.sites+1)]
#         diffee[,,i.site]     <-  ee[,,i.site] - ee[,,(num.sites+1)]
#         meanee[,,i.site]     <- (ee[,,i.site] + ee[,,(num.sites+1)]) / 2
#         diffcc[,,i.site]     <-  cc[,,i.site] - cc[,,(num.sites+1)]
#         meancc[,,i.site]     <- (cc[,,i.site] + cc[,,(num.sites+1)]) / 2
#         diffrho[,,i.site]    <-  rho[,,i.site] - rho[,,(num.sites+1)]
#         meanrho[,,i.site]    <- (rho[,,i.site] + rho[,,(num.sites+1)]) / 2
#         diffccrho[,,i.site]  <-  ccrho[,,i.site] - ccrho[,,(num.sites+1)]
#         meanccrho[,,i.site]  <- (ccrho[,,i.site] + ccrho[,,(num.sites+1)]) / 2
#     }
#     
#     C_mu <- array(NA,c(num.par,num.sites))
#     C_e     <-
#         C_ccrho <-
#         C_cc    <-
#         C_rho   <- array(NA,c(num.par,num.par,num.sites))
#     C_net_e     <- 
#         C_net_ccrho <-
#         C_net_cc    <-
#         C_net_rho   <- array(NA,c(num.par,num.sites))
#     
#     for (i.site in 1 : num.sites) {
#         # Contributions
#         C_mu[,i.site] <- diffpark[,i.site] / parkmean[,(num.sites+1)] * e[,i.site]            # Davison et al. 2013, eq B4
#         C_e[,,i.site]     <- -0.5 *  meanccrho[,,i.site] * diffee[,,i.site] 
#         C_ccrho[,,i.site] <- -0.5 * meanee[,,i.site] * diffccrho[,,i.site]            
#         C_cc[,,i.site]    <- -0.5 * meanee[,,i.site] * meanrho[,,i.site] * diffcc[,,i.site]
#         C_rho[,,i.site]   <- -0.5 * meanee[,,i.site] * meancc[,,i.site] * diffrho[,,i.site]
#         
#         # Net stochastic contributions
#         C_net_e[,i.site]     <- colSums(C_e[,,i.site])
#         C_net_ccrho[,i.site] <- colSums(C_ccrho[,,i.site])
#         C_net_cc[,i.site]    <-colSums(C_cc[,,i.site])
#         C_net_rho[,i.site]   <-colSums(C_rho[,,i.site])
#     }
#     
#     #Name dimensions
#     dimnames(C_mu)        <-
#         dimnames(C_net_e)     <-
#         dimnames(C_net_ccrho) <-
#         dimnames(C_net_cc)    <-
#         dimnames(C_net_rho)  <- list(par=c(1:num.par), site=v.site)
#     
#     # Bind the four types of contributions into the same array
#     C_net <- abind(C_mu, C_net_e, C_net_cc, C_net_rho, along=3)
#     dimnames(C_net)  <- list(par=c(1:num.par), site=v.site, stat=c("mu","e","cc","rho"))
#     cnet <- aperm(C_net,c(1,3,2))
# 
#     # # # #  This part works only with 4 stages # # # # # # # # # 
#     # Reorder parameters: surv, retro, stasis, progr, fec
#     cnet.or <- cnet
#     new.ord <- c(1, 2, 3, 4, # survival
#                  9, 13, 14, 17, 18, 19, # retrogression g12, g13, g23, g14, g24, g34
#                  5, 10, 15, 20,         # stasis g11, g22, g33, g44
#                  6, 7, 8, 11, 12, 16,   # progression g21, g31, g41, g32, g42, g43
#                  21:36)        # fecundity  
#     cnet.or <- cnet.or[new.ord,,]
#     
#     # Aggregate parameters: surv, retro, stasis, progr, fec. Net contribution of each life-cycle component to each population through each statistic
#     cnet <- array(NA,c(5,4,num.sites))
#     for (i.site in 1 : num.sites){
#         cnet[1,,i.site] <- colSums(cnet.or[1:4,,i.site])
#         cnet[2,,i.site] <- colSums(cnet.or[6:10,,i.site])
#         cnet[3,,i.site] <- colSums(cnet.or[11:14,,i.site])
#         cnet[4,,i.site] <- colSums(cnet.or[15:20,,i.site])
#         cnet[5,,i.site] <- colSums(cnet.or[21:36,,i.site])
#     }
#     dimnames(cnet) <- list(par=c("S","G-","G=","G+","F"),
#                            stat=c("mu","e","cc","rho"),
#                            site = v.site)
#     # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#     
#     return(cnet)
# }




