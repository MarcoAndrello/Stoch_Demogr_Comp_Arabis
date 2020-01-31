# Functions for construction and analysis of the matrix population models 
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
    predict(mod, newdat_scaled, type="response", re.form= ~(1 | Site) + (1 + Axis2 | SiteYear),allow.new.levels=T)
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
    # predict(mod, newdat_scaled, type="response", re.form= ~(1| Site) + (1 | SiteYear) + (1 + Nb_Tot.x | SiteQuad),allow.new.levels=T)
    predict(mod, newdat_scaled, type="response", re.form= ~(1| Site) + (1 | SiteYear),allow.new.levels=T) # Predicting at the site level, not at the plot level
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
build.matrices <- function(newdat, site=NULL, year=NULL, nsd_year) { 
    # Calculation of transition and fecundity matrices
    S <- s.x(newdat,site,year)
    G <- h * g.yx(newdat,site,year)
    P <- G 
    for (i in 1 : n) P[,i] = G[,i] * S[i] # growth/survival matrix
    res.f.yx    <- f.yx(newdat,site,year,nsd_year)
    pfruit      <- res.f.yx$pfruit
    log.nfruits <- res.f.yx$log.nfruits
    F           <- res.f.yx$F
    # Fix eviction in P
    # Since the growth function is bounded at 0, there is no eviction at the lower bound but only at the upper bound
    for(i in 1 : n) {
        G[n,i] <- G[n,i] + 1 - sum(G[,i])
        P[,i] <- G[,i] * S[i]
    }
    K <- P + F
    
    # Calculate lambda.det and stable-stage distribution (w)
    eigen.an <- eigen(K)
    lambda.det <- Re(eigen.an$values[1])
    w <- Re(eigen.an$vectors[,1]) / sum(Re(eigen.an$vectors[,1]))
    
    # Calculation of life-cycle components
    # Survival S
    mean.surv <- weighted.mean(colSums(P),w)
    
    # Retrogressive growth G-   
    up.G.mat <- G
    up.G.mat[lower.tri(up.G.mat,diag=T)] <- 0
    mean.retro <- weighted.mean(colSums(up.G.mat),w) 
    
    # Stasis G=
    mean.stasis <- weighted.mean(diag(G),w)
    
    # Progressive growth G+
    lo.G.mat <- G
    lo.G.mat[upper.tri(lo.G.mat,diag=T)] <- 0
    mean.progr <- weighted.mean(colSums(lo.G.mat),w) 
    
    # Reproduction F0
    mean.pfruit <- weighted.mean(pfruit,w)
    
    # Reproductive output F1
    mean.nfruits <- weighted.mean(10^log.nfruits,w) # Number of fruits conditional on fruiting
    
    return(list(G = G,
                P = P,
                F = F,
                K = K,
                pfruit = pfruit,
                log.nfruits = log.nfruits,
                lambda.det = lambda.det,
                w = w,
                mean.surv = mean.surv,
                mean.retro = mean.retro,
                mean.stasis = mean.stasis,
                mean.progr = mean.progr,
                mean.pfruit = mean.pfruit,
                mean.nfruits = mean.nfruits))
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
    
    # Calculate right and left population vectors at each time step
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
    
    # Define sensitivity and elasticity vectors
    # (for growth, it is a matrix)
    elasmat.G    <- 
        elasmat.G.mu <- 
        elasmat.G.si <-
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
    
    # Calculate mean of vital rates over years 
    G.mean <- apply(simplify2array(G), 1:2, mean)
    # F.mean <- apply(simplify2array(F), 1:2, mean)
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
        
        # # Elasticity to total fecundity
        # a <- F.mean
        # elasmat.F.mu <- elasmat.F.mu + ((v[, t + 1] %*% t(w[, t]) * a)/as.numeric((r[t] * t(v[, t + 1]) %*% w[, t + 1])))
        # a <- F[[t]] - F.mean
        # elasmat.F.si <- elasmat.F.si + ((v[, t + 1] %*% t(w[, t]) * a)/as.numeric((r[t] * t(v[, t + 1]) %*% w[, t + 1])))
        
        # Elasticity to reproduction F0 (pfruit)
        for (j in 1 : k) {
            a <- pfruit.mean[j]
            elasmat.F0.mu[j] <- elasmat.F0.mu[j] +   ( ( sum(v[,t+1] * (Fdev0[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
            a <- pfruit[j,t] - pfruit.mean[j]
            elasmat.F0.si[j] <- elasmat.F0.si[j] +   ( ( sum(v[,t+1] * (Fdev0[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
        }
        
        # Elasticity to reproductive output F1 (nfruits)
        for (j in 1 : k) {
            # a <- log.nfruits.mean[j]
            a <- 10^log.nfruits.mean[j]
            elasmat.F1.mu[j] <- elasmat.F1.mu[j] +   ( ( sum(v[,t+1] * (Fdev1[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
            # a <- log.nfruits[j,t] - log.nfruits.mean[j]
            a <- 10^log.nfruits[j,t] - 10^log.nfruits.mean[j]
            elasmat.F1.si[j] <- elasmat.F1.si[j] +   ( ( sum(v[,t+1] * (Fdev1[[t]][,j] * w[j,t])) * a) / as.numeric( (r[t] * t(v[, t + 1]) %*% w[, t + 1]) ) )
        }
        
        # Elasticity to recruit size F2 (nsd)
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
    # elasmat.F.mu <- elasmat.F.mu/tlimit
    # elasmat.F.si <- elasmat.F.si/tlimit
    
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
    # # Total Fecundity
    # E.fec.mu <- sum(elasmat.F.mu)
    # E.fec.si <- sum(elasmat.F.si)
    # Reproduction
    E.F0.mu <- sum(elasmat.F0.mu)
    E.F0.si <- sum(elasmat.F0.si)
    # Reproductive output
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
               # E.fec.mu,
               # E.fec.si,
               E.F0.mu,
               E.F0.si,
               E.F1.mu,
               E.F1.si,
               E.F2.mu,
               E.F2.si)
    wmean <- rowMeans(w[,10:(tlimit+1)]) # The first 10 time steps as burn-in
    list(elast = elast,
          w = wmean)
}



# Customized panel function:
# data contains one replicate of the bootstrap, but it is not used to calculate the correlation
# It is only used to find the [i,j] indices of the couple of variable being considered, by comparing 
# the values contained in "data" with the ones extracted by corrgram and passed to the panel function as "x" and "y" arguments
# the correlation is already stored in "med". its "significance" is assessed using "low" and "high"
panel.fct <- function(x, y, corr = NULL, col.regions, cor.method, digits = 2, cex.cor, data1, med, low, high, ...) {
    # Find i and j
    i <- which(apply(data1,2,function(z){all(x==z)}))
    j <- which(apply(data1,2,function(z){all(y==z)}))
    corr <- med[i,j] #a$estimate
    signif <- 0; if (low[i,j]>0 | high[i,j] <0) signif <- 1
    auto <- missing(cex.cor)
    usr <- par("usr")
    ncol <- 14
    pal <- col.regions(ncol)
    col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, length.out = ncol + 1), include.lowest = TRUE))
    abscorr <- formatC(abs(corr), digits = digits, format = "f")
    med <- formatC(med, digits = digits, format = "f")
    if (signif==1) col.border <- "black" else col.border <- "lightgray"
    if (auto) cex.cor <- 0.7/strwidth(abscorr)
    rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], border = NA)
    box(col = col.border)
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (signif==1) text(0.5, 0.5, "*", cex = 2.5, col = "black")
}


