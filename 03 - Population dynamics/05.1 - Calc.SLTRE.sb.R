# rm(list=ls())
# 
# library(tidyverse)
# library(gridExtra)
# library(fields)
# library(RColorBrewer)
# library(abind)
# library(corrgram)
# library(popbio)

# Function to reduce matrices: group by 5
# It is used to calculate the SLTRE
reduce.mat.sb <- function(mat.in) {
    P <- mat.in
    P1 <- cbind(rowMeans(P[,2:6]),
                rowMeans(P[,7:11]),
                rowMeans(P[,12:16]),
                rowMeans(P[,17:dim(P)[2]])
    )
    A <- rbind(colSums(P1[2:6,]),
               colSums(P1[7:11,]),
               colSums(P1[12:16,]),
               colSums(P1[17:dim(P)[1],])
    )
    A <- rbind(P1[1,],A)
    B <- rbind(P[1,1],
               sum(P[2:6,1]),
               sum(P[7:11,1]),
               sum(P[12:16,1]),
               sum(P[17:dim(P)[1],1]))
    res <- cbind(B,A)
    res
}

# Function to calculate elasticities with a seed bank stage
calc.elast.vital.rates.sb <- function(surv, growth, fec, germ){
    k <- length(surv)  # Determine number of stages and matrix dimension
    dim.mat <- k + 1
    
    F <- matrix(fec,dim.mat,k) # Redefine matrix F
    F <- cbind(rep(0,dim.mat),F)
    nsd <- F[2:dim.mat,2] / sum(F[2:dim.mat,2])  # Newborn size distribution
    nsdprime <- c(0,nsd)
    phiprime <- colSums(F) # Total fecundity
    g_sb <- c((1-germ),germ*nsd) # Redefine matrix G
    G <- matrix(growth,k,k) 
    G <- rbind(rep(0,k),G)
    G <- cbind(g_sb,G)
    survMat <- matrix(surv,k,k,byrow=T) # Redefine matrix survMat
    survMat <- rbind(rep(0,k),survMat)
    survMat <- cbind(1,survMat)
    sens <- sensitivity((G*survMat + F)) # Sensitivity and lambda
    lambda <- Re(eigen((G*survMat + F))$values[1])
    # Elasticity to G: only the submatrix containing the true 16 parameters of G
    EG <- ( sens[2:dim.mat,2:dim.mat] * survMat[2:dim.mat,2:dim.mat] ) * G[2:dim.mat,2:dim.mat] / lambda # Elasticity to G
    # Elasticity to F: only the submatrix containing the true 20 parameters of F
    EF <- sens[,2:dim.mat] * F[,2:dim.mat]/lambda
    # Elasticity to S: only the submatrix containing the true 4 parameters of S (we do not consider survival of the 1 stage i.e. the seed bank)
    ES <- colSums( sens[2:dim.mat,2:dim.mat] * G[2:dim.mat,2:dim.mat] ) * surv/lambda
    # Elasticity to germ
    Egerm <- ( -sens[1,1] +
                   sum(-sens[1,] * phiprime) +
                   sum(sens[,1] * nsdprime) +
                   sum(sens * matrix(nsdprime,dim.mat,dim.mat) * matrix(phiprime,dim.mat,dim.mat,byrow=T)) ) * germ / lambda
    r = log(lambda)
    elast = c(ES, as.vector(EG), as.vector(EF), Egerm)
    list(r = r, elast = elast)
}


# load("Matrices_for_seedbank_SLTRE.RData")

# Stochastic Life-Table Response Experiment (SLTRE) (Davison et al. 2013 Am Nat)
# Including a seed bank stage
calc.SLTRE.sb <- function(matrices.K, matrices.G, matrices.P, matrices.F){
    num.sites <- length(matrices.K)
    num.years <- length(matrices.K[[1]])
    
    # Transform lists into arrays and reduce matrices
    # A <- P <- G <- F <- array(NA,c(4, 4, num.sites, num.years))
    A <- P <- G <- F <- array(NA,c(5, 5, num.sites, num.years))
    for (i.site in 1 : num.sites){
        for(i.year in 1 : num.years) {
            A[,,i.site,i.year] <- reduce.mat.sb(matrices.K[[i.site]][[i.year]])
            P[,,i.site,i.year] <- reduce.mat.sb(matrices.P[[i.site]][[i.year]])
            G[,,i.site,i.year] <- reduce.mat.sb(matrices.G[[i.site]][[i.year]])
            F[,,i.site,i.year] <- reduce.mat.sb(matrices.F[[i.site]][[i.year]])
        }
    }
    rm(matrices.K, matrices.G, matrices.P, matrices.F)
    dim.mat <- dim(A)[1]
    k <-dim.mat-1 # minus one because we do not count the seed bank stage
    
    # Calculate parameters
    num.par <- 1 + k + (k*k) + (dim.mat*k) # germination + stage-specific surv(j), gamma(i,j) and fec(i,j)
    surv <- array(NA,c(k,num.sites,num.years)) # k, year, site
    g <- array(NA,c((k*k),num.sites,num.years))
    f <- array(NA,c((dim.mat*k),num.sites,num.years))
    germ <- 1 - G[1,1,,]
    for (i.site in 1 : num.sites){
        for(i.year in 1 : num.years) {
            surv[,i.site,i.year] <- colSums(P[2:dim.mat, 2:dim.mat, i.site, i.year])
            g[,i.site,i.year] <- as.vector(G[2:dim.mat, 2:dim.mat, i.site, i.year]) # g11, g21, g31, g41, g12, ... g44
            f[,i.site,i.year] <- as.vector(F[, 2:dim.mat, i.site, i.year]) # same
        }
    }
    ### erase matrices here?
    
    # Arrange parameters in the same array
    park <- abind(surv, g, f, germ, along=1)
    
    # Calculate parameters of the reference site (mean site) (for each year?)
    survR <- apply(surv, c(1,3), mean)
    gR <- apply(g, c(1,3), mean)
    fR <- apply(f, c(1,3), mean)
    germR <- apply(germ, 2, mean)
    parkR <- array(abind(survR, gR, fR, germR, along=1), c(num.par, 1, num.years))
    
    # Add it to the array
    park <- abind(park, parkR, along=2)
    dimnames(park) <- list(par=c(1:num.par),site=c(1:(num.sites+1)),year=c(1:num.years))
    rm(surv,g,f,survR,gR,fR,parkR)
    
    # Calculate means over years
    parkmean <- apply(park,c(1,2),mean)
    dimnames(parkmean) <- list(par=c(1:num.par),site=c(1:(num.sites+1)))
    
    # Calculate intrinsic growth rates (r) and elasticities in the mean population
    r <- array(NA,(num.sites+1))
    e <- array(NA,c(num.par,(num.sites+1)))
    
    # Calculate elasticities
    for (i.site in 1 : (num.sites+1)) {
        res <- calc.elast.vital.rates.sb(surv   = parkmean[1:k,i.site],
                                      growth = parkmean[(k+1) : (k+(k*k)), i.site],
                                      fec    = parkmean[((k+(k*k))+1) : (num.par-1), i.site],
                                      germ   = parkmean[num.par, i.site])
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
    new.ord <- c(1, 2, 3, 4, # survival
                 9, 13, 14, 17, 18, 19, # retrogression g12, g13, g23, g14, g24, g34
                 5, 10, 15, 20,         # stasis g11, g22, g33, g44
                 6, 7, 8, 11, 12, 16,   # progression g21, g31, g41, g32, g42, g43
                 21:40,        # fecundity  
                 41)            # germination
    cnet.or <- cnet.or[new.ord,,]
    
    # Aggregate parameters: surv, retro, stasis, progr, fec. Net contribution of each life-cycle component to each population through each statistic
    cnet <- array(NA,c(6,4,num.sites)) # 6 life-cycle components
    for (i.site in 1 : num.sites){
        cnet[1,,i.site] <- colSums(cnet.or[1:4,,i.site])
        cnet[2,,i.site] <- colSums(cnet.or[6:10,,i.site])
        cnet[3,,i.site] <- colSums(cnet.or[11:14,,i.site])
        cnet[4,,i.site] <- colSums(cnet.or[15:20,,i.site])
        cnet[5,,i.site] <- colSums(cnet.or[21:40,,i.site])
        cnet[6,,i.site] <- cnet.or[41,,i.site]
    }
    dimnames(cnet) <- list(par=c("S","G-","G=","G+","F","gamma"),
                           stat=c("mu","e","cc","rho"),
                           site = v.site)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    return(cnet)
}
