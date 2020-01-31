# Functions to calculate the deterministic intrinsic growth rate and deterministic elasticities to lower-level vital rates 
# They areused in the SLTRE
# There is a Main version and a Seed bank version

# Main version
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



# SEED BANK model
calc.elast.vital.rates_SeedBank <- function(surv, growth, F0, F1, F2, germ){
  
  n <- length(surv)  # Determine number of stages and matrix dimension
  dim.mat <- n + 1
  
  # Since the input F1 is the log.nfruits
  F1 <- 10^F1
  
  # Reshape growth into a matrix
  growthMat <- matrix(growth,nrow=n)
  
  # Define coefficient of conversion from fruits to recruits
  epsilon <- 0.02
  
  # Construct matrix A', following the equations in the Appendix of the manuscript
  A <- matrix(NA,nrow=dim.mat, ncol=dim.mat)
  A[1,1] <- (1-germ)
  for (k in 2 : dim.mat) {
    A[k,1] <- germ * F2[k-1]
  }
  for (l in 2 : dim.mat) {
    A[1,l] <- (1-germ) * epsilon * F0[l-1] * F1[l-1]
  }
  for (k in 2 : dim.mat) {
    for (l in 2 : dim.mat) {
      A[k,l] <- surv[l-1]*growthMat[k-1,l-1] + germ * epsilon * F0[l-1] * F1[l-1] * F2[k-1]
    }
  }
  
  # Calculating sensitivities, lambda and log-lambda of the A' matrix
  sens <- sensitivity(A)
  lambda <- Re(eigen(A)$values[1])
  r = log(lambda)
  
  # Calculating sensitivities to vital rates following the equations in the Appendix
  # Sensitivity to surv
  Sens_surv <- rep(0,n)
  for (j in 1 : n){
    for (k in 2 : dim.mat){
      Sens_surv[j] <- Sens_surv[j] + sens[k,j+1]*growthMat[k-1,j]
    }
  }
  
  # Sensitivity to growth
  Sens_growthMat <- matrix(0,nrow=n,ncol=n)
  for (i in 1 : n) { 
    for (j in 1 : n){
      Sens_growthMat[i,j] <- sens[i+1,j+1]*surv[j]
    }
  }
  Sens_growth <- as.vector(Sens_growthMat) #because growth was vectoriez in the same way (g11, g21, g31 etc.)
  
  # Sensitivity to reproduction
  Sens_F0 <- rep(0,n)
  for (j in 1 : n){
    Sens_F0[j] <- sens[1,j+1] * (1-germ) * epsilon *F1[j]
    for (k in 2 : dim.mat){
      Sens_F0[j] <- Sens_F0[j] + sens[k,j+1] * germ * epsilon * F1[j] * F2[k-1]
    }
  }
  
  # Sensitivity to reproductive output
  Sens_F1 <- rep(0,n)
  for (j in 1 : n){
    Sens_F1[j] <- sens[1,j+1] * (1-germ) * epsilon *F0[j]
    for (k in 2 : dim.mat){
      Sens_F1[j] <- Sens_F1[j] + sens[k,j+1] * germ * epsilon * F0[j] * F2[k-1]
    }
  }
  
  # Sensitivity to recruit size
  Sens_F2 <- rep(0,n)
  for (i in 1 : n){
    Sens_F2[i] <- sens[i+1,1] * germ
    for (l in 2 : dim.mat) {
      Sens_F2[i] <- Sens_F2[i] + sens[i+1,l] * germ * epsilon * F0[l-1] * F1[l-1]
    }
  }
  
  # Sensitivity to germination
  Sens_germ <- -sens[1,1]
  for (k in 2 : dim.mat){
    Sens_germ <- Sens_germ + sens[k,1] * F2[k-1]
  }
  for (l in 2 : dim.mat){
    Sens_germ <- Sens_germ - sens[1,l] * epsilon * F0[l-1] * F1[l-1] # Note the "minus" sign !
  }
  for (k in 2 : dim.mat){
    for (l in 2 : dim.mat) {
      Sens_germ <- Sens_germ + sens[k,l] * epsilon * F0[l-1] * F1[l-1] * F2[k-1]
    }
  }
  
  # Elasticities
  ES <- surv / lambda * Sens_surv
  EG <- growth / lambda * Sens_growth
  EF0 <- F0 / lambda * Sens_F0 # Note that EF0 and EF1 are equal (but the sensitivities are not)
  EF1 <- F1 / lambda * Sens_F1 # It is fine, it can be shown analytically 
  EF2 <- F2 / lambda * Sens_F2
  Egerm <- germ / lambda * Sens_germ
  
  # Formatting results for output
  elast = c(ES, EG, EF0, EF1, EF2, Egerm)
  list(r = r, elast = elast)
}
