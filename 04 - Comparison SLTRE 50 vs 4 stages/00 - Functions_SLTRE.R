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

# for nsd, we take the sum because it is a frequency distribution
reduce.F2 <- function(F.in) {
  c(sum(F.in[1:5]),
    sum(F.in[6:10]),
    sum(F.in[11:15]),
    sum(F.in[16:length(F.in)])
  )
}


# calculate elasticities to lower level paramters
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



