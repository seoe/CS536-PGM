# Data: Feb.28.2019
# ST 562 Probabilistic Graphical Model Project
# Member: Laurel Hopkins, Eugene Seo, Chuan Tian
# Code description: Simulation data for N-mixture model

library("unmarked")
library("sigmoid")

# fixed parameters
nFeatures <- 10 # number of features (Assuming it's same for both habitat and detection probability)
nRepeat <- 50

# experimental parameters
nSites <- 200 # number of total sites (200~500)
nVisits <- 1 # number of total visits (1~5)
nK <- 100 # max number of observations
#case <- 1

cat("nSites:", nSites, "nVisits:", nVisits, "nK:", nK, "\n")
#set.seed("1234")

# clear the previous result files
if (file.exists("result_case_1.csv")) {
  for (i in 1:3) {
    file.remove(paste("s", nSites, "_v", nVisits, "_k", nK, "_c", i, ".csv", sep=""))
  }
}

for (i in 1:nRepeat) {
  start_time <- Sys.time()
  cat("Repeat", i, "\n")
  
  # 1. generate parameters for habitat features
  a1 <- runif(nFeatures + 1, 0, 0.8) 
  a2 <- runif(nFeatures + 1, 0, 0.8) 
  
  # 2. generate parameters for detection features
  b1 <- rnorm(nFeatures + 1) 
  b2 <- rnorm(nFeatures + 1)
   
  # 3. generate habitat features
  z1 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures)
  z2 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures)
  
  # 4. generate detection features according to the degree of overlap with habitat features
  for (case in 1:3){
    result <- array(0, c(1,6))
    cat(" Generating data (", case, ") >>> ")
    if ( case == 1 ) { # w = z
      w1_dummy <- z1 
      w2_dummy <- z2
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
        w1[[i]] <- as.matrix(w1_dummy[,i])
        w2[[i]] <- as.matrix(w2_dummy[,i])
      }
      
     } else if (case == 2) { # w partially = z
      w1_dummy <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
      w2_dummy <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
      sub.idx <- sample(1:nFeatures, round(nFeatures/2))
      w1_dummy[,sub.idx] <- z1[,sub.idx] 
      sub.idx <- sample(1:nFeatures, round(nFeatures/2))
      w2_dummy[,sub.idx] <- z2[,sub.idx]
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
       w1[[i]] <- as.matrix(w1_dummy[,i])
       w2[[i]] <- as.matrix(w1_dummy[,i])
      }
      
    } else { # w != z
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
        w1[[i]] <- matrix(runif(nSites*nVisits), nSites, nVisits) 
        w2[[i]] <- matrix(runif(nSites*nVisits), nSites, nVisits)
      }
    }

    # 5. generate Poisson parameter and true abundance of species 
    lambda1 <- exp(cbind(1,z1) %*% a1) 
    N1 <- matrix(rpois(nSites * nVisits, lambda1), nSites, nVisits)
    
    lambda2 <- exp(cbind(1,z2) %*% a2)
    N2 <- matrix(rpois(nSites * nVisits, lambda2), nSites, nVisits)
    
    # 6. generate detection probability
    ones = list()
    ones[[1]] = matrix(1, nSites, nVisits)
    W1 = c(ones, w1)
    prod = matrix(0, nSites, nVisits)
    for (i in 1:nFeatures + 1) {
      prod = prod + W1[[i]] * b1[i]
    }
    p <- sigmoid(prod) # detection probability
    
    W2 = c(ones, w2)
    prod = matrix(0, nSites, nVisits)
    for (i in 1:nFeatures + 1) {
      prod = prod + W2[[i]] * b2[i]
    }
    fp <- sigmoid(prod)  # false positive rate


    # 7. generating observed counts Y ================================================
    Y <- matrix(NA, nSites, nVisits) # observed counts for A
    A <- matrix(NA, nSites, nVisits) # true observed counts from A
    B <- matrix(NA, nSites, nVisits) # false observed counts from B (mis-identification)
    for (i in 1:nSites) {
      for (j in 1:nVisits) {
        A[i,j] <- rbinom(1, N1[i, j], p[i, j])  # true observed counts from A
        B[i,j] <- rbinom(1, N2[i, j], fp[i, j]) # false observed counts from B (mis-identification)
        Y[i,j] <- A[i,j] + B[i,j] # total observed counts recorded as A
      }
    }
    cat("Max Y:", max(Y), ">>> ")
    result[5] <- max(Y)
    result[6] <- mean(Y)    
    
    # 8. learning models
    if ( max(Y) < nK ) {
      cat("Fitting models >>> \n")
      # Fit a model with A
      z1.1 <- data.frame(x1=z1)
      NAME <- paste0("w1.", 1:nFeatures)
      names(w1) <- NAME
      umf1 <- unmarkedFramePCount(y=A, siteCovs=z1.1,obsCovs=w1)
      fm1 <- pcount(~w1.1 + w1.2 + w1.3 + w1.4 + w1.5 + w1.6 + w1.7 + w1.8  + w1.9 + w1.10
                    ~x1.1 + x1.2 + x1.3 + x1.4 + x1.5 + x1.6 + x1.7 + x1.8 + x1.9 + x1.10, 
                    umf1, K=nK)
      b1_hat <- coef(fm1, type = "det")
      a1_hat <- coef(fm1, type = "state")
      
      prod <- matrix(0, nSites, nVisits)
      for (i in 1:nFeatures + 1) {
        prod = prod + W1[[i]] * b1_hat[i]
      }
      p_hat <- sigmoid(prod) # estimated detection probability
      lambda1_hat <- exp(cbind(1, z1) %*% a1_hat)
      result[1] <- sqrt(mean((p_hat - p)^2))
      result[2] <- sqrt(mean(((lambda1_hat - lambda1) / lambda1) ^ 2)) 

      # Fit a model with B
      z2.1 <- data.frame(x2=z2)
      NAME <- paste0("w2.", 1:nFeatures)
      names(w2) <- NAME
      umf2 <- unmarkedFramePCount(y=B, siteCovs=z2.1,obsCovs=w2)
      fm2 <- pcount(~w2.1 + w2.2 + w2.3 + w2.4 + w2.5 + w2.6 + w2.7 + w2.8  + w2.9 + w2.10
                    ~x2.1 + x2.2 + x2.3 + x2.4 + x2.5 + x2.6 + x2.7 + x2.8 + x2.9 + x2.10, 
                    umf2, K=nK)
      b2_hat <- coef(fm2, type = "det")
      a2_hat <- coef(fm2, type = "state")
      
      prod = matrix(0, nSites, nVisits)
      for (i in 1:nFeatures + 1) {
        prod <- prod + W2[[i]] * b2_hat[i]
      }
      fp_hat <- sigmoid(prod) # estimated detection probability
      lambda2_hat <- exp(cbind(1, z2) %*% a2_hat)
      result[3] <- sqrt(mean((fp_hat - fp)^2))
      result[4] <- sqrt(mean(((lambda2_hat - lambda2) / lambda2) ^ 2))
      
      write.table(result, paste("s", nSites, "_v", nVisits, "_k", nK, "_c", i, ".csv", sep=""), col.names=F, row.names=F,sep=",", append=TRUE)    
    }
  }
  end_time <- Sys.time()
  print(end_time - start_time)  
}
print(colMeans(result))
