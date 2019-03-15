# Data: Feb.28.2019
# ST 562 Probabilistic Graphical Model Project
# Member: Laurel Hopkins, Eugene Seo, Chuan Tian
# Code description: Simulation data for N-mixture model

library("unmarked")
library("sigmoid")
nSites <- 10 #500 # number of total sites 
nVisits <- 1 # number of total visits
nFeatures <- 10 # number of features (Assuming it's same for both habitat and detection probability)
nRepeat <- 1 #50
nK <- 200
case <-1
cat("nSites:", nSites, ", nVisits:", nVisits, ", nFeatures:", nFeatures, "\n")

a1 <- runif(nFeatures + 1) # generate random parameters for habitat features (z1)
a2 <- runif(nFeatures + 1) # generate random parameters for habitat features (z2)
b1 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w1)
b2 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w2)


#if (file.exists("result_case_1.csv")) {
#  file.remove("result_case_1.csv")
#  file.remove("result_case_2.csv")
#  file.remove("result_case_3.csv")  
#}

for(i in 1:nRepeat) {
  start_time <- Sys.time()
  cat("Repeat", i, "\n")  
  # Generating simulation data for a species A ==================================
  z1 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures) # generate random habitat features
  
  # Generating simulation data for a species B ==================================
  z2 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures) # generate random habitat features
  
  #for(case in 1:3){
    result <- array(0, c(1,4))
    cat(" Generating data ( case ", case, ") >>> ")
    if ( case == 1 ) {
      w1_dummy <- z1 # w1 = z1
      w2_dummy <- z2 # w2 = z2
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
        w1[[i]] <- as.matrix(w1_dummy[,i])
        w2[[1]] <- as.matrix(w2_dummy[,i])
      }
      w1_dummyXb1 = sweep(cbind(1,w1_dummy),MARGIN=2,b1,`*`)
      p <- sigmoid(w1_dummyXb1) # detection probability
      w2_dummyXb2 = sweep(cbind(1,w2_dummy),MARGIN=2,b2,`*`)
      fp <- sigmoid(w2_dummyXb2) # false positive rate  
      
    } else if (case == 2) {
      w1_dummy <- matrix(runif(nSites*nFeatures), nSites, nFeatures)       
      w2_dummy <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
      sub.idx <- sample(1:nFeatures, round(nFeatures/2))
      w1_dummy[,sub.idx] <- z1[,sub.idx] # w1 partially = z1
      w2_dummy[,sub.idx] <- z2[,sub.idx] # w2 partially = z2  ### DOES THIS MAKE SENSE? 
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      
      for (i in 1:nFeatures){
        w1[[i]] <- as.matrix(w1_dummy[,i])
        w2[[i]] <- as.matrix(w2_dummy[,i])
      }
      
      w1_dummyXb1 = sweep(cbind(1,w1_dummy),MARGIN=2,b1,`*`)
      p <- sigmoid(w1_dummyXb1) # detection probability
      w2_dummyXb2 = sweep(cbind(1,w2_dummy),MARGIN=2,b2,`*`)
      fp <- sigmoid(w2_dummyXb2) # false positive rate  

    } else {  
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (j in 1:nFeatures){
        w1[[i]] <- matrix(runif(nSites*nVisits), nSites, nVisits) # w1 != z1
        w2[[i]] <- matrix(runif(nSites*nVisits), nSites, nVisits) # w2 != z2
      }
    }

    lambda1 <- exp(cbind(1,z1) %*% a1) # parameter for the true abundance of A 
    N1 <- rpois(nSites, lambda1) # true abundance for species A

    lambda2 <- exp(cbind(1,z2) %*% a2) # parameter for the true abundance of B
    N2 <- rpois(nSites, lambda2) # true abundance for species B

    # Generating observed counts Y ================================================
    Y <- matrix(NA, nSites, nVisits) # observed counts for A
    A <- matrix(NA, nSites, nVisits) # true observed counts from A
    B <- matrix(NA, nSites, nVisits) # false observed counts from B (mis-identification)
    for(j in 1:nVisits) {
      A[,j] <- rbinom(nSites, N1, p[j])  # true observed counts from A
      B[,j] <- rbinom(nSites, N2, fp[j]) # false observed counts from B (mis-identification)
      Y[,j] <- A[,j] + B[,j] # total observed counts recorded as A
    }
    cat("Max Y:", max(Y), ">>> ")
    
    if ( max(Y) < 200 ) {
      cat("Fitting models >>> \n")
      # Fit a model with A
      umf1 <- unmarkedFramePCount(y=A, siteCovs=data.frame(x=z1),obsCovs=w1)
      fm1 <- pcount(~V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8  + V9 + V10 
                    ~x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10, 
                    umf1, K=nK)
      b1_hat <- coef(fm1, type="det")
      W1 <- cbind(rep(1,nSites), w1)
      
      p_hat <- plogis(W1 %*% b1_hat)
      N1_hat <- round(bup(ranef(fm1))) # Estimated population size
      result[1] <- sqrt(mean((p_hat - p)^2)) # USE RMSE
      result[2] <- sqrt(mean(abs(N1_hat - N1) / N1))  ## CHUAN TO FIX # USE RMSE 

      # Fit a model with B
      visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)
      umf2 <- unmarkedFramePCount(y=B, siteCovs=data.frame(x=z2),obsCovs=w2)
      fm1 <- pcount(~V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8  + V9 + V10 
                    ~x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10, 
                    umf2, K=nK)
      b2_hat = coef(fm2, type="det")
      W2 = cbind(rep(1,nSites), w2)
      fp_hat = plogis(W2 %*% b2_hat)
      N2_hat = round(bup(ranef(fm2))) # Estimated population size
      result[3] <- sqrt(mean((fp_hat - fp)^2))  # USE
      result[4] <- sqrt(mean(abs(N1_hat - N2) / N2))  ## CHUAN TO FIX 
      write.table(result, paste("result_case_", case, ".csv", sep=""), col.names=F, row.names=F,sep=",", append=TRUE)    
    }
  #}
  end_time <- Sys.time()
  print(end_time - start_time)  
}
#print(colMeans(result))