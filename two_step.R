# Date: Mar.18.2019
# ST 536 Probabilistic Graphical Model Project
# Member: Laurel Hopkins, Eugene Seo, Chuan Tian
# Code description: two-step model

library("unmarked")
library("sigmoid")
nSites <- 200 #500 # number of total sites 
nVisits <- 1 # number of total visits    In this model, nVisits always = 1
nFeatures <- 10 # number of features (Assuming it's same for both habitat and detection probability)
nRepeat <- 10 #50
nK <- 500
case <-3   # In this model, case always = 3
cat("nSites:", nSites, ", nVisits:", nVisits, ", nFeatures:", nFeatures, "\n")
result <- matrix(NA, nRepeat, 7)
#set.seed("1234")

# a1 <- runif(nFeatures+1, 0, 0.8) # generate random parameters for habitat features (z1)
# a2 <- runif(nFeatures + 1, 0, 0.8) # generate random parameters for habitat features (z2)
# b1 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w1)
# b2 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w2)


for(k in 1:nRepeat) {
  a1 <- runif(nFeatures+1, 0, 1) # generate random parameters for habitat features (z1)
  a2 <- runif(nFeatures + 1, 0, 1) # generate random parameters for habitat features (z2)
  b1 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w1)
  b2 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w2)
  start_time <- Sys.time()
  cat("Repeat", k, "\n")  
  # Generating simulation data for a species A ==================================
  z1 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures) # generate random habitat features
  
  # Generating simulation data for a species B ==================================
  z2 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures) # generate random habitat features
  
  #for(case in 1:3){
    # result <- array(0, c(1,5))
    cat(" Generating data ( case ", case, ") >>> ")
    if ( case == 1 ) {
      w1_dummy <- z1 # w1 = z1
      w2_dummy <- z2 # w2 = z2
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
        w1[[i]] <- as.matrix(w1_dummy[,i])
        w2[[i]] <- as.matrix(w2_dummy[,i])
      }
      
     } else if (case == 2) {
      w1_dummy <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
      w2_dummy <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
      sub.idx <- sample(1:nFeatures, round(nFeatures/2))
      w1_dummy[,sub.idx] <- z1[,sub.idx] # w1 partially = z1
      sub.idx <- sample(1:nFeatures, round(nFeatures/2))
      w2_dummy[,sub.idx] <- z2[,sub.idx] # w2 partially = z2
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
       w1[[i]] <- as.matrix(w1_dummy[,i])
       w2[[i]] <- as.matrix(w1_dummy[,i])
      }
      
    } else {  
      w1 <- vector("list", nFeatures)
      w2 <- vector("list", nFeatures)
      for (i in 1:nFeatures){
        w1[[i]] <- matrix(runif(nSites*nVisits), nSites, nVisits) # w1 != z1
        w2[[i]] <- matrix(runif(nSites*nVisits), nSites, nVisits) # w2 != z1
      }
    }

    lambda1 <- exp(cbind(1,z1) %*% a1) # parameter for the true abundance of A 
    N1 <- matrix(rpois(nSites * nVisits, lambda1), nSites, nVisits) # true abundance for species A -fixed for nVisits > 1 --Chuan
    ones = list()
    ones[[1]] = matrix(1, nSites, nVisits)
    W1 = c(ones, w1)
    prod = matrix(0, nSites, nVisits)
    for (i in 1:nFeatures + 1) {
      prod = prod + W1[[i]] * b1[i]
    }
    p <- sigmoid(prod) # detection probability  ## FIXED--Chuan


    lambda2 <- exp(cbind(1,z2) %*% a2) # parameter for the true abundance of B
    N2 <- matrix(rpois(nSites * nVisits, lambda2), nSites, nVisits) # true abundance for species B -fixed for nVisits > 1 --Chuan
    W2 = c(ones, w2)
    prod = matrix(0, nSites, nVisits)
    for (i in 1:nFeatures + 1) {
      prod = prod + W2[[i]] * b2[i]
    }
    fp <- sigmoid(prod)  # false positive rate  ## FIXED--Chuan


    # Generating observed counts Y ================================================
    Y <- matrix(NA, nSites, nVisits) # observed counts for A
    A <- matrix(NA, nSites, nVisits) # true observed counts from A
    B <- matrix(NA, nSites, nVisits) # false observed counts from B (mis-identification)
    for (i in 1:nSites) { # Fixed to accomondate nVisits > 1 -- Chuan
      for(j in 1:nVisits) {
        A[i,j] <- rbinom(1, N1[i, j], p[i, j])  # true observed counts from A
        B[i,j] <- rbinom(1, N2[i, j], fp[i, j]) # false observed counts from B (mis-identification)
        Y[i,j] <- A[i,j] + B[i,j] # total observed counts recorded as A
      }
    }
    cat("Max Y:", max(Y), ">>> ")
    
    if ( max(Y) < nK ) {
      # 1st step: Estimate Y1, Y2 using E-M based on two indepdent Poission models (not considering the identifiability between lamda and p yet)
      
      # Initialization
      w1.poi <- matrix(NA, nSites, nFeatures)
      w2.poi <- matrix(NA, nSites, nFeatures)
      for (i in 1:nFeatures) {
        w1.poi[, i] <- w1[[i]]
        w2.poi[, i] <- w2[[i]]
      }
      
      z <-  cbind(z1, z2)
      w.poi <- cbind(w1.poi, w2.poi)
      my_data = data.frame(y = as.vector(Y), z, w.poi)
      model_total <- glm(y ~. , family = poisson(), data = my_data)
      mu_est <- matrix(model_total$fitted, nSites, nVisits)
      
      mu1_est <- matrix(NA, nSites, nVisits)
      mu2_est <- matrix(NA, nSites, nVisits)
      A_est <- matrix(NA, nSites, nVisits)
      B_est <- matrix(NA, nSites, nVisits)
      
      for (i in 1:nSites) {
        for (j in 1:nVisits) {
          mu1_est[i, j] = runif(n = 1, min = 0, max = mu_est[i, j])
          mu2_est[i, j] = mu_est[i, j] - mu1_est[i, j]
          A_est[i, j] = round(Y[i, j] * mu1_est[i, j] / mu_est[i, j], 0)
          B_est[i, j] = Y[i, j] - A_est[i, j]
        }
      }
      
      my_data_1 <- data.frame(y1_est = as.vector(A_est), z1, w1.poi)
      my_data_2 <- data.frame(y2_est = as.vector(B_est), z2, w2.poi)
      model_1 <- glm(y1_est ~ ., family = poisson(), data = my_data_1)
      model_2 <- glm(y2_est ~ ., family = poisson(), data = my_data_2)
      mu1_est <- as.vector(model_1$fitted)
      mu2_est <- as.vector(model_2$fitted)
      mu_est = mu1_est + mu2_est
      
      coef_1 <- rep(0, 2 * nFeatures + 1)
      coef_2 <- rep(0, 2 * nFeatures + 1)
      
      
      # Iterate until convergence
      while (mean(abs(coef_1 - as.vector(model_1$coef))) > 0.001 | mean(abs(coef_2 - as.vector(model_2$coef))) > 0.001)
      {
        # E-step
        coef_1 <- as.vector(model_1$coef)
        coef_2 <- as.vector(model_2$coef)
        for (i in 1:n) {
          y1_est[i] = round(y[i] * mu1_est[i] / mu_est[i], 0)
          y2_est[i] = y[i] - y1_est[i]
        }
        
        my_data_1$y1_est = y1_est
        my_data_2$y2_est = y2_est
        
        # M-step
        model_1 <- glm(y1_est ~ ., family = poisson(), data = my_data_1)
        model_2 <- glm(y2_est ~ ., family = poisson(), data = my_data_2)
        mu1_est <- as.vector(model_1$fitted)
        mu2_est <- as.vector(model_2$fitted)
        mu_est = mu1_est + mu2_est
      }
      coef_1 <- as.vector(model_1$coef)
      coef_2 <- as.vector(model_2$coef)

    # 2nd step: Feed estimated Y1, Y2 to two independent N-mixure models
    z1.1 <- data.frame(x1=z1)
    NAME <- paste0("w1.", 1:nFeatures)
    names(w1) <- NAME
    umf1 <- unmarkedFramePCount(y= as.matrix(y1_est), siteCovs=z1.1,obsCovs=w1)
    fm1 <- pcount(~w1.1 + w1.2 + w1.3 + w1.4 + w1.5 + w1.6 + w1.7 + w1.8  + w1.9 + w1.10
                  ~x1.1 + x1.2 + x1.3 + x1.4 + x1.5 + x1.6 + x1.7 + x1.8 + x1.9 + x1.10, 
                  umf1, K=nK)
    b1_hat <- coef(fm1, type="det")
    a1_hat <- coef(fm1, type = "state")
    
    prod = matrix(0, nSites, nVisits)
    for (i in 1:nFeatures + 1) {
      prod = prod + W1[[i]] * b1_hat[i]
    }
    p_hat <- sigmoid(prod) # estimated detection probability
    lambda1_hat <- exp(cbind(1, z1) %*% a1_hat)
    
    
    z2.1 <- data.frame(x2=z2)
    NAME <- paste0("w2.", 1:nFeatures)
    names(w2) <- NAME
    umf2 <- unmarkedFramePCount(y=as.matrix(y2_est), siteCovs=z2.1,obsCovs=w2)
    fm2 <- pcount(~w2.1 + w2.2 + w2.3 + w2.4 + w2.5 + w2.6 + w2.7 + w2.8  + w2.9 + w2.10
                  ~x2.1 + x2.2 + x2.3 + x2.4 + x2.5 + x2.6 + x2.7 + x2.8 + x2.9 + x2.10, 
                  umf2, K=nK)
    b2_hat <- coef(fm2, type="det")
    a2_hat <- coef(fm2, type = "state")
    
    prod = matrix(0, nSites, nVisits)
    for (i in 1:nFeatures + 1) {
      prod = prod + W2[[i]] * b2_hat[i]
    }
    fp_hat <- sigmoid(prod) # estimated detection probability
    lambda2_hat <- exp(cbind(1, z2) %*% a2_hat)
    
    # Store results
    mu_est <- p_hat * (matrix(rep(lambda1_hat, nVisits), nSites, nVisits)) + fp_hat * (matrix(rep(lambda2_hat, nVisits), nSites, nVisits))
    mu1_est <- p_hat * (matrix(rep(lambda1_hat, nVisits), nSites, nVisits))
    mu1 <- p * (matrix(rep(lambda1, nVisits), nSites, nVisits))
    mu2_est <- fp_hat * (matrix(rep(lambda2_hat, nVisits), nSites, nVisits))
    mu2 <- fp * (matrix(rep(lambda2, nVisits), nSites, nVisits))
    
    result[k, 1] <- sqrt(mean((p_hat - p)^2)) # USE RMSE
    result[k, 2] <- sqrt(mean(((lambda1_hat - lambda1) / lambda1) ^ 2))  # USE RMSE 
    result[k, 3] <- sqrt(mean((fp_hat - fp)^2))  # USE RMSE
    result[k, 4] <- sqrt(mean(((lambda2_hat - lambda2) / lambda2) ^ 2))  # USE RMSE
    result[k, 5] = mean(abs(A_est - A)) # mean Y1_est - Y1
    result[k, 6] <- sqrt(mean(((mu1_est *  - mu1)/mu1)^2)) # USE RMSE
    result[k, 7] <- sqrt(mean(((mu2_est *  - mu2)/mu2)^2)) # USE RMSE
    # write.table(result, paste("result_case_", case, ".csv", sep=""), col.names=F, row.names=F,sep=",", append=TRUE)    
  #}
    end_time <- Sys.time()
    print(end_time - start_time)
    }
}
print(colMeans(result))
saveRDS(result, paste("/Users/meduse/Downloads/OSU Course Work/2019 Winter/CS 536/Final Project/Code/two_step_nSite_", nSite, "nK_", nK, ".rds", sep=""))
