# Data: Mar.17.2019
# ST 536 Probabilistic Graphical Model Project
# Member: Laurel Hopkins, Eugene Seo, Chuan Tian
# Code description: EM algorithm for N-mixture model

library("unmarked")
library("sigmoid")
nSites <- 200 #500 # number of total sites 
nVisits <- 1 # number of total visits
nFeatures <- 10 # number of features (Assuming it's same for both habitat and detection probability)
nRepeat <- 1 #50
nK <- 300
case <-3
cat("nSites:", nSites, ", nVisits:", nVisits, ", nFeatures:", nFeatures, "\n")
#set.seed("1234")

# a1 <- runif(nFeatures+1, 0, 0.8) # generate random parameters for habitat features (z1)
# a2 <- runif(nFeatures + 1, 0, 0.8) # generate random parameters for habitat features (z2)
# b1 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w1)
# b2 <- rnorm(nFeatures + 1) # generate random parameters for detection prob. features (w2)


#if (file.exists("result_case_1.csv")) {
#  file.remove("result_case_1.csv")
#  file.remove("result_case_2.csv")
#  file.remove("result_case_3.csv")  
#}

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
    result <- array(0, c(1,5))
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
      cat("Fitting models >>> \n")
      
      # Initialization
      z1.1 <- data.frame(x1=z1)
      z2.1 <- data.frame(x2=z2)
      NAME.1 <- paste0("w1.", 1:nFeatures)
      names(w1) <- NAME.1
      NAME.2 <- paste0("w2.", 1:nFeatures)
      names(w2) <- NAME.2
      z.1 <- cbind(z1.1, z2.1)
      w <- c(w1, w2)
      W = c(ones, w)
      umf <- unmarkedFramePCount(y=Y, siteCovs=z.1,obsCovs=w)
      fm <- pcount(~w1.1 + w1.2 + w1.3 + w1.4 + w1.5 + w1.6 + w1.7 + w1.8  + w1.9 + w1.10 +
                     w2.1 + w2.2 + w2.3 + w2.4 + w2.5 + w2.6 + w2.7 + w2.8  + w2.9 + w2.10
                    ~x1.1 + x1.2 + x1.3 + x1.4 + x1.5 + x1.6 + x1.7 + x1.8 + x1.9 + x1.10 + 
                     x2.1 + x2.2 + x2.3 + x2.4 + x2.5 + x2.6 + x2.7 + x2.8 + x2.9 + x2.10,
                    umf, K=nK)
      b_hat <- coef(fm, type="det")
      a_hat <- coef(fm, type = "state")
      
      prod = matrix(0, nSites, nVisits)
      for (i in 1:2 * nFeatures + 1) {
        prod = prod + W[[i]] * b_hat[i]
      }
      p_hat_total <- sigmoid(prod) # estimated detection probability
      z <- cbind(z1, z2)
      lambda_hat <- exp(cbind(1, z) %*% a_hat)
      # result[1] <- sqrt(mean((p_hat - p)^2)) # USE RMSE
      # result[2] <- sqrt(mean(((lambda1_hat - lambda1) / lambda1) ^ 2))  ##FIXED -- Chuan # USE RMSE 
      mu_est <- p_hat_total * (matrix(rep(lambda_hat, nVisits), nSites, nVisits))
      # mu1 <- p * (matrix(rep(lambda1, nVisits), nSites, nVisits))
      # mu2 <- fp * (matrix(rep(lambda2, nVisits), nSites, nVisits))
      mu1_est <- matrix(NA, nSites, nVisits)
      mu2_est <- matrix(NA, nSites, nVisits)
      A_est <- matrix(NA, nSites, nVisits)
      B_est <- matrix(NA, nSites, nVisits)
      
      for (i in 1:nSites) {
        for (j in 1:nVisits) {
          mu1_est[i, j] = runif(n = 1, min = 0, max = mu_est[i, j])
          A_est[i, j] = round(Y[i, j] * mu1_est[i, j] / mu_est[i, j], 0)
          mu2_est[i, j] = mu_est[i, j] - mu1_est[i, j]
          B_est[i, j] = Y[i, j] - A_est[i, j]
        }
      }
      # A_est_0 = matrix(0, nSites, nVisits)
      # B_est_0 = matrix(0, nSites, nVisits)
      # coef_a1 <- rep(0, nFeatures + 1)
      # coef_a2 <- rep(0, nFeatures + 1)
      # coef_b1 <- rep(0, nFeatures + 1)
      # coef_b2 <- rep(0, nFeatures + 1)
      # a1_hat <- rep(5, nFeatures + 1)
      # a2_hat <- rep(5, nFeatures + 1)
      # b1_hat <- rep(5, nFeatures + 1)
      # b2_hat <- rep(5, nFeatures + 1)
      p_hat_0 <- matrix(0, nSites, nVisits)
      fp_hat_0 <- matrix(0, nSites, nVisits)
      lambda1_hat_0 <- matrix(0, nSites, 1)
      lambda2_hat_0 <- matrix(0, nSites, 1)
      p_hat <- matrix(5, nSites, nVisits)
      fp_hat <- matrix(5, nSites, nVisits)
      lambda1_hat <- matrix(5, nSites, 1)
      lambda2_hat <- matrix(5, nSites, 1)
      iter = 0
      
      # E-M 
      # M-step
      while (mean(abs(p_hat_0 - p_hat)) > 0.05 | mean(abs(fp_hat_0 - fp_hat)) > 0.05)
             # |mean(abs(lambda1_hat_0 - lambda1_hat) / lambda1_hat) > 0.5 | mean(abs(lambda2_hat_0 - lambda2_hat) / lambda2_hat) > 0.5)
      { 
        iter = iter + 1
        cat("EM iteration = ", iter, "\n")
        # cat("df_y1 = ", mean(abs(A_est - A_est_0)), "df_y2 = ", mean(abs(B_est - B_est_0)), "\n")
        # df_a1 = mean(abs(coef_a1 - as.vector(a1_hat)))
        # df_a2 = mean(abs(coef_a2 - as.vector(a2_hat)))
        # df_b1 = mean(abs(coef_b1 - as.vector(b1_hat)))
        # df_b2 = mean(abs(coef_b2 - as.vector(b2_hat)))
        # cat("df_a1 = ", df_a1, "df_a2 = ", df_a2, "df_b1 = ", df_b1, "df_b2 = ", df_b2,  "\n")
        # coef_a1 <- as.vector(a1_hat)
        # coef_a2 <- as.vector(a2_hat)
        # coef_b1 <- as.vector(b1_hat)
        # coef_b2 <- as.vector(b2_hat)
        cat("df_p_hat = ", mean(abs(p_hat_0 - p_hat)), "df_fp_hat = ", mean(abs(fp_hat_0 - fp_hat)), "\n")
        cat("df_lambda1_hat = ", mean(abs(lambda1_hat_0 - lambda1_hat) / lambda1_hat), "df_lambda2_hat = ", mean(abs(lambda2_hat_0 - lambda2_hat) / lambda2_hat), "\n")
        p_hat_0 <- p_hat
        fp_hat_0 <- fp_hat
        lambda1_hat_0 <- lambda1_hat
        lambda2_hat_0 <- lambda2_hat
      # A_est_0 = A_est
      # B_est_0 = B_est
      z1.1 <- data.frame(x1=z1)
      NAME <- paste0("w1.", 1:nFeatures)
      names(w1) <- NAME
      umf1 <- unmarkedFramePCount(y=A_est, siteCovs=z1.1,obsCovs=w1)
      fm1 <- pcount(~w1.1 + w1.2 + w1.3 + w1.4 + w1.5 + w1.6 + w1.7 + w1.8  + w1.9 + w1.10
        ~x1.1 + x1.2 + x1.3 + x1.4 + x1.5 + x1.6 + x1.7 + x1.8 + x1.9 + x1.10, 
                    umf1, K=nK)
      b1_hat <- coef(fm1, type="det")
      a1_hat <- coef(fm1, type = "state")
      
      prod = matrix(0, nSites, nVisits)
      for (i in 1:nFeatures + 1) {
        prod = prod + W1[[i]] * b1_hat[i]
      }
      p_hat <- sigmoid(prod) # estimated detection probability  ## FIXED--Chuan
      # N1_hat <- round(bup(ranef(fm1))) # Estimated population size
      lambda1_hat <- exp(cbind(1, z1) %*% a1_hat)
      

      z2.1 <- data.frame(x2=z2)
      NAME <- paste0("w2.", 1:nFeatures)
      names(w2) <- NAME
      # visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)
      umf2 <- unmarkedFramePCount(y=B_est, siteCovs=z2.1,obsCovs=w2)
      fm2 <- pcount(~w2.1 + w2.2 + w2.3 + w2.4 + w2.5 + w2.6 + w2.7 + w2.8  + w2.9 + w2.10
                    ~x2.1 + x2.2 + x2.3 + x2.4 + x2.5 + x2.6 + x2.7 + x2.8 + x2.9 + x2.10, 
                    umf2, K=nK)
      b2_hat <- coef(fm2, type="det")
      a2_hat <- coef(fm2, type = "state")
      
      prod = matrix(0, nSites, nVisits)
      for (i in 1:nFeatures + 1) {
        prod = prod + W2[[i]] * b2_hat[i]
      }
      fp_hat <- sigmoid(prod) # estimated detection probability  ## FIXED--Chuan
      lambda2_hat <- exp(cbind(1, z2) %*% a2_hat)
      
      # E-step
      mu_est <- p_hat * (matrix(rep(lambda1_hat, nVisits), nSites, nVisits)) + fp_hat * (matrix(rep(lambda2_hat, nVisits), nSites, nVisits))
      mu1_est <- p_hat * (matrix(rep(lambda1_hat, nVisits), nSites, nVisits))
      mu2_est <- fp_hat * (matrix(rep(lambda2_hat, nVisits), nSites, nVisits))
      A_est = round(Y * mu1_est / mu_est)
      mu2_est = mu_est - mu1_est
      B_est = Y - A_est
      }
    } 
    cat("final df_p_hat = ", mean(abs(p_hat_0 - p_hat)), "final df_fp_hat = ", mean(abs(fp_hat_0 - fp_hat)), "\n")
    cat("final df_lambda1_hat = ", mean(abs(lambda1_hat_0 - lambda1_hat) / lambda1_hat), 
        "final df_lambda2_hat = ", mean(abs(lambda2_hat_0 - lambda2_hat) / lambda2_hat), "\n")
    result[1] <- result[1] + sqrt(mean((p_hat - p)^2)) # USE RMSE
    result[2] <- sqrt(mean(((lambda1_hat - lambda1) / lambda1) ^ 2))  ##FIXED -- Chuan # USE RMSE 
    result[3] <- sqrt(mean((fp_hat - fp)^2))  # USE RMSE
    result[4] <- sqrt(mean(((lambda2_hat - lambda2) / lambda2) ^ 2))  ## FIXED -- Chuan # USE RMSE
    result[5] = mean(abs(A_est - A)) # mean Y1_est - Y1
    write.table(result, paste("result_case_", case, ".csv", sep=""), col.names=F, row.names=F,sep=",", append=TRUE)    
  #}
    end_time <- Sys.time()
    print(end_time - start_time)  
}
print(colMeans(result))


