n <- 50 #number of samples

for (iter in 1:200) {
  result <- array(0, c(1,10))
  # Generate Y1
  k1 = 3   #number of features
  # regression coefficients
  betas_1 = runif(n = k1 + 1, min = 0.3, max = 1.3)
  # generate covariate values
  X1 <- cbind(rep(1, n), matrix(runif(n=n * k1, min= -1, max=1.5), nrow = n))
  # compute mu's
  mu1 <- as.vector(exp(X1 %*% betas_1))
  # generate y1's
  y1 <- rpois(n=n, lambda=mu1)
  # data set
  my_data_1 <- data.frame(y1 = y1, x1_1 = X1[,2], x1_2 = X1[,3], x1_3 = X1[,4])

  # Generate Y2
  k2 = 4   #number of features
  # regression coefficients
  betas_2 = runif(n = k2 + 1, min = 0.3, max = 1.3)
  # generate covariate values
  X2 <- cbind(rep(1, n), matrix(runif(n=n * k2, min= -1, max=1.5), nrow = n))
  # compute mu's
  mu2 <- as.vector(exp(X2 %*% betas_2))
  # generate y2's
  y2 <- rpois(n=n, lambda=mu2)
  # data set
  my_data_2 <- data.frame(y2 = y2, 
                          x2_1 = X2[,2], x2_2 = X2[,3], x2_3 = X2[,4], x2_4 = X2[,5])
  my_data = cbind(my_data_1, my_data_2)

  model_total <- glm(y1 + y2 ~ x1_1 + x1_2 + x1_3 + x2_1 + x2_2 + x2_3 + x2_4, family = poisson(),
                     data = my_data)

  # Initialization
  y = my_data$y1 + my_data$y2
  cat(iter, max(y), "\n")
  mu_est <- as.vector(model_total$fitted)
  mu1_est <- rep(NA, n)
  mu2_est = rep(NA, n)
  y1_est = rep(NA, n)
  y2_est = rep(NA, n)

  for (i in 1:n) {
    mu1_est[i] = runif(n = 1, min = 0, max = mu_est[i])
    y1_est[i] = round(y[i] * mu1_est[i] / mu_est[i], 0)
    mu2_est[i] = mu_est[i] - mu1_est[i]
    y2_est[i] = y[i] - y1_est[i]
  }

  my_data_1 <- cbind(y1_est, my_data_1)
  my_data_2 <- cbind(y2_est, my_data_2)
  model_1 <- glm(y1_est ~ x1_1 + x1_2 + x1_3, family = poisson(), data = my_data_1)
  model_2 <- glm(y2_est ~ x2_1 + x2_2 + x2_3 + x2_4, family = poisson(), data = my_data_2)
  mu1_est <- as.vector(model_1$fitted)
  mu2_est <- as.vector(model_2$fitted)
  mu_est = mu1_est + mu2_est

  coef_1 <- rep(0, k1 + 1)
  coef_2 <- rep(0, k2 + 1)

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
    model_1 <- glm(y1_est ~ x1_1 + x1_2 + x1_3, family = poisson(), data = my_data_1)
    model_2 <- glm(y2_est ~ x2_1 + x2_2 + x2_3 + x2_4, family = poisson(), data = my_data_2)
    mu1_est <- as.vector(model_1$fitted)
    mu2_est <- as.vector(model_2$fitted)
    mu_est = mu1_est + mu2_est
  }
  result[1] = sqrt(mean((coef_1 - betas_1)^2)) 
  result[2] = sqrt(mean((coef_2 - betas_2)^2)) 
  result[3] = sqrt(mean(((mu1_est - mu1) / mu1) ^ 2)) 
  result[4] = sqrt(mean(((mu2_est - mu2) / mu2) ^ 2))
  result[5] = sqrt(mean((y1_est - y1)^2)) 
  result[6] = sqrt(mean((y2_est - y2)^2)) 
  result[7] = max(y1)
  result[8] = max(y2)
  result[9] = mean(y1)
  result[10] = mean(y2)
  write.table(result, "EM/EM_results.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  
  write.table(t(mu1), "EM/EM_mu1.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(mu2), "EM/EM_mu2.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(mu1_est), "EM/EM_mu1_est.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(mu2_est), "EM/EM_mu2_est.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  
  write.table(t(y1), "EM/EM_y1.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(y2), "EM/EM_y2.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(y1_est), "EM/EM_y1_est.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(y2_est), "EM/EM_y2_est.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  
  write.table(t(betas_1), "EM/EM_betas_1.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(betas_2), "EM/EM_betas_2.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(coef_1), "EM/EM_coef_1.csv", col.names=F, row.names=F, sep=",", append=TRUE)
  write.table(t(coef_2), "EM/EM_coef_2.csv", col.names=F, row.names=F, sep=",", append=TRUE)
}
