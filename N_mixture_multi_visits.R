# Data: Feb.28.2019
# ST 562 Probabilistic Graphical Model Project
# Member: Laurel Hopkins, Eugene Seo, Chuan Tian
# Code description: Simulation data for N-mixture model

remove(list=ls())
library("unmarked")
library("sigmoid")

run_sim <- function(nSites, nVisits) {
  # fixed parameters
  nFeatures <- 10 # number of features (Assuming it's same for both habitat and detection probability)
  nRepeat <- 50

  # experimental parameters
  #nSites <- 200 # number of total sites (200~500)
  #nVisits <- 5 # number of total visits (1~5)  

  # clear the previous result files
  #if (file.exists(paste("s", nSites, "_v", nVisits, "_k", nK, "_c1.csv"))) {
  #  for (i in 1:3) {
  #    file.remove(paste("s", nSites, "_v", nVisits, "_k", nK, "_c", i, ".csv", sep=""))
  #  }
  #}

  for (step in 1:nRepeat) {
    start_time <- Sys.time()
    cat("Repeat", step, "nSites:", nSites, "nVisits:", nVisits, "\n")
    
    # 1. generate parameters for habitat features
    a1 <- runif(nFeatures + 1) 
    a2 <- runif(nFeatures + 1) 
    
    # 2. generate parameters for detection features
    b1 <- rnorm(nFeatures + 1) 
    b2 <- rnorm(nFeatures + 1)
     
    # 3. generate habitat features
    z1 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures)
    z2 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures)
    
    # 4. generate detection features according to the degree of overlap with habitat features
    if (nVisits == 1) init_case = 1 else init_case = 3
    for (case in init_case:3) {      
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
      nK <- max(Y) + 100 # max number of observations
      result[5] <- max(Y)
      result[6] <- mean(Y)    
      
      # 8. learning models
      #if ( FALSE ) {
      if ( max(Y) < 1000 ) {
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
        
        write.table(result, paste("s", nSites, "_v", nVisits, "_k", 0, "_c", case, ".csv", sep=""), col.names=F, row.names=F,sep=",", append=TRUE)    
      }
    }
    end_time <- Sys.time()
    print(end_time - start_time)  
  }
  print(colMeans(result))
}


merge_results <- function() {
  nSites_list = c(200,300,400,500)
  nVisits_list = c(1,5,10)
  nK_list = c(0) #,100,500,1000)
  case_list = c(1,2,3)
  combined_result <- data.frame()
  for (s in nSites_list) {
    for (v in nVisits_list) {
      for (k in nK_list) {
        for (c in case_list) {
          filename = paste("csv/s", s, "_v", v, "_k", k, "_c", c, ".csv", sep="")          
          if (file.exists(filename)) {
            print(filename)
            myResult = read.csv(filename, header=FALSE)
            colnames(myResult) = c("RMSE_p", "RMSE_lambda1", "RMSE_fp", "RMSE_lambda2", "maxY", "meanY")
            Case <- rep(c,nrow(myResult))
            myResult <- cbind(Case, myResult)
            nK <- rep(k,nrow(myResult))
            myResult <- cbind(nK, myResult)
            nVisits <- rep(v,nrow(myResult))
            myResult <- cbind(nVisits, myResult)
            nSites <- rep(s,nrow(myResult))
            myResult <- cbind(nSites, myResult)
            combined_result <- rbind(combined_result, myResult)
            print(dim(combined_result))
          }
        }
      } 
    }
  }
  return(combined_result)
}

library(FSA)
library(ggplot2)
library(gridExtra)
draw_plot <- function() {  
  result = merge_results()
  print(max(result$maxY))
  print(mean(result$maxY))
  result$YmaxGroup <- cut(result$maxY, breaks=c(0,50,100,150,200))
  result$YmeanGroup <- cut(result$meanY, breaks=c(0,15,30,45,60))
  
  group_list = unique(result$YmaxGroup)
  for (i in group_list) {
    sub_result = result[result$nVisits == 1 & result$YmaxGroup == i,] 
    sub_result$Case = factor(sub_result$Case)
    # nSites vs RMSE (nVisits = 1)
    Sum <- Summarize(RMSE_p ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p1 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_p") 
    
    Sum <- Summarize(RMSE_lambda1 ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p2 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda1")
    
    Sum <- Summarize(RMSE_fp ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p3 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_fp") 
    
    Sum <- Summarize(RMSE_lambda2 ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p4 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda2")
    
    png(paste("figures/nSites-RMSE (maxY-",i,").png", sep=""))
    grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, respect = T, top = paste("nVisits = 1, Range of max Y = ", i, sep=""))
    dev.off()
    
    
    # nVisits vs RMSE (Case = 3)
    sub_result = result[result$Case == 3 & result$YmaxGroup == i,] 
    sub_result$nSites = factor(sub_result$nSites)
    Sum <- Summarize(RMSE_p ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p1 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_p") 
    
    Sum <- Summarize(RMSE_lambda1 ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p2 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda1")
    
    Sum <- Summarize(RMSE_fp ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p3 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_fp") 
    
    Sum <- Summarize(RMSE_lambda2 ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p4 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda2")
    
    png(paste("figures/nVisits-RMSE (maxY-",i,").png", sep=""))
    grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, respect = T, top = paste("Case = 3, Range of max Y = ", i, sep=""))
    dev.off()  
  }
  
  group_list = unique(result$YmeanGroup)
  for (i in group_list) {
    sub_result = result[result$nVisits == 1 & result$YmeanGroup == i,] 
    sub_result$Case = factor(sub_result$Case)
    # nSites vs RMSE (nVisits = 1)
    Sum <- Summarize(RMSE_p ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p1 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_p") 
    
    Sum <- Summarize(RMSE_lambda1 ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p2 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda1")
    
    Sum <- Summarize(RMSE_fp ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p3 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_fp") 
    
    Sum <- Summarize(RMSE_lambda2 ~ nSites + Case, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p4 <- ggplot(Sum, aes(x=nSites, y=mean, group=Case, color=Case)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda2")
    
    png(paste("figures/nSites-RMSE (meanY-",i,").png", sep=""))
    grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, respect = T, top = paste("nVisits = 1, Range of mean Y = ", i, sep=""))
    dev.off()
    
    # nVisits vs RMSE (Case = 3)
    sub_result = result[result$Case == 3 & result$YmeanGroup == i,] 
    sub_result$nSites = factor(sub_result$nSites)
    Sum <- Summarize(RMSE_p ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p1 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_p") 
    
    Sum <- Summarize(RMSE_lambda1 ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p2 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda1")
    
    Sum <- Summarize(RMSE_fp ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p3 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_fp") 
    
    Sum <- Summarize(RMSE_lambda2 ~ nVisits + nSites, data=sub_result)
    Sum$se <- Sum$sd / sqrt(Sum$n)
    p4 <- ggplot(Sum, aes(x=nVisits, y=mean, group=nSites, color=nSites)) + geom_line()+
    geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda2")
    
    png(paste("figures/nVisits-RMSE (meanY-",i,").png", sep=""))
    grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, respect = T, top = paste("Case = 3, Range of mean Y = ", i, sep=""))
    dev.off()  
  }  
  
  # MaxY vs RMSE (nSites = 500, nVisits = 1)
  result$YmaxGroup <- cut(result$maxY, breaks=c(0,50,100,150,200))
  sub_result = result[result$nSites == 200 & result$nVisits == 1,]
  level_order <- c('(0,50]', '(50,100]', '(100,150]', '(150,200]')
  
  Sum <- Summarize(RMSE_p ~ YmaxGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p1 <- ggplot(Sum, aes(x=factor(YmaxGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_p") + xlab("max Y")
  
  Sum <- Summarize(RMSE_lambda1 ~ YmaxGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p2 <- ggplot(Sum, aes(x=factor(YmaxGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda1") + xlab("max Y")
  
  Sum <- Summarize(RMSE_fp ~ YmaxGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p3 <- ggplot(Sum, aes(x=factor(YmaxGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_fp") + xlab("max Y")
  
  Sum <- Summarize(RMSE_lambda2 ~ YmaxGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p4 <- ggplot(Sum, aes(x=factor(YmaxGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda2") + xlab("max Y")

  png("figures/maxY-RMSE.png")
  grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, respect = T, top = "nSites = 500 & nVisits == 1")
  dev.off()
  
  # meanY vs RMSE (nSites = 500, nVisits = 1)
  result$YmeanGroup <- cut(result$meanY, breaks=c(0,15,30,45,60))
  sub_result = result[result$nSites == 200 & result$nVisits == 1,]
  level_order <- c('(0,15]', '(15,30]', '(30,45]', '(45,60]')
  
  Sum <- Summarize(RMSE_p ~ YmeanGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p1 <- ggplot(Sum, aes(x=factor(YmeanGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_p") + xlab("mean Y")
  
  Sum <- Summarize(RMSE_lambda1 ~ YmeanGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p2 <- ggplot(Sum, aes(x=factor(YmeanGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda1") + xlab("mean Y")
  
  Sum <- Summarize(RMSE_fp ~ YmeanGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p3 <- ggplot(Sum, aes(x=factor(YmeanGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_fp") + xlab("mean Y")
  
  Sum <- Summarize(RMSE_lambda2 ~ YmeanGroup, data=sub_result)
  Sum$se <- Sum$sd / sqrt(Sum$n)
  p4 <- ggplot(Sum, aes(x=factor(YmeanGroup, level=level_order), y=mean, group = 1)) + geom_line()+
  geom_pointrange(aes(ymin=mean-se, ymax=mean+se)) + ylab("RMSE_lambda2") + xlab("mean Y")

  png("figures/meanY-RMSE.png")
  grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2, respect = T, top = "nSites = 500 & nVisits == 1")
  dev.off()    
  
  #ggplot(Sum, aes(x = nSites, y = mean, color = Case)) +
  #  geom_point(shape = 15, size  = 4) +
  #  geom_errorbar(aes(ymin  = mean - se, ymax  = mean + se),
  #                    width = 0.2, size  = 0.7) +
  #  theme_bw() + theme(axis.title = element_text(face = "bold")) +
  #  ylab("RMSE on p") 
  
  
  #library(dplyr)
  #sub_result = combined_result[combined_result$nVisits == 1,]
  #sub_summary <- 
  #  sub_result %>% 
  #  group_by(nSites, Case) %>%
  #  summarise(RMSE_p = mean(RMSE_p, na.rm = TRUE), RMSE_lambda1 = mean(RMSE_lambda1, na.rm = TRUE), 
  #            RMSE_fp = mean(RMSE_fp, na.rm = TRUE), RMSE_lambda2 = mean(RMSE_lambda2, na.rm = TRUE))
  #sub_summary = data.frame(sub_summary)  
  #sub_summary$Case = factor(sub_summary$Case)
  #ggplot(data=sub_summary, aes(x=nSites, y=RMSE_p, group = Case, colour = Case)) +
  #  geom_line()
}