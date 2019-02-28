# Data: Feb.27.2019
# Author: Eugene Seo
# Description: [ST562] Simulation data for N-mixture model

library("sigmoid")
nSites <- 5 # number of total sites 
nVisits <- 1 # number of total visits
nFeatures <- 10 # number of features (Assuming it's same for both habitat and detection probability)

# Generating simulation data for a species A ==================================
z1 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures) # generate random habitat features
a1 <- runif(nFeatures) # generate random parameters for habitat features (z1)

case = 3 # experiment cases (1: z1 = w1, 2: z1 partially = w1, 3: z1 != w1)
if ( case == 1 ) {
  w1 = z1 # w1 = z1
} else if (case == 2) {
  w1 <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
  sub.idx <- sample(1:nFeatures, round(nFeatures/2))
  w1[,sub.idx] <- z1[,sub.idx] # w1 partially = z1
} else {
  w1 <- matrix(runif(nSites*nFeatures), nSites, nFeatures) # w1 != z1
}
b1 <- rnorm(nFeatures) # generate random parameters for detection prob. features (w1)

lambda1 <- exp(z1 %*% a1) # parameter for the true abundance of A
N1 <- rpois(nSites, lambda1) # true abundance for species A
p <- sigmoid(w1 %*% b1) # detection probability
N1
p

# Generating simulation data for a species B ==================================
z2 <-  matrix(runif(nSites*nFeatures), nSites, nFeatures) # generate random habitat features
a2 <- runif(nFeatures) # generate random parameters for habitat features (z2)

case = 3 # experiment cases (1: z2 = w2, 2: z2 partially = w2, 3: z2 != w2)
if ( case == 1 ) {
  w2 <- z2 # w2 = z2
} else if (case == 2) {
  w2 <- matrix(runif(nSites*nFeatures), nSites, nFeatures) 
  sub.idx <- sample(1:nFeatures, round(nFeatures/2))
  w2[,sub.idx] <- z2[,sub.idx] # w2 partially = z2
} else {
  w2 <- matrix(runif(nSites*nFeatures), nSites, nFeatures) # w2 != z2
}
b2 <- rnorm(nFeatures) # generate random parameters for detection prob. features (w2)

lambda2 <- exp(z2 %*% a2) # parameter for the true abundance of B
N2 <- rpois(nSites, lambda2) # true abundance for species B
fp <- sigmoid(w2 %*% b2) # detection probability
N2
fp


# Generating observed counts Y ================================================
Y <- matrix(NA, nSites, nVisits) # observed counts for A
A <- matrix(NA, nSites, nVisits) # true observed counts from A
B <- matrix(NA, nSites, nVisits) # false observed counts from B (mis-identification)
for(j in 1:nVisits) {
  A[,j] <- rbinom(nSites, N1, p[j])  # true observed counts from A
  B[,j] <- rbinom(nSites, N2, fp[j]) # false observed counts from B (mis-identification)
  Y[,j] <- A[,j] + B[,j] # total observed counts recorded as A
}
A
B
y

# Fit a model with A
visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)
umf1 <- unmarkedFramePCount(y=A, siteCovs=data.frame(x=z1),obsCovs=list(visit=visitMat))
fm1 <- pcount(~w1 ~z1, umf1, K=300)
plogis(coef(fm1, type="det")) # Should be close to p
p
sum(bup(ranef(fm1))) # Estimated population size
sum(N1) # Actual population size

# Fit a model with B
visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)
umf2 <- unmarkedFramePCount(y=B, siteCovs=data.frame(x=z2),obsCovs=list(visit=visitMat))
# Fit a model
fm2 <- pcount(~w2 ~z2, umf2, K=300)
plogis(coef(fm2, type="det")) # Should be close to fp
fp
sum(bup(ranef(fm2))) # Estimated population size
sum(N2) # Actual population size