# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This script simulates data using parameter estimates from:
#
# Riecke TV, Ravussin P-A, Longchamp L, Trolliet D, Gibson D, & Schaub M
# A hierarchical population model for the estimation of latent prey abundance
# and demographic rates of a nomadic predator
#
# and analyses the data using the data-generating model, adequately recovering
# parameter estimates.
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(IPMbook)
library(latex2exp)
library(jagsUI)
library(vioplot)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulation settings and arrays to store estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nYears <- 25
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of simulations and arrays to store model output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSims <- 250
alpha.est <- array(NA, dim = c(nSims, 7, 7))
beta.est <- array(NA, dim = c(nSims, 5, 7))
zeta.est <- array(NA, dim = c(nSims, 2, 7))
kappa.est <- matrix(NA, nSims, 7)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note that readers can also explore simpler parmeterizations, 
# e.g., simply modeling latent variation in prey abundance measured as
# a function of nest initiation date and prey captures, 
# by deleting portions of the model to simplify the model structure.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("IPM.jags")
cat("
model {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # priors
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # intercepts
  alpha[1] ~ dnorm(3, 1)       # nu: rodent remains
  alpha[2] ~ dnorm(100, 0.1)   # delta: laying date
  alpha[3] ~ dnorm(1, 1)       # xi: clutch size
  alpha[4] ~ dnorm(0, 1)       # psi: egg-to-fledging probability
  alpha[5] ~ dnorm(0, 1)       # phi1: juvenile survival
  alpha[6] ~ dnorm(1, 1)       # phi2: adult survival
  alpha[7] ~ dnorm(1, 1)       # omega: immigrants

  # mean detection and secondary sex ratio (sex ratio at hatch)
  mu.p ~ dbeta(1,1)
  sr = 0.5
  
  # factor loadings for latent variables
  zeta[1] <- 1
  zeta[2] ~ dnorm(0,0.01)
  
  # effects of latent variable on demographic parameters
  beta[1] ~ dnorm(0,0.1)
  beta[2] ~ dnorm(0,0.1)
  beta[3] ~ dnorm(0,0.1)
  beta[4] ~ dnorm(0,0.1)
  beta[5] ~ dnorm(0,0.1)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Structural equation model to estimate latent 'breeding conditions' (eta)
  # as a function of prey capture rates and breeding phenology
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sigma.eta ~ dgamma(1,1)
  sigma.delta ~ dgamma(1,1)
  tau.delta <- pow(sigma.delta, -2)
  
  for (t in 1:(nyears)){
    # latent breeding conditions
    eta.star[t] ~ dnorm(0, 1)
    eta[t] = eta.star[t] * sigma.eta
    # number of voles captured in an average nest
    log(nu[t]) <- alpha[1] + zeta[1] * eta[t]
    # mean nest initiation date
    delta[t] <- alpha[2] + zeta[2] * eta[t]
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # data likelihoods for prey capture and nest initiation date
  # observed variables
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i in 1:nS){
    d[i] ~ dnorm(delta[loop[i]], tau.delta)
    x[i] ~ dpois(nu[loop[i]])    
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # END SEM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # IPM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  for (t in 1:nyears){
    log(xi[t]) <- alpha[3] + beta[1] * eta[t]
    logit(psi[t]) <- alpha[4] + beta[2] * eta[t]   
  }

  for (i in 1:nS){
    n[i] ~ dpois(psi[loop[i]] * xi[loop[i]])  # nestlings
    c[i] ~ dpois(xi[loop[i]])                 # clutch
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Population model
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (t in 1:(nyears)){
    J[t] ~ dpois(N[t] * xi[t] * psi[t] * sr)   # fledged juvenile females
    y[t] ~ dpois(N[t])
  }

  # initial population size
  S[1] ~ dpois(15)
  A[1] ~ dpois(15)
  N[1] = S[1] + A[1]

  for (t in 1:(nyears-1)){
    S[t+1] ~ dbin(phi[t,1], J[t])                     # fledglings surviving and recruiting
    A[t+1] ~ dbin(phi[t,2], N[t])                     # surviving adults
    N[t+1] = S[t+1] + A[t+1] + I[t]                   # breeding population size
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # immigration
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (t in 1:(nyears-1)){
    I[t] ~ dpois(omega[t])
    log(omega[t]) <- alpha[7] + beta[5] * eta[t+1] + kappa * zT[t]
    epsilon.star.omega[t] ~ dnorm(0, 1)
    epsilon.omega[t] = sigma.omega * epsilon.star.omega[t]
  }
  
  kappa ~ dnorm(0,0.1)
  tau.omega <- pow(sigma.omega,-2)  
  sigma.omega ~ dgamma(1,1)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # survival and detection probability
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (t in 1:(nyears-1)){
    logit(phi[t,1]) = alpha[5] + beta[3] * (eta[t+1]-eta[t])
    logit(phi[t,2]) = alpha[6] + beta[4] * (eta[t+1]-eta[t])
    p[t] = mu.p

  }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # state-transition and observation matrices    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # state-transition matrices; t = [1, ..., T-1]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    for (t in 1:(nyears-1)){
      ps[1,t,1] = 0
      ps[1,t,2] = phi[t,1]
      ps[2,t,1] = 0
      ps[2,t,2] = phi[t,2]
    }  

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # state-transition matrices; t = [1, ..., T-1]
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    for (t in 1:(nyears-1)){
      po[1,t] = 0
      po[2,t] = p[t]
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # marray likelihood
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # From here onwards, no changes needed regardless of which model is fitted
    # Calculate probability of non-encounter (dq) and reshape 
    # the array for the encounter probabilities
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (t in 1:(nyears-1)){
      for (s in 1:ns){
        dp[s,t,s] <- po[s,t]
        dq[s,t,s] <- 1-po[s,t]
      } #s
      for (s in 1:(ns-1)){
        for (m in (s+1):ns){
          dp[s,t,m] <- 0
          dq[s,t,m] <- 0
        } #s
      } #m
      for (s in 2:ns){
        for (m in 1:(s-1)){
          dp[s,t,m] <- 0
          dq[s,t,m] <- 0
        } #s
      } #m
    } #t
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define the multinomial likelihood
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (t in 1:((nyears-1)*ns)){
    marr[t,1:(nyears *ns-(ns-1))] ~ dmulti(pi[t,], rel[t])
  }
  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Define cell probabilities of the multistate m-array
    # Matrix U: product of probabilities of state-transition and non-encounter (needed because
    # there is no product function for matrix multiplication in JAGS)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    for (t in 1:(nyears-2)){
      U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- ones
      for (j in (t+1):(nyears-1)){
        U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-2)*ns+(1:ns)] %*% ps[,j-1,] %*% dq[,j-1,]
      } #j
    } #t
    U[(nyears-2)*ns+(1:ns), (nyears-2)*ns+(1:ns)] <- ones
  
  
    for (t in 1:(nyears-2)){
      
      # Diagonal
      pi[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% ps[,t,] %*% dp[,t,]
      
      # Above main diagonal
      for (j in (t+1):(nyears-1)){
        pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] %*% ps[,j,] %*% dp[,j,]
      } #j
      
    } #t
    
    pi[(nyears-2)*ns+(1:ns),(nyears-2)*ns+(1:ns)] <- ps[,nyears-1,] %*% dp[,nyears-1,]
   
    
    # Below main diagonal
    for (t in 2:(nyears-1)){
      for (j in 1:(t-1)){

        pi[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero

      } #j
    } #t
  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Last column: probability of non-recapture
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (t in 1:((nyears-1)*ns)){
      pi[t,(nyears*ns-(ns-1))] <- 1-sum(pi[t,1:((nyears-1)*ns)])
    } #t


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # END IPM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


}
",fill = TRUE)
sink()






for (ii in 1:nSims){
  
  print(ii)
  print(Sys.time())
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # define parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  alpha <- NULL
  zeta <- NULL
  beta <- NULL
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate latent prey abundance (non-AR)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # rho.eta <- 0.5
  eta <- NULL
  # eta[1] <- 0
  # sigma.eta <- 1
  # for (t in 2:n){
  #   eta[t] <- rnorm(1, eta[t-1]*rho.eta, sigma.eta * sqrt(1-rho.eta^2))
  # }
  eta <- rnorm(nYears, 0, 1)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # plot latent prey abundance
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # plot(eta ~ seq(1:nYears), type = 'b', las = 1, cex.lab = 1.5,
  #      ylab = TeX("Latent prey abundance ($\\eta$)"),
  #      xlab = TeX("Year"))
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate population, begin at N = 1
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  S <- NULL
  A <- NULL
  I <- NULL
  J <- NULL
  N <- NULL
  
  zT <- as.numeric(scale(1:nYears))
  
  nu <- NULL
  delta <- NULL
  xi <- NULL
  psi <- NULL
  phi <- matrix(NA, nYears, 2)
  omega <- NULL
  
  # initial population size
  S[1] <- rpois(1, 15)
  A[1] <- rpois(1, 15)
  N[1] = S[1] + A[1]
  
  alpha <- c(3.5, 102, 1.55, -0.3, -1.3, 0.7, 1)
  zeta <- c(1, -16)
  beta <- c(0.3, 1.4, 1.0, 1.356, 0.5)
  kappa <- -0.2
  pi = 0.5
  
  c <- NULL
  d <- NULL
  x <- NULL
  
  # pr <- array(0, c(nYears-1,nYears,2))
  # marr <- array(NA, c(nYears-1,nYears,2))
  z <- matrix(NA, N[1] + S[1], nYears)
  p <- c(0.4,0.4)
  q <- rep(1-p[2], nYears-1)
  
  z[1:S[1],1] <- 2
  z[(S[1]+1):(N[1]+S[1]),1] <- 3
  
  # t <- 1
  for (t in 1:(nYears-1)){
    
    # simulate demographic parameters
    nu[t] <- exp(alpha[1] + zeta[1] * eta[t]) 
    delta[t] <- alpha[2] + zeta[2] * eta[t]
    xi[t] <- exp(alpha[3] + beta[1] * eta[t])
    psi[t] <- plogis(alpha[4] + beta[2] * eta[t])
    phi[t,1] <- plogis(alpha[5] + beta[3] * (eta[t+1] - eta[t]))
    phi[t,2] <- plogis(alpha[6] + beta[4] * (eta[t+1] - eta[t]))
    omega[t] <- exp(alpha[7] + beta[5] * eta[t+1] + kappa * zT[t])
    
    
    # add juveniles (fledglings)
    J[t] <- rpois(1, N[t] * xi[t] * psi[t] * pi)
    zJ <- matrix(rep(c(rep(0, t-1), # before now, the animal was not here
                       1, # now at occasion t, animal is here
                       rep(NA, nYears-t)), # rest of row is NA
                     J[t]), byrow = TRUE, nrow = J[t], ncol = nYears)
    z <- rbind(z, zJ)    


    # survival process for juveniles
    z[which(z[,t] == 1), t+1] <- ifelse(rbinom(J[t], 1, phi[t,1]) == 1, 2, 0) 

    # survival process for second-years (first-time breeders)
    z[which(z[,t] == 2), t+1] <- ifelse(rbinom(length(which(z[,t] == 2)), 1, phi[t,2]) == 1, 3, 0)  
      
    # survival process for adults
    z[which(z[,t] == 3), t+1] <- ifelse(rbinom(length(which(z[,t] == 3)), 1, phi[t,2]) == 1, 3, 0)        
          
    # add immigrants (as second-years, in this case it's irrelevant)
    I[t] <- rpois(1, omega[t])
    zI <- matrix(rep(c(rep(0, t), # before now, the animal was not here
                       2, # now at occasion t, animal is here
                       rep(NA, nYears-t-1)), # rest of row is NA
                     I[t]), byrow = TRUE, nrow = I[t], ncol = nYears)
    z <- rbind(z, zI)  
    
    N[t+1] = length(which(z[,t+1] >= 2))

  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate count and breeding data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  c <- NULL
  n <- NULL
  d <- NULL
  y <- NULL
  tmp <- NULL
  r <- matrix(NA, nYears-1, 2)
  
  for (t in 1:nYears){
    y[t] <- rpois(1, N[t])
  }
  
  for (t in 1:(nYears-1)){
    tmp[t] <- rbinom(1, y[t], p[1])
    
    x <- append(x, rpois(tmp[t], nu[t]))
    d <- append(d, rpois(tmp[t], delta[t]))
    
    c <- append(c, rpois(tmp[t], xi[t]))
    n <- append(n, rpois(tmp[t], xi[t] * psi[t]))
    r[t,1] <- rbinom(1, rpois(1, y[t] * xi[t] * psi[t]), p[1])
    r[t,2] <- rbinom(1, y[t], p[1])
  }
  
  loop <- rep(1:(nYears-1), times = tmp)
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate capture-recapture data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ch <- matrix(NA, nrow(z), ncol(z))
  for (i in 1:nrow(ch)){
    for (t in 1:ncol(ch)){
      
      if (!is.na(z[i,t]) & z[i,t] == 1){
        tmp <- rbinom(1, 1, p[1])
        ch[i,t] <- ifelse(tmp == 1, 1, 0)
      }
      
      # treat SYs and ASYs as equivalent
      if (!is.na(z[i,t]) & z[i,t] >= 2){
        tmp <- rbinom(1, 1, p[2])
        ch[i,t] <- ifelse(tmp == 1, 2, 0)
      }
      
      if (is.na(z[i,t]) | z[i,t] == 0){
        ch[i,t] <- 0
      }
    }
  }
  # table(y)
  ch <- subset(ch, rowSums(ch) > 0)
  marr <- marray(ch)
  rel <- rowSums(marr)

  # Bundle data
  ns <- 2
  jags.data <- list(nyears = nYears, nS = length(loop),
                    c = c, n = n, d = d, x = x, loop = loop, 
                    marr = marr, rel = rel, ns = ns,
                    ones = diag(ns), zero = matrix(0, ns, ns),
                    y = y, zT = zT)
  
  # Initial values
  inits <- function(){list()}  
  
  # Parameters monitored
  parameters <- c('alpha','beta','sigma','p','mu.p','sigma.omega','zeta','kappa',
                  'N','S','A','I','phi','xi','psi','nu','delta','eta','omega')
  
  
  # MCMC settings
  ni <- 25000
  nt <- 10
  nb <- 15000
  nc <- 4
  
  # runtime ~3 minutes for 25k iterations, 4 parallel chains, results may vary
  Sys.time()
  m <- jags(jags.data, inits, parameters, "IPM.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
            parallel = T, n.cores = nc)
  Sys.time()  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # save parameter estimates
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (j in 1:7){
    alpha.est[ii,j,] <- c(m$q50$alpha[j],m$mean$alpha[j],m$q2.5$alpha[j],m$q97.5$alpha[j],
                          m$f$alpha[j],m$n.eff$alpha[j],m$Rhat$alpha[j])
  }
  
  for (j in 1:5){
    beta.est[ii,j,] <- c(m$q50$beta[j],m$mean$beta[j],m$q2.5$beta[j],m$q97.5$beta[j],
                         m$f$beta[j],m$n.eff$beta[j],m$Rhat$beta[j])
  }
  
  for (j in 1:2){
    zeta.est[ii,j,] <- c(m$q50$zeta[j],m$mean$zeta[j],m$q2.5$zeta[j],m$q97.5$zeta[j],
                         m$f$zeta[j],m$n.eff$zeta[j],m$Rhat$zeta[j])
  }
  
  kappa.est[ii,] <- c(m$q50$kappa,m$mean$kappa,m$q2.5$kappa,m$q97.5$kappa,
                      m$f$kappa,m$n.eff$kappa,m$Rhat$kappa)
}







save.image("E:/O_Final/Tengmalm_owl/ECOMOD_Revision_Jan_2025/sim_results_250sims.RData")
load("E:/O_Final/Tengmalm_owl/ECOMOD_Revision_Jan_2025/sim_results_250sims.RData")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# retain converged models as R-hat < 1.2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
keep <- NULL
for (ii in 1:nSims){
  keep[ii] <- ifelse(all(alpha.est[ii,,7] < 1.2) & all(beta.est[ii,,7] < 1.2), 1, 0)
}
table(keep)
keep <- which(keep == 1)
nKeep <- length(keep)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# quantiles of medians for violin plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qalpha <- matrix(NA, 7, 5)
qbeta <- matrix(NA, 5, 5)
for (j in 1:7){
  qalpha[j,] <- quantile(alpha.est[keep,j,1], c(0.025,0.25,0.5,0.75,0.975))
}
for (j in 1:5){
  qbeta[j,] <- quantile(beta.est[keep,j,1], c(0.025,0.25,0.5,0.75,0.975))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a figure of parameter estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf("E:/O_Final/Tengmalm_owl/ECOMOD_Revision_Jan_2025/Figure7.pdf", width = 9, height = 5)
par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))

vioplot(alpha.est[keep,,1], xaxt = 'n', ylim = c(-2,4), drawRect = F)
axis(side = 1, at = c(1:7),
     labels = c(TeX("$\\nu$"), TeX("$\\delta$"),
                TeX("$\\xi$"),TeX("$\\psi$"),TeX("$\\phi_{F}$"),
                TeX("$\\phi_{A}$"),TeX("$\\omega$")))
arrows(1:7, qalpha[,1], 1:7, qalpha[,5], length= 0, col = 'white', lwd = 1)
arrows(1:7, qalpha[,2], 1:7, qalpha[,4], length= 0, col = 'white', lwd = 3)
points(qalpha[,3], cex = 2, pch = 21, bg = 'white')
points(alpha, pch = 1, col = 'red', cex = 2)
mtext(TeX("Intercept ($\\alpha$)"), side = 2, line = 2.5, cex = 2)


vioplot(beta.est[,,1], xaxt = 'n', ylim = c(-2,3.5), drawRect = F)
axis(side = 1, at = c(1:5),
     labels = c(TeX("$\\xi$"),TeX("$\\psi$"),TeX("$\\phi_{F}$"),TeX("$\\phi_{A}$"),TeX("$\\omega$")))
arrows(1:5, qbeta[,1], 1:5, qbeta[,5], length= 0, col = 'white', lwd = 1)
arrows(1:5, qbeta[,2], 1:5, qbeta[,4], length= 0, col = 'white', lwd = 3)
points(qbeta[,3], cex = 2, pch = 21, bg = 'white')
points(beta, pch = 1, col = 'red', cex = 2)
mtext(TeX("Slope ($\\beta$)"), side = 2, line = 2.5, cex = 2)
mtext(TeX("Parameter"), side = 1, line = -1.5, cex = 1.5, outer = T)

dev.off()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter constancy: 
# here this is measured as the mean of the difference between the median
# of the posterior distribution and truth across all simulations that 
# adequately converged
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(alpha.est[keep,1,1] - alpha[1])
mean(alpha.est[keep,2,1] - alpha[2])
mean(alpha.est[keep,3,1] - alpha[3])
mean(alpha.est[keep,4,1] - alpha[4])
mean(alpha.est[keep,5,1] - alpha[5])
mean(alpha.est[keep,6,1] - alpha[6])
mean(alpha.est[keep,7,1] - alpha[7])


mean(beta.est[keep,1,1] - beta[1])
mean(beta.est[keep,2,1] - beta[2])
mean(beta.est[keep,3,1] - beta[3])
mean(beta.est[keep,4,1] - beta[4])
mean(beta.est[keep,5,1] - beta[5])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assess coverage:
# the proportion of simulations that adequately converged in which
# the truth parameter value was contained within the 95% Bayesian credible intervals
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(which(alpha.est[keep,1,3] < alpha[1] & alpha.est[keep,1,4] > alpha[1]))/nKeep
length(which(alpha.est[keep,2,3] < alpha[2] & alpha.est[keep,2,4] > alpha[2]))/nKeep
length(which(alpha.est[keep,3,3] < alpha[3] & alpha.est[keep,3,4] > alpha[3]))/nKeep
length(which(alpha.est[keep,4,3] < alpha[4] & alpha.est[keep,4,4] > alpha[4]))/nKeep
length(which(alpha.est[keep,5,3] < alpha[5] & alpha.est[keep,5,4] > alpha[5]))/nKeep
length(which(alpha.est[keep,6,3] < alpha[6] & alpha.est[keep,6,4] > alpha[6]))/nKeep
length(which(alpha.est[keep,7,3] < alpha[7] & alpha.est[keep,7,4] > alpha[7]))/nKeep


length(which(beta.est[keep,1,3] < beta[1] & beta.est[keep,1,4] > beta[1]))/nKeep
length(which(beta.est[keep,2,3] < beta[2] & beta.est[keep,2,4] > beta[2]))/nKeep
length(which(beta.est[keep,3,3] < beta[3] & beta.est[keep,3,4] > beta[3]))/nKeep
length(which(beta.est[keep,4,3] < beta[4] & beta.est[keep,4,4] > beta[4]))/nKeep
length(which(beta.est[keep,5,3] < beta[5] & beta.est[keep,5,4] > beta[5]))/nKeep






