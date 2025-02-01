# load input
load("E:\\O_Final\\Tengmalm_owl\\results\\input.RData")
library(jagsUI)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X) Integrated Population Model for Tengmalm's owl
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("owl_IPM.jags")
cat("
model {


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Structural equation model to estimate latent 'breeding conditions' (eta)
  # as a function of prey capture rates and breeding phenology
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  zeta[1] <- 1
  zeta[2] ~ dnorm(0,0.01)

  sigma.eta ~ dgamma(1,1)
  tau.eta <- pow(sigma.eta, -2)
  
  sigma.delta ~ dgamma(1,1)
  tau.delta <- pow(sigma.delta, -2)
  
  alpha.delta ~ dnorm(100,0.1)
  alpha.nu ~ dnorm(3, 0.1)


  for (t in 1:(n.occasions)){
  
    # latent breeding conditions
    eta.star[t] ~ dnorm(0, 1)
    eta[t] = eta.star[t] * sigma.eta
    
    # number of voles captured in an average nest
    log(nu[t]) <- alpha.nu + eta[t] * zeta[1]
    
    # mean nest initiation date
    delta[t] <- alpha.delta + eta[t] * zeta[2]
    
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # data likelihoods for prey capture and nest initiation date
  # observed variables
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i in known.d){
    date[i] ~ dnorm(delta[brood.year[i]], tau.delta)
  }
  
  for (i in 1:n.prey){
    x[i] ~ dpois(nu[prey.year[i]])
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # END SEM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Fecundity model
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sigma.xi ~ dgamma(1,1)
  tau.xi <- pow(sigma.xi, -2)
  alpha.xi ~ dnorm(1.5,0.1)

  sigma.psi ~ dgamma(1,1)
  tau.psi <- pow(sigma.psi, -2)
  alpha.psi ~ dnorm(0,0.444)

  beta.xi ~ dnorm(0,0.01)


  beta.psi ~ dnorm(0,0.01)

  sex.ratio <- 0.5

  for (t in 1:n.occasions){
  
    epsilon.star.xi[t] ~ dnorm(0, 1)
    epsilon.xi[t] = epsilon.star.xi[t] * sigma.xi
    log(xi[t]) <- alpha.xi + beta.xi * eta[t] + epsilon.xi[t] 

    epsilon.star.psi[t] ~ dnorm(0, 1)
    epsilon.psi[t] = epsilon.star.psi[t] * sigma.psi
    logit(psi[t]) <- alpha.psi + beta.psi * eta[t] + epsilon.psi[t]    

  }


  
  for (i in 1:n.broods){
    f[i] ~ dpois(psi[brood.year[i]] * xi[brood.year[i]] + 0.0000001)
  }

  for (i in known.c){
    c[i] ~ dpois(xi[brood.year[i]])
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # Population model
  #
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (t in 1:(n.occasions)){
    J[t] ~ dpois(N[t] * xi[t] * psi[t] * sex.ratio)   # fledged juvenile females
    y.nests[t] ~ dpois(N[t])
  }

  # initial population size
  S[1] ~ dpois(15)
  A[1] ~ dpois(15)
  N[1] = S[1] + A[1]

  for (t in 1:(n.occasions-1)){
    S[t+1] ~ dbin(phi.j[t], J[t])                     # fledglings surviving and recruiting
    A[t+1] ~ dbin(phi.a[t], N[t])                     # surviving adults
    N[t+1] = S[t+1] + A[t+1] + I[t]                   # breeding population size
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # immigration
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (t in 1:(n.occasions-1)){
    I[t] ~ dpois(omega[t])
    log(omega[t]) <- alpha.omega + beta.omega * eta[t+1] + kappa.omega * zT[t] + epsilon.omega[t] #
    epsilon.star.omega[t] ~ dnorm(0, 1)
    epsilon.omega[t] = sigma.omega * epsilon.star.omega[t]
  }
  
  kappa.omega ~ dnorm(0,0.1)
  beta.omega ~ dnorm(0,0.1)
  alpha.omega ~ dnorm(0,0.01)
  tau.omega <- pow(sigma.omega,-2)  
  sigma.omega ~ dgamma(1,1)


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # survival and breeding probability
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (t in 1:(n.occasions-1)){
    logit(phi.j[t]) = alpha.phi[1] + beta.phi[1] * (eta[t+1]-eta[t]) + epsilon.phi[t]
    logit(phi.a[t]) = alpha.phi[2] + beta.phi[2] * (eta[t+1]-eta[t]) + epsilon.phi[t]
    logit(p[t,1]) = alpha.p[1] + epsilon.p[t]  
    logit(p[t,2]) = alpha.p[2] + epsilon.p[t]  
    epsilon.star.phi[t] ~ dnorm(0, 1)
    epsilon.phi[t] = epsilon.star.phi[t] * sigma.phi
    
    epsilon.star.p[t] ~ dnorm(0, 1)
    epsilon.p[t] = epsilon.star.p[t] * sigma.p
  }

  beta.phi[1] ~ dnorm(0,0.1)
  beta.phi[2] ~ dnorm(0,0.1)

  alpha.p[1] ~ dnorm(0, 0.444)
  alpha.p[2] ~ dnorm(0, 0.444)
  alpha.phi[1] ~ dnorm(0, 0.444)
  alpha.phi[2] ~ dnorm(0, 0.444)
  
  sigma.p ~ dgamma(1,1)  
  sigma.phi ~ dgamma(1,1)
  # sigma.phi.a ~ dgamma(1,1)
  # tau.phi.a <- pow(sigma.phi.a, -2)
  # sigma.phi.j ~ dgamma(1,1)
  # tau.phi.j <- pow(sigma.phi.j, -2)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # adult m-array
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
   for (t in 1:(n.occasions-1)){
     marr.a[t,1:n.occasions] ~ dmulti(pr.a[t, ], rel.a[t])
   }

   for (t in 1:(n.occasions-1)){
     q.a[t] <- 1 - p[t,2]            
     pr.a[t,t] <- phi.a[t] * p[t,2]
   
   for (j in (t+1):(n.occasions-1)){
     pr.a[t,j] <- prod(phi.a[t:j]) * prod(q.a[t:(j-1)]) * p[j,2]
   }}

   for (t in 1:(n.occasions-1)){
     pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
   }



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # juvenile m-array
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   for (t in 1:(n.occasions-1)){
     marr.j[t,1:n.occasions] ~ dmulti(pr.j[t, ], rel.j[t])
   }

   for (t in 1:(n.occasions-1)){
     pr.j[t,t] <- (phi.j[t] * p[t,1]) * sex.ratio
     q.j[t] <- 1 - p[t,1]
   }
   
   for (t in 1:(n.occasions-2)){
     pr.j[t,t+1] <- (phi.j[t] * q.j[t] * phi.a[t+1] * p[t+1,2]) * sex.ratio
   }

   for (t in 1:(n.occasions-3)){
     for (j in (t+2):(n.occasions-1)){
       pr.j[t,j] <- (phi.j[t] * q.j[t] * prod(phi.a[(t+1):j]) * prod(q.a[(t+1):(j-1)]) * p[j,2]) * sex.ratio
     }
   }
   
   for (t in 1:(n.occasions-1)){
     pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])   # last column
   }
   
   

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # below main diagonals
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   for (t in 2:(n.occasions-1)){
     for (j in 1:(t-1)){
       pr.a[t,j] <- 0
       pr.j[t,j] <- 0
   }}
   
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # goodness-of-fit tests for capture-recapture data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (t in 1:(n.occasions-1)){
    for (j in 1:n.occasions){

      expmarr.j[t,j] <- rel.j[t]*pr.j[t,j]
      expmarr.a[t,j] <- rel.a[t]*pr.a[t,j]

      E.org.j[t,j] <- pow((pow(marr.j[t,j], 0.5)-pow(expmarr.j[t,j], 0.5)), 2)       
      E.org.a[t,j] <- pow((pow(marr.a[t,j], 0.5)-pow(expmarr.a[t,j], 0.5)), 2)
     
      }
   }
   

   for (t in 1:(n.occasions-1)){
     marr.new.j[t,1:n.occasions] ~ dmulti(pr.j[t, ], rel.j[t])  
     marr.new.a[t,1:n.occasions] ~ dmulti(pr.a[t, ], rel.a[t])
     for (j in 1:n.occasions){
       E.new.j[t,j] <- pow((pow(marr.new.j[t,j], 0.5)-pow(expmarr.j[t,j], 0.5)), 2)    
       E.new.a[t,j] <- pow((pow(marr.new.a[t,j], 0.5)-pow(expmarr.a[t,j], 0.5)), 2)
      }
   }
   
   fit[1] <- sum(E.org.j[,])
   fit[2] <- sum(E.org.a[,])   
   fit.new[1] <- sum(E.new.j[,])
   fit.new[2] <- sum(E.new.a[,])   



}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(n.occasions = dim(marr.a)[2], 
                  c = c, f = f, date = date, brood.year = brood.year, n.broods = length(f), known.c = known.c,
                  x = x, territory = territory, prey.year = prey.year, 
                  n.territories = max(territory), n.prey = length(x),
                  age = brood.age-1, vole = vole, mouse = mouse, known.d = known.d,
                  marr.a = marr.a, rel.a = rel.a,
                  marr.j = marr.j, rel.j = rel.j,
                  y.nests = y.nests, zT = zT)

# Initial values
inits <- function(){list(epsilon.star.omega = rep(0, n.years - 1),
                         alpha.p = c(-2,-1), alpha.phi.a = 1, 
                         alpha.phi.j = -2, alpha.omega = 2)}  

# Parameters monitored
parameters <- c('phi.a', 'sigma.phi.a','beta.phi','alpha.phi',
                'phi.j', 'sigma.phi.j','sigma.phi',
                'gamma', 'beta.gamma', 'alpha.gamma','kappa.gamma',
                'omega', 'beta.omega', 'alpha.omega','zeta', 'kappa.omega',
                'beta.phi','alpha.delta','alpha.nu',
                'sigma.omega',
                'xi','alpha.xi','beta.xi','sigma.xi','kappa.xi',
                'psi', 'alpha.psi','beta.psi','sigma.psi','kappa.psi',
                'delta','sigma.delta','beta.p',
                'eta','theta','sigma.eta',
                'mu','alpha.mu',
                'nu','alpha.nu',
                'sigma.d','zeta',
                'p', 'sigma.p','alpha.p',
                'pi','tau',
                'I','B','A','J','S','bS','bA','N',
                'fit','fit.new')


# MCMC settings
ni <- 250000
nt <- 100
nb <- 50000
nc <- 4

# runtime ~1h hours with 4 chains, 250k iterations
Sys.time()
m <- jags(jags.data, inits, parameters, "owl_IPM.jags", 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
          parallel = T, n.cores = nc)
Sys.time()

print(m, digits = 3)
m$Rhat
hist(m$sims.list$sigma.eta)

save.image("E:\\O_Final\\Tengmalm_owl\\results\\joint_SEM_250k_phirandom.RData")










dim(m$sims.list$p)
dim(m$sims.list$phi.j)

library(vioplot)
vioplot(m$sims.list$p[,,1])
vioplot(m$sims.list$p[,,2])

vioplot(m$sims.list$phi.j[,])
vioplot(m$sims.list$phi.a[,])

vioplot(m$sims.list$beta.phi)


hist(m$q50$eta[2:n.years] - m$q50$eta[1:(n.years-1)])


# m$Rhat
# save.image("O:\\Tengmalm_owl\\results\\joint_SEM_no_gamma_500k.RData")
# load("O:\\Tengmalm_owl\\results\\joint_SEM_no_gamma_500k.RData")

dim(m$sims.list$xi)
dim(m$sims.list$nu)

tmp <- m$sims.list$nu/m$sims.list$xi[,,2]
dim(tmp)

vioplot(tmp, drawRect = F)
plot(m$q50$psi[,2] ~ colMeans(tmp))
plot(m$q50$psi[] ~ m$q50$nu)



m$q50$omega
m$q50$N
y.nests

pdf("O:/Tengmalm_owl/JAnE_revision/reviewer_comments/nu_divided_xi.pdf",
    height = 6, width = 6)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(colMeans(tmp) ~ m$q50$eta, 
     ylab = expression(nu/xi), xlab = expression(eta))
dev.off()

par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
plot(m$sims.list$fit.new[,1] ~ m$sims.list$fit[,1], xlim = c(0,50), ylim = c(0,50))
abline(0,1)

plot(m$sims.list$fit.new[,2] ~ m$sims.list$fit[,2], xlim = c(0,50), ylim = c(0,50))
abline(0,1)

mean(m$sims.list$fit.new[,1] > m$sims.list$fit[,1])
mean(m$sims.list$fit.new[,2] > m$sims.list$fit[,2])

hist(m$sims.list$beta.p[,1])
hist(m$sims.list$beta.p[,2])
hist(m$sims.list$sigma.p)


length(m$q50$eta)
plot(m$q50$p[,1] ~ m$q50$eta[1:30])
plot(m$q50$p[,2] ~ m$q50$eta[1:30])
summary(lm(m$q50$p[,2] ~ m$q50$eta[1:30]))

