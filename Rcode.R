############################################################
# Additional Robustness Simulations
# dCov vs Tau* vs HSIC
############################################################

library(energy)
library(TauStar)
library(dHSIC)

set.seed(12345)

############################################################
# One simulation replicate
############################################################

one.rep <- function(n = 100,
                    scenario = "heavy_tail",
                    alpha = 0.05){

  X1 <- runif(n,-1,1)
  X2 <- runif(n,-1,1)

  ##########################################################
  # Error generation
  ##########################################################

  if(scenario=="heavy_tail"){

    eps <- rt(n,df=3)

  }

  if(scenario=="laplace"){

    U <- runif(n)
    eps <- ifelse(U<0.5,
                  log(2*U),
                  -log(2*(1-U)))

  }

  if(scenario=="skewed"){

    eps <- rchisq(n,df=3)-3

  }

  if(scenario=="heteroscedastic"){

    eta <- rnorm(n)
    eps <- (1+0.5*abs(X1))*eta

  }

  if(scenario=="contaminated"){

    ind <- rbinom(n,1,0.05)

    eps <- (1-ind)*rnorm(n,0,1) +
      ind*rnorm(n,0,5)

  }

  ##########################################################
  # Alternatives
  ##########################################################

  if(scenario=="oscillatory"){

    eps <- 0.3*sin(2*pi*X1)+rnorm(n)

  }

  if(scenario=="quadratic"){

    eps <- 0.3*(X1^2-1)+rnorm(n)

  }

  ##########################################################
  # Covariates
  ##########################################################

  W <- cbind(X1,X2)

  ##########################################################
  # dCov
  ##########################################################

  dcov.p <- dcov.test(eps,W,R=199)$p.value

  ##########################################################
  # Tau*
  ##########################################################

  tau.p <- TauStar::taustar.test(X1,eps)$p.value

  ##########################################################
  # HSIC
  ##########################################################

  hsic.p <- dhsic.test(X=W,
                       Y=matrix(eps,ncol=1),
                       method="permutation",
                       B=199)$p.value

  c(dcov.p < alpha,
    tau.p < alpha,
    hsic.p < alpha)

}

############################################################
# Empirical power
############################################################

power.study <- function(B=500,
                        n=100,
                        scenario){

  out <- replicate(B,
                   one.rep(n=n,
                           scenario=scenario))

  rowMeans(out)

}

############################################################
# Scenarios
############################################################

scenarios <- c(
  "heavy_tail",
  "laplace",
  "skewed",
  "heteroscedastic",
  "contaminated",
  "oscillatory",
  "quadratic"
)

############################################################
# Run simulation
############################################################

results <- matrix(NA,
                  nrow=length(scenarios),
                  ncol=3)

rownames(results) <- scenarios

colnames(results) <- c(
  "dCov",
  "TauStar",
  "HSIC"
)

for(i in 1:length(scenarios)){

  cat("Running:",
      scenarios[i],"\n")

  results[i,] <-
    power.study(
      B=500,
      n=100,
      scenario=scenarios[i]
    )

}

############################################################
# Table
############################################################

results <- round(results,3)

results
