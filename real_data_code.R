## read data
load("GA_data.Rdata")
load("GA_map.Rdata")
data_use$county <- as.character(data_use$county)
data_use$House_market_price <- data_use$House_market_price/1000
data_use$Population <- data_use$Population/1000
names(data_use) <- c("county", "Housing.Cost", "Unemployment.Rate",
                     "Property.Tax", "House.Market.Price", "White.Race", 
                     "Age", "Population")
library(ggplot2)
GAmap <- fortify(mapdata)

library(dplyr)
long1 <- aggregate(GAmap$long, by = list(subregion = GAmap$subregion), mean)
lat1 <- aggregate(GAmap$lat, by = list(subregion = GAmap$subregion), mean)
GAdata <- data.frame("SR_ID" = c(1:159), "Latitude" = lat1$x, "Longitud" = long1$x)
GAdata <- cbind(GAdata, data_use)

## MGWR-INLA-code
library(INLA)
inla.rgeneric.beta.model <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
           theta = NULL) {
    envir = parent.env(environment())
    
    interpret.theta <- function() { # transform the parameters to make them unbounded
      return(list(prec = exp(theta[1L]),
                  rho = 1/(1+exp(-theta[2L])),
                  mu0 = theta[3L]))
    }
    graph <- function() {
      return (inla.as.sparse(matrix(1, n, n)))
    }
    Q <- function() { # precision matrix
      require(Matrix)
      param <- interpret.theta()
      cov.mat <- 1/param$prec * exp(-D/param$rho)
      prec.mat <- solve(cov.mat)
      return(inla.as.sparse(prec.mat))
    }
    mu <- function() {
      return (rep(interpret.theta()$mu, n))
    }
    log.norm.const <- function() {
      return (numeric())
    }
    log.prior <- function() { # log priors for the parameters of interest
      param = interpret.theta()
      res <- dgamma(param$prec, 0.01, 0.01, log = T) + log(param$prec) +
        log(1) + log(param$rho) + log(1 - param$rho) +
        dnorm(param$mu0, 0, 1, log = T)
      return(res)
    }
    initial <- function() {
      return(c(0, 0, 0))
    }
    quit <- function() {
      return(invisible())
    }
    
    if (!length(theta) || is.null(theta)) theta = initial()
    
    val <- do.call(match.arg(cmd), args = list())
    return(val)
  }

data.inla <- GAdata
# scale the variables to make the bandwidths comparable
data.inla$y <- as.vector(scale(data.inla$Housing.Cost))
data.inla$x1 <- as.vector(scale(data.inla$Unemployment.Rate))
data.inla$x2 <- as.vector(scale(data.inla$Property.Tax))
data.inla$x3 <- as.vector(scale(data.inla$House.Market.Price))
data.inla$x4 <- as.vector(scale(data.inla$White.Race))
data.inla$x5 <- as.vector(scale(data.inla$Age))
data.inla$x6 <- as.vector(scale(data.inla$Population))

Dist0 <- as.matrix(dist(data.frame(data.inla$Longitud, data.inla$Latitude)))
Dmax <- max(Dist0); Dmax
# scale the distance matrix with the maximum distance
Dist01 <- Dist0/Dmax
D <- Dist01
NN <- nrow(D); NN

# generate the self-defined function
beta.model <- inla.rgeneric.define(inla.rgeneric.beta.model,
                                   D = D,
                                   n = NN)
# model formula
f.beta <- y ~ -1 + 
  f(idx1, x1, model = beta.model, n = NN) + 
  f(idx2, x2, model = beta.model, n = NN) + f(idx3, x3, model = beta.model, n = NN) + 
  f(idx4, x4, model = beta.model, n = NN) + f(idx5, x5, model = beta.model, n = NN) +
  f(idx6, x6, model = beta.model, n = NN) 
# run the model to get the results
m.beta <- inla(f.beta,
               safe = FALSE, 
               data = data.frame(data.inla, idx1 = 1:NN, idx2 = 1:NN, idx3 = 1:NN
                                 , idx4 = 1:NN, idx5 = 1:NN, idx6 = 1:NN),
               family = "gaussian",
               verbose = T)
