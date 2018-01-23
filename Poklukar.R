## ----copyright, echo=FALSE, results='hide', message=FALSE----------------
## ----------------------------------------------------------- ##
## COPYRIGHT:
##
## Author: Petra Poklukar, Stockholm University
## Course: Bayesian Methods HT2017 at the Department of Mathematics
##         University of Stockholm (course page available on
##         http://kurser.math.su.se/course/view.php?id=586)
##
## Description: R file corresponding to the Hand in exercise 
##         available at http://kurser.math.su.se/pluginfile.php/48029/mod_resource  
##         /content/4/lab3.pdf. This file contains computational
##         solutions to Exercises 1, 3 and 4.
##
## Date: 01-11-2017
## ----------------------------------------------------------- ##

## ----setup, include=FALSE, cache=FALSE-----------------------------------
## ----------------------------------------------------------- ##
## LIBRARIES:

# - knitr
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', fig.show='hold', 
               size='footnotesize')
options(formatR.arrow=TRUE,width=90)

# - rjags
library("rjags")
set.seed(12345)
## ----------------------------------------------------------- ##

## ----ex1.parameters, echo=FALSE------------------------------------------
# ------------------------------------------------------------- #
# EXERCISE 1
# Comparing the results of Bayesian inference for different 
# sampling schemes and priors.
# ------------------------------------------------------------- #
#
# Plot parameters:
# xa - the number of defected elements
# xb - the number of perfect products preceding the 3rd defect
# rb - the number of observed defects in part b)
# theta - the probability of a defected product
# ------------------------------------------------------------- #
xa <- 3
xb <- 27 
rb <- 3
theta <- seq(0, 1, length.out=1000)

# a) 
# Likelihood ~ Binomial(30, theta)
# Prior - uniform, Pi(theta) = 1 
# Posterior ~ Beta(4, 28)

## ----Ex1plots, cache=TRUE, echo=FALSE, messages=FALSE, fig.cap="Plots of the four posterior distributions from parts \\eqref{a}, \\eqref{b}, \\eqref{c} and \\eqref{dd}.", fig.height=6----
plot(theta, dbeta(theta, xa+1, 31-xa), type="l", col="red", lwd = 2, 
     xlab=expression(theta), ylab="Posterior", ylim=c(0, 8.5), xlim=c(0, 0.5))

# b) 
# Likelihood ~ NegBin(rb, theta)
# Prior - uniform, Pi(theta) = 1 
# Posterior ~ Beta(4, 28) 

lines(theta, dbeta(theta, rb+1, xb+1), type="l", col="firebrick2", lwd = 2) 

# c) 
# Likelihood ~ Binomial(30, theta)
# Prior - Jeffrey's, Pi(theta) ~ Beta(0.5, 0.5)
# Posterior ~ Beta(3.5, 27.5) 

lines(theta, dbeta(theta, xa+1/2, 61/2-xa), type="l", col="darkorange1", lwd = 2) 

# d) 
# Likelihood ~ NegBin(rb, theta)
# Prior - Jeffrey's, Pi(theta) ~ Beta(0, 0.5)
# Posterior ~ Beta(3, 27.5) 

lines(theta, dbeta(theta, rb, xb+1/2), type="l", col="dodgerblue3", lwd = 2)

legend("topright",  inset=.05, legend=c("a,b) Posterior Beta(4, 28)\nusing uniform prior", "c) Posterior Beta(3.5, 27.5)\nusing Jeffrey's prior", "d) Posterior Beta(3, 27.5)\nusing Jeffrey's prior"),
       col=c("firebrick2", "darkorange1", "dodgerblue3"), lwd=2, bty="n",
       cex=1.1, pt.cex = 1, y.intersp = 1.8)

## ----EXERCISE3, echo=FALSE-----------------------------------------------
# ------------------------------------------------------------- #
# EXERCISE 3
# Empirical Bayes inference in the Poisson-Gamma model.
# ------------------------------------------------------------- #

## ----mloglikFun, echo=TRUE-----------------------------------------------
# ------------------------------------------------------------- #
# mloglik: returns the log marginal likelihood for the 
# Poisson-Gamma model
# 
# Parameters:
# eta - vector of hyperparameters (a, b) for prior Gamma distribution
# y - n-dimensional vector of observations 
# t - n-dimensional correction vector in Poisson distribution 
# ------------------------------------------------------------- #
mloglik <- function(eta, y, t) {
  sum(dnbinom(y, size = eta[1], prob = 1/(eta[2] * t + 1), log = TRUE))
}

## ----psummaryFun, tidy=TRUE, highlight=TRUE------------------------------
# ------------------------------------------------------------- #
# psummary: for given Poisson-Gamma model with parameters 'lambda', 
# 'e', 'a' and 'b', this function returns a matrix of dimensions 
# n x 4 where n is the number of observations. For each pair (yi, ei), 
# the ith row consists of posterior mean of lambda, an equitailed 
# 100x(1-alpha)% credibility region for lambda, and the posteriror
# probability of the hypothesis H0: lambda <= 1
# 
# Parameters:
# y - vector of observations (n-dimensional)
# t - correction vector in Poisson distribution (n-dimensional)
# a - shape parameter of the prior Gamma distribution
# b - scale parameter of the prior Gamma distribution
# alpha - credibility region parameter
# ------------------------------------------------------------- #
psummary <- function(y, t, a, b, alpha = 0.05) {
    n <- length(y)
    res <- matrix(ncol = 4, nrow = n, 
                  dimnames=list(NULL, c("posterior mean ", 
                                        "2.5% tail ", 
                                        "97.5% tail ", 
                                        "H_0: lambda <= 1")))

    for (i in 1:n) {
        # posterior mean is the weighted average
        res[i, 1] <- (a + y[i]) * (b/(b * t[i] + 1))
        
        # 95% credibility interval2
        res[i, 2:3] <- qgamma(c(alpha/2, 1 - alpha/2), 
                              shape = y[i] + a, 
                              scale = b/(b * t[i] + 1))
        
        # probability that lambda <= 1
        res[i, 4] <- pgamma(1, shape = y[i] + a, 
                            scale = b/(b * t[i] + 1))
    }
    return(res)
}

## ----psummary.altFun, message=F, warning=FALSE---------------------------
# ------------------------------------------------------------- #
# psummary.alt: returns exactly the same as psummary. The 
# only difference is that it uses JAGS library for
# computations instead of exact formulas from parts a, b, c.
#
# Parameters: see psummary description
# ------------------------------------------------------------- #
psummary.alt <- function(y, t, a, b, alpha = 0.05) {
    n <- length(y)
    res <- matrix(ncol = 4, nrow = n, 
                  dimnames=list(NULL, c("posterior mean ", 
                                        "2.5% tail ", 
                                        "97.5% tail ",
                                        "H0: lambda <= 1")))
    
    for (i in 1:n) {
        res[i, ] <- createJagsModel(y[i], t[i], a, b, alpha = 0.05)
    }
    res
}

# ------------------------------------------------------------- #
# createJagsModel: compiles a jags model and returns a list 
# consisting of the posterior mean of lambda, an equitailed 
# 100x(1-alpha)% credibility region for lambda, and the 
# posterior probability of the hypothesis H0: lambda <= 1. 
# 
# Parameters:
# yi- the ith observation
# ti - the ith correction scalar
# a - shape parameter of the prior Gamma distribution
# b - scale parameter of the prior Gamma distribution
# alpha - credibility region parameter
# ------------------------------------------------------------- #
createJagsModel <- function(yi, ti, a, b, alpha = 0.05) {
    PLmodel.string <- "model{
        mu <- ti*lambda
        yi ~ dpois(mu)
    
        # Prior
        lambda ~ dgamma(a, 1/b)
    }"
    
    PLmodel <- jags.model(textConnection(PLmodel.string),
                          data = list('yi' = yi, 'ti' = ti,
                                      'a' = a, 'b' = b), 
                          n.chains = 1, n.adapt = 100, quiet = TRUE)
    
    samples.i <- window(coda.samples(PLmodel, 
                                     variable.names = c("lambda"), 
                                     n.iter = 20000), start = 9000)
    lambda <- as.numeric(unlist(samples.i))
    
    postMean.i <- mean(lambda)
    ci95.i <- quantile(lambda, c(alpha/2, 1 - alpha/2))
    h0.i <- mean(lambda <= 1)
    
    return(c("posterior mean" = postMean.i, ci95.i, h0.i))
}

## ----max.lik.optim, cache=TRUE-------------------------------------------
# optimizing  marginal log likelihood with given values
y <- c(13, 5, 36)
e <- c(5.7219, 8.9395, 40.8851) # we use e instead of t
max.lik <- optim(c(1, 1), mloglik, y = y, t = e, 
                 control = list(fnscale = -1))$par

## ----psummaryResults, eval=TRUE, cache=TRUE------------------------------
# output of psummary with optimized parameters
psummary.res <- psummary(y, e, max.lik[1], max.lik[2])
psummary.alt.res <- psummary.alt(y, e, max.lik[1], max.lik[2])

## ----table.summary, echo=FALSE, eval=TRUE, cache=TRUE--------------------
# displaying results of psummary
library(kableExtra)
options(knitr.table.format = "latex")
kable(psummary.res, booktabs = T,
      caption = "Output of the \\texttt{psummary} function.", align = 'c') %>%
kable_styling(position = "center") %>%
row_spec(0, bold = T)

## ----table.summary.alt, echo=FALSE, eval=TRUE, cache=TRUE----------------
# displaying results of psummary.alt
kable(psummary.alt.res, booktabs = T,
      caption = "Output of the \\texttt{psummary.alt} function.", align = 'c') %>%
kable_styling(position = "center") %>%
row_spec(0, bold = T)

## ----pair1.densities, echo=FALSE, fig.cap="Plots of  prior, likelihood and posterior densities as functions of $\\lambda$ for the pair $(y_1, t_1)$.", cache=TRUE, fig.height=6----
# plotting prior, likelihood and posterior for pair (y1, t1)
a <- max.lik[1]
b <- max.lik[2]
i <- 1
lambda.grid <- seq(0, 5, length.out = 1000)

# prior
plot(lambda.grid, dgamma(lambda.grid, shape = a, scale = b), 
     type = "l", lwd = 2, col = "darkorange", xlab = expression(lambda), 
     ylab = "density", lty = 3, ylim = c(0, 1.1))

# posterior
lines(lambda.grid, dgamma(lambda.grid, shape = a + y[i], 
                          scale = b/(b * e[i] + 1)), 
      type = "l", lwd = 2, col = "darkred", xlab = expression(lambda))

# likelihood
lines(lambda.grid, dpois(y[i], e[i] * lambda.grid), type = "l", 
      lwd = 2, col = "darkblue", lty = 2, xlab = expression(lambda))

legend("topright", inset = .02, legend = c("Prior density", 
                                           "Likelihood density", 
                                           "Posterior density"),
       col = c("darkorange", "darkblue", "darkred"), lty = c(3, 2, 1), 
       lwd = 2, title = paste("For pair (y1, t1) = (", y[1], ",", e[1], ")"),
       cex = 1.1, pt.cex = 1, y.intersp = 1.1)

## ----mse.lambdapeeb, echo=TRUE-------------------------------------------
# ------------------------------------------------------------- #
# mse.lambdapeeb: returns MSE for the posterior mean
# 
# parameters:
# a, b - hyperparameters of the prior Gamma distribution
# t - correction scalar for the Poisson likelihood distribution
# lambda - observed value of the parameter
# ------------------------------------------------------------- #
mse.lambdapeeb <- function(a, b, t, lambda) {
  (b/(b * t + 1))^2 * ((lambda/b)^2 + lambda * (t - 2 * a/b) + a^2)
}

## ----MSE.plots, echo=FALSE, fig.cap="Plots of MSE for the posterior mean and the MLE for $\\lambda$.", cache=TRUE, fig.height=6----
# plotting MSE for PM and MLE
lambda.grid <- seq(0, 2, length.out = 1000)

plot(lambda.grid, mse.lambdapeeb(5, 1/4, 10, lambda.grid),
     type = "l", lwd = 2, col = "darkorange1", 
     xlab = expression(lambda), ylab = "MSE", lty = 5, ylim = c(0, 0.2))

lines(lambda.grid, lambda.grid/10, type = "l", lwd = 2, col = "dodgerblue3")

legend("topleft", inset = .05, legend = c("PEEB", "ML"),
       col = c("darkorange1", "dodgerblue3"), lwd = 2, bty = "n",
       cex = 1.1, pt.cex = 1, y.intersp = 1.8, lty = c(5, 1))

## ----EXERCISE4, echo=FALSE-----------------------------------------------
# ------------------------------------------------------------- #
# EXERCISE 4
# Carrying out a Bayesian hierarchical models analysis from the 
# paper 'Methods for estimating and interpreting provider-specific
# standardized mortality ratios', Liu et al., 2003, using 
# dialysis facility report data FY2017.
# ------------------------------------------------------------- #

## ----data.processing, eval=FALSE-----------------------------------------
## # reading the file, converting extracted values to numeric type
## dfrdata <- read.csv("DFR_Data_FY2017.csv", header = TRUE)
## dataNY <- subset(dfrdata, state == "NY")
## 
## dia <- subset(dataNY, select=c(dyy4_f, exdy4_f, deay4_f, smry4_f))
## colnames(dia) <- c("yrs_at_risk", "mu", "y", "rho_ml")
## 
## dia$yrs_at_risk <- as.numeric(paste(dia$yrs_at_risk))
## dia$mu <- as.numeric(paste(dia$mu))
## dia$y <- as.numeric(paste(dia$y))
## dia$rho_ml <- as.numeric(paste(dia$rho_ml))
## 
## # recalculating missing rhos
## missingRho <- which(is.na(dia$rho_ml))
## for (i in 1:length(missingRho)) {
##   dia$rho_ml[missingRho[i]] <- dia$y[missingRho[i]]/dia$mu[missingRho[i]]
## }
## # eliminating NA values
## dia <- dia[complete.cases(dia), ]
## write.csv(dia, file = "dia.csv")

## ----readingDIA----------------------------------------------------------
dia <- read.csv("dia.csv", row.names = 1, header = TRUE)

## ----table1, echo=FALSE--------------------------------------------------
# ------------------------------------------------------------- #
# percentageProviderStrata: returns a list containing the 
# percentage of providers in chosen strata with MLE-estimated 
# percentiles in lower 25%, middle 50% and upper 25%, and the
# number of providers in the strata. 
# 
# Parameters:
# lsize - lower limit for the strata
# usize - upper limit for the strata
# qDia - vector of global quantiles that we wish to compare to
# ------------------------------------------------------------- #
percentageProviderStrata <- function(lsize, usize, qDia) {
  data <- dia[dia$yrs_at_risk >= lsize & dia$yrs_at_risk < usize, ]
  numProviders <- nrow(data)
  rhoMlProviders <- data$rho_ml

  lower25 <- sum(rhoMlProviders < qDia[1]) * 100/numProviders
  middle50 <- sum(rhoMlProviders >= qDia[1] & 
                    rhoMlProviders < qDia[2]) * 100/numProviders
  upper25 <- sum(rhoMlProviders >= qDia[2]) * 100/numProviders
  return(round(c(lower25, middle50, upper25, numProviders)))
}

# ------------------------------------------------------------- #
# percentageProviders: returns a data.frame corresponding to Table 1
# 
# Parameters: 
# strataPoints - vector of points determining strata sizes
# ------------------------------------------------------------- #
percentageProviders <- function(strataPoints) {
  qDia <- quantile(dia$rho_ml, c(0.25, 0.75)) # global quantiles
  n <- length(strataPoints)
  tableDF <- data.frame(matrix(nrow = n + 1, ncol = 5), 
                        stringsAsFactors = FALSE)
  colnames(tableDF) <- c("Provider size", "Lower 25%", "Middle 50%", 
                         "Upper 25%", "No. of providers")
  strataPoints <- c(0, strataPoints, Inf)
  for (i in 1:(n + 1)) {
    tableDF[i, 2:5] <- percentageProviderStrata(strataPoints[i], 
                                                strataPoints[i + 1], qDia)
  }
  
  # For displaying the results
  strata <- function(lp, up) {
    if (lp == 0) {
      paste("<", toString(up))
    } else if (up == Inf) {
      paste(">=", toString(lp))
    } else {
      paste(toString(lp), "-", toString(up))
    }
  }
  
  tableDF[, 1] <- mapply(strata, strataPoints[1:(n + 1)], 
                         strataPoints[2:(n + 2)])
  return(tableDF)
}

## ----MLE.SMR.table, cache=TRUE, echo=FALSE-------------------------------
# Displaying percentage of providers with SMR in 3 different categories
kable(percentageProviders(c(25, 50, 100, 150)), booktabs = T,
      caption = "Percentage of providers in small, moderate and 
      large MLE SMR category compared to the whole \\textit{dia} dataset.", 
      align = 'c') %>% 
  kable_styling(position = "center") %>% 
  row_spec(0, bold = T)

## ----createSMRmodel------------------------------------------------------
# ------------------------------------------------------------- #
# createSMRmodel: returns mcmc.list of samples of SMRs for the 
# Poisson-Normal-Gamma model from the paper.
# 
# Parameters:
# sigma2 - squared variance of zeta hyperparameter
# a - the Gamma parameter for lambda hyperparameter
# nChains - desired number of parallel chains
# nSamples - desired number of samples 
# burnin - desired size of the burnin
# ------------------------------------------------------------- #
createSMRmodel <- function(sigma2, a, nChains = 1, nSamples, burnin) {
  SMRmodel.string <- "model{
    for (i in 1:n) {
      par[i] <- mu[i] * rho[i]
      y[i] ~ dpois(par[i])
      rho[i] <- exp(theta[i])
      theta[i] ~ dnorm(zeta, lambda)
    }
    
    # Priors
    zeta ~ dnorm(0, 1/sigma2)
    lambda ~ dgamma(a, a)
  }"
  
    
  SMRmodel <- jags.model(textConnection(SMRmodel.string), 
                         data = list('n' = nrow(dia), 'y' = dia$y, 
                                     'mu' = dia$mu, 'sigma2' = sigma2,
                                     'a' = a), 
                         n.chains = nChains, n.adapt = 100, quiet = TRUE)

  samples <- window(coda.samples(SMRmodel,
                                 variable.names = c('rho'), 
                                 n.iter = nSamples + burnin), 
                    start = burnin + 101)
  return(samples)
}

## ----runDiagnostic-------------------------------------------------------
# ------------------------------------------------------------- #
# runDiagnostic: returns gelman.diag object evaluated on the 
# mcmc.list from the above JAGS model.
# 
# Parameters:
# sigma2 - squared variance of zeta hyperparameter
# a - the Gamma parameter for lambda hyperparameter
# nChains - desired number of parallel chains
# nSamples - desired number of samples 
# burnin - desired size of the burnin
# ------------------------------------------------------------- #
runDiagnostic <- function(sigma2, a, nChains = 1, nSamples = 2500, burnin) {
  mysamples <- createSMRmodel(sigma2, a, nChains, nSamples, burnin)
  # Attempt of using trace and gelman plots - cumbersome hence not used.
  # gr <- gelman.plot(mysamples)
  # shrink <- gr$shrink[1, , ] 
  # psrf <- gelman.diag(window(mysamples, end=gp$last.iter[1]), multivariate=FALSE) 
  # medMean <- mean(shrink[, 1] - psrf$psrf[, 1])
  # uppCIMean <- mean(shrink[, 1] - psrf$psrf[, 1])
  
  gd <- gelman.diag(mysamples)
  return(gd)
}

## ----convergenceDiagnostic, cache=TRUE-----------------------------------
# Testing the convergence of chains
testConv1 <- runDiagnostic(10^6, 10^(-4), nChains = 2, 
                           nSamples = 2500, burnin = 1000)
testConv2 <- runDiagnostic(10^6, 10^(-4), nChains = 2, 
                           nSamples = 2500, burnin = 5000)
testConv3 <- runDiagnostic(10^6, 10^(-4), nChains = 3, 
                           nSamples = 2500, burnin = 1000)
testConv4 <- runDiagnostic(10^6, 10^(-4), nChains = 3, 
                           nSamples = 2500, burnin=5000)

## ----rhoSamples, cache=TRUE----------------------------------------------
# Samples with accepted convergence
rhoSamples <- createSMRmodel(10^6, 10^(-4), 3, 2500, 5000)

## ----rhoSamples.Table1, cache=TRUE, fig.cap="Posterior means of SMR and their $95\\%$ credibility intervals.", fig.height=6----
# plotting PM and their 95% CI
xaxis <- seq(1, nrow(dia))
posteriorMeans <- summary(rhoSamples)$statistics[, "Mean"]
posteriorCI <- summary(rhoSamples)$quantiles[, c(1, 5)]
posteriorMeansCI <- data.frame(posteriorMeans, posteriorCI)
posteriorMeansCI <- posteriorMeansCI[order(posteriorMeans), ]

plot(xaxis,posteriorMeansCI$posteriorMeans, ylim = c(0, 3), 
     pch = 18, ylab = "PM SMR", xlab = "Provider") # means
arrows(xaxis, posteriorMeansCI$X2.5., xaxis, posteriorMeansCI$X97.5., 
       length = 0.05, angle = 90, code = 3) # CIs
abline(h = 1)

## ----MLESMRs.Figure1, cache=TRUE, fig.cap="MLE SMRs and their $95\\%$ confidence intervals.", fig.height=6----
# plotting MLE and their 95% CI
mlLowCI <- c()
mlUppCI <- c()
for (i in 1:nrow(dia)) {
  PoissonTest <- poisson.test(dia$y[i], T = dia$mu[i], r = dia$rho_ml[i], 
                              alternative = "two.sided")$conf.int
  mlLowCI[i] <- PoissonTest[1]
  mlUppCI[i] <- PoissonTest[2]
}

mlMeansCI <- data.frame(mlMeans = dia$rho_ml, mlLowCI, mlUppCI)
mlMeansCI <- mlMeansCI[order(mlMeansCI$mlMeans), ]

plot(xaxis,mlMeansCI$mlMeans, pch = 18, ylab = "MLE SMR", 
     xlab = "Provider", ylim = c(0, 12)) # means
arrows(xaxis, mlMeansCI$mlLowCI, xaxis, mlMeansCI$mlUppCI, 
       length = 0.05, angle = 90, code = 3) # CIs
abline(h = 1)

## ----smallMLESMRs, echo=FALSE--------------------------------------------
# yrs_at_risk for providers with small SMR
smallMLESMRs <- which(dia$rho_ml < 0.01)
smallMLESMRs.yrs <- mean(dia$yrs_at_risk[smallMLESMRs])

## ----provnames1, eval=FALSE----------------------------------------------
## # Extracting provider names
## dfrdata <- read.csv("DFR_Data_FY2017.csv", header = TRUE)
## dfrdataProviders <- subset(dfrdata, state == "NY", select = c(provname) )
## dfrdataProviders$provname <- as.character(dfrdataProviders$provname)
## colnames(dfrdataProviders) <- c("providers")
## write.csv(dfrdataProviders, file = "provnames.csv")

## ----provnames2, cache=TRUE----------------------------------------------
dfrdataProviders <- read.csv("provnames.csv", row.names = 1, 
                             header = TRUE)

# Extracting data for the Westchester provider
westIndex <- grep("WESTCHESTER CENTER FOR RENAL CARE", 
                  dfrdataProviders$providers)
westProvider <- dia[westIndex, ]

## ----westMeans1----------------------------------------------------------
# computing ML and PM of the Westchester provider
westML <- westProvider$rho_ml # MLE
westPM <- posteriorMeans[westIndex] # PM

## ----westEBmean----------------------------------------------------------
# Computing EBPM for the Westchester provider
priorA <- optim(c(1,1), mloglik, y = dia$y, t = dia$mu, 
                control = list(fnscale = -1))$par
westEBmean <- (priorA[1] + westProvider$y)*
  (priorA[2]/(westProvider$mu * priorA[2] + 1))

## ----westProvider--------------------------------------------------------
print(westProvider)

## ----westDiffMeans, cache=TRUE, echo=FALSE-------------------------------
# Displaying MLE, PM and PMEB for the Westchester provider
westM <- data.frame(matrix(nrow = 1, ncol = 3))
colnames(westM) <- c("MLE", "PM", "PMEB")
westM[1, ] <- c(westML, westPM, westEBmean)
kable(westM, booktabs = T, caption = "MLE, PM and PMEB for the Westchester Center for Renal Care facility.", 
      align = 'c') %>% 
  kable_styling(position = "center") %>% 
  row_spec(0, bold = T)

## ----posteriorHypothesis-------------------------------------------------
# Posterior probability for the hypothesis
posteriorHyp <- mean(unlist(rhoSamples[ , westIndex]) > 1)

## ----priorEB-------------------------------------------------------------
# EB prior and posterior probabilities of the hypothesis
priorEB <- integrate(function(x) 
  dgamma(x, shape = priorA[1], scale = priorA[2]), 1, Inf)

postEB <- integrate(function(x) 
  dgamma(x, shape = (priorA[1] + westProvider$y),
         scale = (priorA[2]/(westProvider$mu * priorA[2] + 1))), 1, Inf)

## ----priorPostWest, cache=TRUE, echo=FALSE-------------------------------
# Displaying prior and posterior probabilities of the hypothesis
PrPoWest <- data.frame(matrix(nrow = 2, ncol = 2))
colnames(PrPoWest) <- c("Prior rho_211 > 1", "Posterior rho_211 > 1")
rownames(PrPoWest) <- c("Liu et al.", "Empirical Bayes")
PrPoWest[1, ] <- c(0.5, posteriorHyp)
PrPoWest[2, ] <- c(round(priorEB$value, 2), postEB$value)
kable(PrPoWest, booktabs = T, caption = "Prior and posterior probabilities for the hypothesis $\\rho_{211} > 1$ under the two models.", 
      align = 'c') %>% 
  kable_styling(position = "center") %>% 
  row_spec(0, bold = T)

## ----extendedJags--------------------------------------------------------
# ------------------------------------------------------------- #
# createRankmodel: returns mcmc.list of samples of ranks for the 
# Poisson-Normal-Gamma model from the paper.
# 
# Parameters:
# sigma2 - squared variance of zeta hyperparameter
# a - the Gamma parameter for lambda hyperparameter
# nChains - desired number of parallel chains
# nSamples - desired number of samples 
# burnin - desired size of the burnin
# ------------------------------------------------------------- #
createRankmodel <- function(sigma2, a, nChains = 1, nSamples, burnin) {
  Rankmodel.string <- "model{
    for (i in 1:n) {
      par[i] <- mu[i] * rho[i]
      y[i] ~ dpois(par[i])
      rho[i] <- exp(theta[i])
      theta[i] ~ dnorm(zeta, lambda)
    }
  
    r <- rank(rho[]) 
    pr <- r*100/n
  
    # Priors
    zeta ~ dnorm(0, 1/sigma2)
    lambda ~ dgamma(a, a)
  }"
  
  Rankmodel <- jags.model(textConnection(Rankmodel.string),
                         data = list('n' = nrow(dia), 'y' = dia$y, 
                                     'mu' = dia$mu, 'sigma2'=sigma2,
                                     'a' = a), 
                         n.chains = nChains, n.adapt = 100, quiet = TRUE)
  
  samples <- window(coda.samples(Rankmodel,
                                 variable.names = c('r'), 
                                 n.iter = nSamples + burnin), 
                    start = burnin + 101)
  return(samples)
}

## ----top3, cache=TRUE----------------------------------------------------
# Computing top3 facilities from the updated JAGS program
rankSamples <- createRankmodel(10^6, 10^(-4), 3, 2500, 5000)

top3Probs <- c()
for (i in 1:255) {
  top3Probs[i] <- mean(unlist(rankSamples[, i]) <= 3)
}
# indices of top3 facilities in dia.csv
top3ProbIndex <- order(top3Probs, decreasing = TRUE)[1:3]

# determining the names of top3 facilities
diaRows <- as.numeric(rownames(dia)[top3ProbIndex])
top3Names <- dfrdataProviders[paste(diaRows[1:3]), ]

# top3 posterior probabilities and MLEs
top3Prob <- top3Probs[top3ProbIndex]
top3MLE <- dia$rho_ml[top3ProbIndex]

# generating data.frame
top3 <- data.frame(top3Prob, top3MLE, stringsAsFactors = FALSE)
colnames(top3) <- c("Probability of being top3", "MLE")
topp3 <- cbind("provname" = top3Names, top3)

## ----top3.table, cache=TRUE, echo=FALSE----------------------------------
# Displaying the obtained top3 facilities
kable(topp3, booktabs = T,
      caption = "Top 3 facilities.", 
      align = 'c') %>% 
  kable_styling(position = "center", font_size = 8) %>% 
  row_spec(0, bold = T)

## ----top3Provinfo, size="scriptsize"-------------------------------------
print(cbind("provname" = top3Names, dia[top3ProbIndex, ]))

