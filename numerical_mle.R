# Numerical estimation of the MLE of a multivariate model

# simulate exponential regression
n = 100
p = 3

# First we generate the 'true' data, i.e. the covariates, the beta
# (true parameter), and the mean
X = matrix(rnorm(n*p), nrow = n, ncol = p)
X = cbind(1, X)

betaVec = c(1, 1, 0.5, -1)

# each Y_i has its own mean, mu_i, which we get by combining each x_i with beta
# and taking the exp. I.e. Y_i = exp(beta_1 + x_i2*beta_2 + x+i3*beta_3 + x_i4*beta_4)
mu = exp(X %*% betaVec)


# Now we generate the observations y from an exponential distribution with mean
# of the true data. However, `rexp` expects a 'rate' argument with is 1/mean
y = rexp(n, 1/mu)

# plot the response to the first covariate (x_2)
# since the mean is an exponential, we need to plot the log of y
# to see a linear dependency in the plot
# plotting the second feature, corresponding to beta=1 will give a positive
# dependency, while plotting the fourth feature (beta=-1) will give a negative
# dependency
plot(X[,4], log(y), pch=19)

# now we want to estimate the parameters - beta - from the observations
# and the covariates. This is done with the 'optim' 
# define the log likelihood function, to be used by optim
logLikeRegExp = function(beta, y, X) {
  # since rate is 1/mean
  mu = exp(X %*% beta)
  # use the log density function, and return the sum
  return (sum(dexp(y, rate = 1/mu, log = TRUE)))
}

# now optimize to find the MLE
# find initial values by linear regression
lmFit = lm(log(y) ~ 0 + X)
initVal = lmFit$coefficients

optRes = optim(initVal, fn = logLikeRegExp, gr = NULL, method = "BFGS", y = y, X = X, 
               control = list(fnscale = -1), hessian = TRUE)
betaHat = optRes$par

# Let's find the standard error
# The negative hessian at betaHat is equivalent to the observed Fisher information
# and the inverse of the negative hessian is an estimator of the asymptotic
# covariance matrix (https://stats.stackexchange.com/a/68095/494223)
covBetaHat = -solve(optRes$hessian)
stdErr = sqrt(diag(covBetaHat))

# to calculate the 95% confidence interval for betaHat we need to multiply
# the standard error with 1.96. This is the value z for which P(betaHat > z) = (1-0.95)/2
# (see http://www.stat.yale.edu/Courses/1997-98/101/confint.htm)
# the number can be retrieved by qnorm(1 - (1-0.95)/2)
# So the interval is (betaHat-1.96*stdErr, betaHat+1.96*stdErr)
confInterval = r(betaHat-1.96*stdErr, betaHat+1.96*stdErr)
confInterval


