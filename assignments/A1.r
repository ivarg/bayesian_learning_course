# Problem 1
# Let y1..yn | theta ~ Bern(theta), and assume that you have obtained a sample 
# with s=14 successes in n=20 trials. Assume a Beta(alfa_0, beta_0) prior for 
# theta and let alfa_0 = beta_0 = 2.

n = 20
s = 14
f = n-s
alfa_0 = beta_0 = 2

# a)
# Draw random numbers from the posterior theta | y ~ Beta(alfa_0 + s, beta_0 + f),
# where y = (y_1, ... , y_n) and verify graphically that the Monte Carlo (MC) 
# estimates of the posterior mean and standard deviation converges to the true 
# values as the number of random draws grows large.

# Posterior parameters
alfa_p = alfa_0 + s
beta_p = beta_0 + f

# Posterior mean and variance
postMean = alfa_p / (alfa_p + beta_p)
postVar = sqrt((alfa_p * beta_p)/((alfa_p + beta_p)^2 * (alfa_p + beta_p + 1)))

numSamples = 10000
post = rbeta(numSamples, alfa_p, beta_p)

sampleMeans = seq(1, numSamples)
sampleVars = seq(1, numSamples)
for (i in 1:numSamples) {
  sample = rbeta(i, alfa_p, beta_p)
  sampleMeans[i] = sum(sample)/i
  sampleVars[i] = sqrt(var(sample))
}

par(mfcol=c(1,2))

plot(sampleMeans, type="l", yaxt="n", xlab="Sample size", ylab="", col="steelblue", main="Sample mean")
abline(a=postMean, b=0, lty="dashed")
axis(2, at=c(postMean), labels=c("Posterior mean"))
#axis(2, at=c(min(sampleMeans), postMean, max(sampleMeans)),
#     labels=c(round(min(sampleMeans), digits=2), "Posterior mean", round(max(sampleMeans), digits=2)))

plot(sampleVars, type="l", xlab="Sample size", ylab="", col="indianred", yaxt="n", main="Sample sd")
abline(a=postVar, b=0, lty="dashed")
axis(2, at=c(postVar), labels=c("Posterior sd"))

par(mfcol=c(1,1))

