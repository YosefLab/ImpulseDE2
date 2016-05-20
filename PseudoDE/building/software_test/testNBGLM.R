library(MASS)
data <- c(rnbinom(n=10,mu=20,size=1),rnbinom(n=10,mu=50,size=1),rnbinom(n=10,mu=30,size=1))
data
groups <- c(rep("cluster1",10),rep("cluster2",10),rep("cluster3",10))
groups
thetahat <- 1.5
coefs_mu <- c(10,10,10)
fit <- glm.nb(data ~ groups, init.theta = thetahat, start=coefs_mu)
fitted=fit$fitted.values
theta=fit$theta
coefs=fit$coefficients
exp(coefs[1])
exp(coefs[1]+coefs[2])
exp(coefs[1]+coefs[3])

