# Modest smoke simulation, not a reproduction of the paper's 10,000-run tables.
library(spikedGLHT)
set.seed(20261018)
p <- 50; N <- 100
X <- rbind(c(rep(1,40),rep(0,60)),c(rep(0,40),rep(1,60)))
C <- matrix(c(1,-1),2,1)
repetitions <- 100L
rejections <- replicate(repetitions, {
    Y <- matrix(rnorm(p*N),p,N)
    Y[1:3,] <- sweep(Y[1:3,],1,sqrt(c(20,12,8)),`*`)
    spiked_glht(Y,X,C,K=3,calibration='bootstrap',B=499)$reject
})
print(rowMeans(rejections))
cat('Monte Carlo standard error near 5%:',sqrt(0.05*0.95/repetitions),'\n')
