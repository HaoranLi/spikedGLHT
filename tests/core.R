library(spikedGLHT)
near <- function(a, b, tol = 1e-8) stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error=identity), 'error'))
# Deterministic synthetic data; no dependency on research-script output.
set.seed(20261018)
p <- 40; N <- 80
X <- rbind(1, rep(c(-1,1),each=N/2), rnorm(N))
C <- diag(3)[,2:3,drop=FALSE]
Y <- matrix(rnorm(p*N),p,N)
Y[1:2,] <- Y[1:2,] * c(6,4)

# Independent direct implementation of paper Algorithm 1 and Corollary 3.1.
p <- nrow(Y); N <- ncol(Y); m <- nrow(X); n <- N-m; q <- ncol(C)
S <- Y %*% (diag(N)-t(X)%*%solve(X%*%t(X),X)) %*% t(Y)/n
es <- eigen(S,symmetric=TRUE)
bulk <- mean(es$values[3:p]); gamma <- p/n
root <- function(v,s) { z<-v/s; b<-gamma-1-z; (-b+sqrt(pmax(b*b-4*z,0)))*s/2 }
l1 <- root(es$values[1:2],bulk)
sigma <- bulk*(1+gamma*(2+sum(1/(l1/bulk-1)))/(p-2))
ell <- root(es$values[1:2],sigma)
ell[es$values[1:2]/sigma<(1+sqrt(gamma))^2] <- sigma*(1+sqrt(gamma))
b <- spiked_glht(Y,X,C,K=2)
near(b$sigma_sq,sigma); near(b$spikes,ell)
z2 <- pmax(0,(1-gamma/(ell/sigma-1)^2)/(1+gamma/(ell/sigma-1)))
lambda <- b$weights
T1 <- (sum(lambda[-1]*(z2*ell+(1-z2)*sigma))+lambda[1]*(sum(ell)+(p-2)*sigma))/p
T2 <- (sum(lambda[-1]^2*(z2*ell+(1-z2)*sigma)^2 +
       2*lambda[1]*lambda[-1]*(z2*ell^2+(1-z2)*sigma^2))+
       lambda[1]^2*(sum(ell^2)+(p-2)*sigma^2))/p
near(b$theta1,T1); near(b$theta2,T2)
stopifnot(abs(sum(b$weights^2)-1)<1e-10,
          abs(sum(b$least.favorable.prior)-1)<1e-10,
          min(b$least.favorable.prior)>-1e-8,b$least.favorable.prior[1]<=0.5+1e-8)

# Rotation invariance, and nonzero D equivalent to removing a compatible mean.
set.seed(41)
U <- qr.Q(qr(matrix(rnorm(p*p),p,p)))
c <- spiked_glht(U%*%Y,X,C,K=2)
near(b$statistic,c$statistic)
Beta <- matrix(rnorm(p*m),p,m)
d <- spiked_glht(Y+Beta%*%X,X,C,D=Beta%*%C,K=2)
near(b$statistic,d$statistic)

# Single-spike, isotropic, p>N, selection, boundary and bootstrap reproducibility.
for (kk in c(0,1)) {
  z <- spiked_glht(Y,X,C,K=kk)
  stopifnot(all(is.finite(z$statistic)),length(z$weights)==kk+1)
}
set.seed(52)
wide <- matrix(rnorm(120*N),120,N); wide[1,] <- wide[1,]*8
w <- spiked_glht(wide,X,C,K=1)
stopifnot(all(is.finite(w$statistic)))
for (sel in c('PLY','KN')) {
 if(sel=='KN' && !requireNamespace('RMTstat',quietly=TRUE)) next
 s <- spiked_glht(Y,X,C,method=sel,K_max=3,delta=0.1)
 stopifnot(s$K>=0,s$K<=3)
}
zero <- spiked_glht(Y,X,C,K_max=0)
stopifnot(zero$K==0)
for (kk in c(0,1,2)) {
 set.seed(61)
 h <- spiked_glht(Y,X,C,K=kk,calibration='bootstrap',B=99)
 set.seed(61)
 j <- spiked_glht(Y,X,C,K=kk,calibration='bootstrap',B=99)
 near(h$bootstrap.statistics,j$bootstrap.statistics)
 stopifnot(all(h$p.value>=0.01),all(h$p.value<=1),
           identical(dim(h$bootstrap.statistics),c(99L,3L)))
 # Sample bootstrap means agree with independently computed centering.
 set.seed(62)
 h <- spiked_glht(Y,X,C,K=kk,calibration='bootstrap',B=3000)
 stopifnot(abs(mean(h$bootstrap.statistics[,'LH'])-q*h$theta1) <
           5*sd(h$bootstrap.statistics[,'LH'])/sqrt(3000))
 lh <- h$bootstrap.statistics[,'LH']
 variance_se <- sqrt((mean((lh-mean(lh))^4)-var(lh)^2)/length(lh))
 stopifnot(abs(var(lh)-2*q*h$theta2/p)<5*variance_se)
}
# Exactly zero overlaps must not make the prior optimization singular.
g <- spikedGLHT:::.weight_geometry(c(1+sqrt(p/n),8),1,n,p)
z <- spikedGLHT:::.minimax(g,0.2)
stopifnot(all(is.finite(z$weights)),abs(sum(z$prior)-1)<1e-8,z$prior[1]<=0.2+1e-8)

fails(spiked_glht(Y,X,C,K=-1))
fails(spiked_glht(Y,X,C,K=1.5))
fails(spiked_glht(Y,X,C,rho=0))
fails(spiked_glht(Y,X,C,alpha=1))
fails(spiked_glht(Y,X,C,D=matrix(0,q,p)))
fails(spiked_glht(Y,rbind(X,X[1,]),rbind(C,C[1,])))
fails(spiked_glht(Y,X,cbind(C,C[,1])))
fails(spiked_glht(Y,X,C,K=1,calibration='bootstrap',B=0))
fails(spiked_glht(Y,X,C,K_max=n))
fails(spiked_glht(replace(Y,1,NA),X,C))
cat('All numerical and validation checks passed.\n')

# With a single contrast the raw transforms are monotone, hence bootstrap ranks agree.
set.seed(77)
one <- spiked_glht(Y,X,C[,1,drop=FALSE],K=1,calibration='bootstrap',B=99)
stopifnot(length(unique(one$p.value))==1L,length(unique(one$reject))==1L)
