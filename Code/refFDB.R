library(refund)
library(FDboost)
library(R.matlab)

setwd('/Users/mjm556/Documents/MATLAB')
Yt  <- readMat('Y.mat')
Y   <- Yt$Y

Xt  <- readMat('simX.mat')
X   <- Xt$simX

modelR  <- pffr(Y ~ ff(X, limits = "s<t",
                       splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)),
                                         k = c(10, 10))) - 1)

coefs   <- coef(modelR)$smterms$ff$value
se      <- coef(modelR)$smterms$ff$se
bhmat   <- matrix(coefs, nrow = 40, ncol = 40)
semat   <- matrix(se, nrow = 40, ncol = 40)
lower   <- bhmat - qnorm(0.975)*semat
upper   <- bhmat + qnorm(0.975)*semat

writeMat('bhmatR.mat', bhmat = bhmat)
writeMat('lowerR.mat', lower = lower)
writeMat('upperR.mat', upper = upper)

dat <- list(Y = Y, X = X, t = 1:ncol(Y), s = 1:ncol(X))

modelF <- FDboost(Y ~ 1 + bhist(x = X, s = s, time = t), timeformula = ~ bbs(t, knots = 10), data = dat)

bhatF 	<- coef(modelF)$smterms$bhist$value

writeMat('bhmatF.mat', bhmatF = bhatF)


