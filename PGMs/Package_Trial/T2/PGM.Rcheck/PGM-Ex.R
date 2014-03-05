pkgname <- "PGM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PGM')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Bsublin")
### * Bsublin

flush(stderr()); flush(stdout())

### Name: Bsublin
### Title: Sublinear function
### Aliases: Bsublin
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, R0 = 0) 
{
    Bx = X
    Bx[X > R] = (R + R0)/2
    ind = X > R0 & X <= R
    Bx[ind] = (-X[ind]^2 + 2 * R * X[ind] - R0^2)/(2 * (R - R0))
    return(Bx)
  }



cleanEx()
nameEx("Copula.Norm.Pois")
### * Copula.Norm.Pois

flush(stderr()); flush(stdout())

### Name: Copula.Norm.Pois
### Title: Copula transform a matrix from normal to Poisson
### Aliases: Copula.Norm.Pois
### Keywords: ~kwd1 ~kwd2

### ** Examples

		X <- matrix(rnorm(20), nrow=5, ncol=4)
		transX <- Copula.Norm.Pois(X, lambda=1)

## The function is currently defined as
function (X, lambda) 
{
    n = nrow(X)
    p = ncol(X)
    val = 0
    dcuts = NULL
    cnt = 0
    while (val < max(0.9999, 1 - 2/(n * p))) {
        val = ppois(cnt, lambda)
        cnt = cnt + 1
        dcuts = c(dcuts, val)
    }
    Y = matrix(0, n, p)
    oldval = min(X)
    for (i in 1:length(dcuts)) {
        val = quantile(X, dcuts[i])
        Y[which(X < val & X >= oldval)] = i - 1
        oldval = val
    }
    Y[X == max(X)] = max(Y) + 1
    return(Y)
  }



cleanEx()
nameEx("LPGM.select")
### * LPGM.select

flush(stderr()); flush(stdout())

### Name: LPGM.select
### Title: Log-Linear Graphical Model based on Pair-wise Poisson Markov
###   Network
### Aliases: LPGM.select
### Keywords: ~kwd1 ~kwd2

### ** Examples

library(PGM)
library(huge)
n = 200
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace=TRUE), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

#-# Run LPGM
lpgm.path.all = LPGM.select(t(simData), nlams=20, N=10, beta=0.05, parallel=FALSE)
str(lpgm.path.all)

## The function is currently defined as
function (X, method = "LPGM", N = 100, beta = 0.05, lmin = 0.01, 
    nlams = 20, lambda.path = NULL, parallel = T, nCpus = 4) 
{
    if (is.null(lambda.path)) {
        lmax = myglmnet.max(X)
        lambda.path = exp(seq(log(lmax), log(lmin), l = nlams))
    }
    if (parallel == T) {
        b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
        ghat = list()
        ghat.path = list()
        ghat.path$path = vector("list", length(lambda.path))
        v = c()
        for (i in 1:N) {
            cat(paste(method, ": Conducting sampling ... in progress: ", 
                floor(100 * (i/N)), "%", collapse = ""), "\r")
            flush.console()
            glmpois.good <- 1
            while (glmpois.good) {
                good <- 1
                while (good) {
                  index = sample(1:ncol(X), b, replace = F)
                  if (sum(apply(X[, index], 1, function(x) length(unique(x)) == 
                    1)) == 0) {
                    good <- 0
                  }
                }
                tryCatch({
                  ghat.path$raw = glmpois(X[, index], lambda = lambda.path, 
                    parallel = T, nCpus = nCpus)
                  glmpois.good <- 0
                }, error = function(e) {
                  cat("glmnet returns empty model. Try again.")
                })
            }
            for (j in 1:length(lambda.path)) {
                tmp = ghat.path$raw[, , j]
                tmp[abs(tmp) < 1e-06] = 0
                tmp[abs(tmp) > 1e-06] = 1
                diag(tmp) = 0
                if (is.null(ghat.path$path[[j]])) {
                  ghat.path$path[[j]] = tmp
                }
                else {
                  ghat.path$path[[j]] = ghat.path$path[[j]] + 
                    tmp
                }
            }
        }
        for (i in 1:length(lambda.path)) {
            D = ghat.path$path[[i]]
            D = D/N
            D = 2 * D * (1 - D)
            v = c(v, mean(D[upper.tri(D)]))
        }
        v = cummax(v)
        ghat$v = v
        ghat$lambda.path = lambda.path
        ghat$opt.lambda = lambda.path[which(v == max(v[v < beta]))]
        ghat$network = glmpois(X, lambda = lambda.path, parallel = T, 
            nCpus = nCpus)
        ghat$network = lapply(1:nlams, function(r) {
            return(ghat$network[, , r])
        })
        ghat$opt.index = which(v == max(v[v < beta]))
        cat(paste("\n", method, " Completed.", "\n", sep = ""))
        return(ghat)
    }
    if (parallel == F) {
        b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
        ghat = list()
        v = c()
        for (j in 1:length(lambda.path)) {
            cat(paste(method, ": Conducting sampling ... in progress: ", 
                floor(100 * (i/N)), "%", collapse = ""), "\r")
            flush.console()
            D = matrix(0, nrow = nrow(X), ncol = nrow(X))
            for (i in 1:N) {
                glmpois.good <- 1
                while (glmpois.good) {
                  good <- 1
                  while (good) {
                    index = sample(1:ncol(X), b, replace = F)
                    if (sum(apply(X[, index], 1, function(x) length(unique(x)) == 
                      1)) == 0) {
                      good <- 0
                    }
                  }
                  tryCatch({
                    tmp = glmpois(X[, index], lambda = lambda.path[j], 
                      parallel = F)
                    glmpois.good <- 0
                  }, error = function(e) {
                    cat("glmnet returns empty model. Try again.\n")
                  })
                }
                tmp[abs(tmp) < 1e-06] = 0
                tmp[abs(tmp) > 1e-06] = 1
                D = D + tmp
            }
            D = D/N
            D = 2 * D * (1 - D)
            v = c(v, mean(D[upper.tri(D)]))
        }
        v = cummax(v)
        ghat$v = v
        ghat$lambda.path = lambda.path
        ghat$opt.lambda = lambda.path[which(v == max(v[v < beta]))]
        ghat$network = glmpois(X, lambda = lambda.path, parallel = parallel, 
            nCpus = nCpus)
        ghat$network = lapply(1:nlams, function(r) {
            return(ghat$network[, , r])
        })
        ghat$opt.index = which(v == max(v[v < beta]))
        cat(paste("\n", method, " Completed.", "\n", sep = ""))
        return(ghat)
    }
  }



cleanEx()
nameEx("SPGM.select")
### * SPGM.select

flush(stderr()); flush(stdout())

### Name: SPGM.select
### Title: Log-Linear Graphical Model based on Pair-wise Sub-linear
###   truncated Poisson Markov Network
### Aliases: SPGM.select
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, R0 = 0, N = 100, beta = 0.05, lmin = 0.01, nlams = 20, 
    lambda.path = NULL, parallel = T, nCpus = 4) 
{
    if (R < 0) {
        cat("ERROR: Truncating threshold R should be positive. \n")
        ghat = NULL
        return(ghat)
    }
    Xorig <- X
    X <- round(Bsublin(X, R, R0))
    return(LPGM.select(X, method = "SPGM", N = N, beta = beta, 
        lmin = lmin, nlams = nlams, lambda.path = lambda.path, 
        parallel = parallel, nCpus = nCpus))
  }



cleanEx()
nameEx("TPGM.select")
### * TPGM.select

flush(stderr()); flush(stdout())

### Name: TPGM.select
### Title: Log-Linear Graphical Model based on Pair-wise truncated Poisson
###   Markov Network
### Aliases: TPGM.select
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, N = 100, beta = 0.05, lmin = 0.01, nlams = 20, 
    lambda.path = NULL, parallel = T, nCpus = 4) 
{
    if (R < 0) {
        cat("ERROR: Truncating threshold R should be positive. \n")
        ghat = NULL
        return(ghat)
    }
    Xorig <- X
    X[X > R] <- R
    return(LPGM.select(X, method = "TPGM", N = N, beta = beta, 
        lmin = lmin, nlams = nlams, lambda.path = lambda.path, 
        parallel = parallel, nCpus = nCpus))
  }



cleanEx()
nameEx("WPGM.neighborhood")
### * WPGM.neighborhood

flush(stderr()); flush(stdout())

### Name: WPGM.neighborhood
### Title: WPGM neighborhood
### Aliases: WPGM.neighborhood
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, R, lam, startb = 0) 
{
    n = nrow(X)
    p = ncol(X)
    thr = 1e-08
    maxit = 1e+06
    Xt = cbind(t(t(rep(1, n))), X)
    if (sum(startb) == 0) {
        bhat = matrix(rnorm(p + 1) * 0.01, p + 1, 1)
    }
    else {
        bhat = startb
    }
    step = 0.1
    ind = 1
    iter = 1
    while (thr < ind & iter < maxit) {
        oldb = bhat
        t = 1
        grad = wpgmGrad(Xt, Y, R, oldb)
        oldobj = wpgmObj(Xt, Y, R, oldb)
        tmp = oldb - t * grad
        bhat[1] = tmp[1]
        bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - lam * 
            t, max, 0)
        while (wpgmObj(Xt, Y, R, bhat) > oldobj - t(grad) %*% 
            (oldb - bhat) + sum((oldb - bhat)^2)/(2 * t)) {
            t = t * step
            tmp = oldb - t * grad
            bhat[1] = tmp[1]
            bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                lam * t, max, 0)
        }
        iter = iter + 1
        ind = sum((oldb - bhat)^2)/sum(oldb^2)
    }
    return(list(alpha = bhat[1], beta = bhat[-1]))
  }



cleanEx()
nameEx("WPGM.network")
### * WPGM.network

flush(stderr()); flush(stdout())

### Name: WPGM.network
### Title: Poisson network
### Aliases: WPGM.network
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, nlams, lmin = 0.001, lambda = NULL, parallel = T, 
    ncores = 4) 
{
    if (is.null(lambda)) {
        lmax = lambdaMax(t(X))
        lambda = exp(seq(log(lmax), log(lmin), l = nlams))
    }
    if (nlams != length(lambda)) {
        print("nlams is not equal to lams")
    }
    ghat = c()
    if (nlams > 0) {
        ghat = array(0, dim = c(nrow(X), nrow(X), length(lambda)))
    }
    wrapper <- function(i) {
        fit = WPGM.path.neighborhood(t(X[-i, ]), X[i, ], R, nlams, 
            lambda = lambda, 0)
        fit$beta = as.matrix(fit$Bmat)
        if (i == 1) {
            ghat[i, 2:nrow(X), ] = fit$beta
        }
        else if (i == nrow(X)) {
            ghat[i, 1:(nrow(X) - 1), ] = fit$beta
        }
        else {
            ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), ]
            ghat[i, (i + 1):nrow(X), ] = fit$beta[i:nrow(fit$beta), 
                ]
        }
        return(ghat[i, , ])
    }
    ghat2 = c()
    if (parallel) {
        library(multicore)
        ghat2 = mclapply(1:nrow(X), wrapper, mc.cores = ncores)
    }
    else {
        ghat2 = lapply(1:nrow(X), wrapper)
    }
    for (i in 1:nrow(X)) {
        ghat[i, , ] = ghat2[[i]]
    }
    ghat = lapply(1:nlams, function(r) {
        return(ghat[, , r])
    })
    return(ghat)
  }



cleanEx()
nameEx("WPGM.path.neighborhood")
### * WPGM.path.neighborhood

flush(stderr()); flush(stdout())

### Name: WPGM.path.neighborhood
### Title: WPGM neighborhood over a regularization path
### Aliases: WPGM.path.neighborhood
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, R, nlams, lmin = 0.01, lambda = NULL, startb = 0) 
{
    n = nrow(X)
    p = ncol(X)
    if (is.null(lambda)) {
        lmax = lambdaMax(t(X))
        lambda = exp(seq(log(lmax), log(lmin), l = nlams))
    }
    if (nlams == 1 & is.null(lambda)) {
        lambda = lmax
    }
    thr = 1e-08
    maxit = 1e+06
    Xt = cbind(t(t(rep(1, n))), X)
    if (sum(startb) == 0) {
        bhat = matrix(rnorm(p + 1)/p, p + 1, 1)
    }
    else {
        bhat = startb
    }
    alphas = 0
    Bmat = matrix(0, p, nlams)
    step = 0.1
    for (i in 1:nlams) {
        ind = 1
        iter = 1
        while (thr < ind & iter < maxit) {
            oldb = bhat
            t = 1
            grad = wpgmGrad(Xt, Y, R, oldb)
            oldobj = wpgmObj(Xt, Y, R, oldb)
            tmp = oldb - t * grad
            bhat[1] = tmp[1]
            bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                lambda[i] * t, max, 0)
            newobj = wpgmObj(Xt, Y, R, bhat)
            while (newobj > 9999999 | is.na(newobj) | is.na(newobj)) {
                t = t/p
                tmp = oldb - t * grad
                bhat[1] = tmp[1]
                bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                  lambda[i] * t, max, 0)
                newobj = wpgmObj(Xt, Y, R, bhat)
            }
            while (newobj > oldobj - t(grad) %*% (oldb - bhat) + 
                sum((oldb - bhat)^2)/(2 * t)) {
                t = t * step
                tmp = oldb - t * grad
                bhat[1] = tmp[1]
                bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                  lambda[i] * t, max, 0)
                newobj = wpgmObj(Xt, Y, R, bhat)
            }
            iter = iter + 1
            ind = sum((oldb - bhat)^2)
        }
        alphas[i] = bhat[1]
        Bmat[, i] = bhat[-1]
    }
    return(list(alpha = alphas, Bmat = Bmat, lambda = lambda))
  }



cleanEx()
nameEx("WPGM.select")
### * WPGM.select

flush(stderr()); flush(stdout())

### Name: WPGM.select
### Title: Winsorized Poisson Graphical Model (WPGM)
### Aliases: WPGM.select
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R = max(X), method = "star", N = 100, beta = 0.05, 
    lambda.path = NULL, nlams = 20, ncores = 4, parallel = F) 
{
    if (is.null(lambda.path)) {
        lmax = lambdaMax(t(X))
        lambda.path = exp(seq(log(lmax), log(1e-04), l = nlams))
    }
    b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
    ghat = list()
    ghat.path = list()
    ghat.path$path = vector("list", length(lambda.path))
    v = c()
    for (i in 1:N) {
        cat(paste("WPGM: Conducting sampling ... in progress: ", 
            floor(100 * (i/N)), "%", collapse = ""), "\r")
        flush.console()
        index = sample(1:ncol(X), b, replace = F)
        ghat.path$raw = WPGM.network(X[, index], R, nlams = length(lambda.path), 
            lambda = lambda.path, parallel = parallel, ncores = ncores)
        for (j in 1:length(lambda.path)) {
            tmp = ghat.path$raw[[j]]
            tmp[abs(tmp) < 1e-06] = 0
            tmp[abs(tmp) > 1e-06] = 1
            diag(tmp) = 0
            if (is.null(ghat.path$path[[j]])) {
                ghat.path$path[[j]] = tmp
            }
            else {
                ghat.path$path[[j]] = ghat.path$path[[j]] + tmp
            }
        }
    }
    for (i in 1:length(lambda.path)) {
        D = ghat.path$path[[i]]
        D = D/N
        D = 2 * D * (1 - D)
        v = c(v, mean(D[upper.tri(D)]))
    }
    v = cummax(v)
    ghat$v = v
    ghat$lambda.path = lambda.path
    ghat$opt.lambda = lambda.path[which(v == max(v[v < beta]))]
    ghat$network = WPGM.network(X, R, nlams = length(lambda.path), 
        lambda = lambda.path, parallel = T)
    ghat$opt.index = which(v == max(v[v < beta]))
    cat("\nWPGM Completed. \n")
    return(ghat)
  }



cleanEx()
nameEx("WPGMSim")
### * WPGMSim

flush(stderr()); flush(stdout())

### Name: WPGMSim
### Title: Winsorized PGM Gibbs Simulator
### Aliases: WPGMSim
### Keywords: ~kwd1 ~kwd2

### ** Examples

wpgm.sim <- WPGMSim(10, 3, 2, rep(0.5, 3), matrix(-1, 3,3))
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, p, R, alpha, Theta, maxit = 10000) 
{
    X = matrix(rpois(n * p, 1), n, p)
    iter = 1
    while (iter < maxit) {
        for (j in 1:p) {
            num = exp(matrix(1, n, 1) %*% t(alpha[j] * c(0:R) - 
                log(factorial(c(0:R)))) + matrix(c(0:R) %x% X[, 
                -j] %*% Theta[-j, j], n, R + 1))
            Pmat = num/matrix(apply(num, 1, sum), n, R + 1)
            X[, j] = apply(apply(Pmat, 1, mymult) == 1, 2, which) - 
                1
        }
        iter = iter + 1
    }
    return(X)
  }



cleanEx()
nameEx("glmpois")
### * glmpois

flush(stderr()); flush(stdout())

### Name: glmpois
### Title: Poisson based neighborhood selection
### Aliases: glmpois
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, lambda, parallel = F, nCpus = 4) 
{
    if (length(lambda) > 1) {
        ghat = array(0, dim = c(nrow(X), nrow(X), length(lambda)))
        if (parallel) {
            wrapper <- function(i) {
                fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                  lambda = lambda)
                fit$beta = as.matrix(fit$beta)
                if (ncol(fit$beta) < length(lambda)) {
                  tmp = matrix(0, nrow = nrow(fit$beta), ncol = length(lambda))
                  tmp[, 1:ncol(fit$beta)] = fit$beta
                  tmp[, ncol(fit$beta):length(lambda)] = fit$beta[, 
                    ncol(fit$beta)]
                  fit$beta = tmp
                }
                if (i == 1) {
                  ghat[i, 2:nrow(X), ] = fit$beta
                }
                else if (i == nrow(X)) {
                  ghat[i, 1:(nrow(X) - 1), ] = fit$beta
                }
                else {
                  ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), 
                    ]
                  ghat[i, (i + 1):nrow(X), ] = fit$beta[i:nrow(fit$beta), 
                    ]
                }
                return(ghat[i, , ])
            }
            library(multicore)
            ghat2 = mclapply(1:nrow(X), wrapper)
            for (i in 1:nrow(X)) {
                ghat[i, , ] = ghat2[[i]]
            }
            return(ghat)
        }
        if (parallel == F) {
            wrapper <- function(i) {
                fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                  lambda = lambda)
                fit$beta = as.matrix(fit$beta)
                if (ncol(fit$beta) < length(lambda)) {
                  tmp = matrix(0, nrow = nrow(fit$beta), ncol = length(lambda))
                  tmp[, 1:ncol(fit$beta)] = fit$beta
                  tmp[, ncol(fit$beta):length(lambda)] = fit$beta[, 
                    ncol(fit$beta)]
                  fit$beta = tmp
                }
                if (i == 1) {
                  ghat[i, 2:nrow(X), ] = fit$beta
                }
                else if (i == nrow(X)) {
                  ghat[i, 1:(nrow(X) - 1), ] = fit$bet
                }
                else {
                  ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), 
                    ]
                  ghat[i, (i + 1):nrow(X), ] = fit$beta[i:nrow(fit$beta), 
                    ]
                }
                return(ghat[i, , ])
            }
            ghat2 = lapply(1:nrow(X), wrapper)
            for (i in 1:nrow(X)) {
                ghat[i, , ] = ghat2[[i]]
            }
            return(ghat)
        }
    }
    if (length(lambda) == 1) {
        ghat = matrix(0, nrow = nrow(X), ncol = nrow(X))
        if (parallel) {
            library(snowfall)
            sfInit(cpus = nCpus)
            sfExport("X", local = T)
            sfExport("ghat", local = T)
            sfLibrary(glmnet)
            wrapper <- function(i) {
                fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                  lambda = lambda)
                fit$beta = as.numeric(fit$beta)
                if (i == 1) {
                  ghat[i, 2:nrow(X)] = fit$beta
                }
                else if (i == nrow(X)) {
                  ghat[i, 1:(nrow(X) - 1)] = fit$beta
                }
                else {
                  ghat[i, 1:(i - 1)] = fit$beta[1:(i - 1)]
                  ghat[i, (i + 1):nrow(X)] = c(fit$beta[i:length(fit$beta)])
                }
                return(ghat[i, ])
            }
            sfExport("wrapper")
            ghat = sfSapply(1:nrow(X), wrapper)
            sfStop()
            return(ghat)
        }
        for (i in 1:nrow(X)) {
            fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                lambda = lambda)
            fit$beta = as.numeric(fit$beta)
            if (i == 1) {
                ghat[i, 2:nrow(X)] = fit$beta
            }
            else if (i == nrow(X)) {
                ghat[i, 1:(nrow(X) - 1)] = fit$beta
            }
            else {
                ghat[i, 1:(i - 1)] = fit$beta[1:(i - 1)]
                ghat[i, (i + 1):nrow(X)] = c(fit$beta[i:length(fit$beta)])
            }
        }
        return(ghat)
    }
  }



cleanEx()
nameEx("lambdaMax")
### * lambdaMax

flush(stderr()); flush(stdout())

### Name: lambdaMax
### Title: Maximum lambda
### Aliases: lambdaMax
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X) 
{
    tmp = t(X) %*% X
    return(max(tmp[upper.tri(tmp)]))
  }



cleanEx()
nameEx("myglmnet.max")
### * myglmnet.max

flush(stderr()); flush(stdout())

### Name: myglmnet.max
### Title: Maximum lambda from binary search
### Aliases: myglmnet.max
### Keywords: ~kwd1

### ** Examples

library(PGM)
library(huge)
library(glmnet)
n = 200
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace =TRUE), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

lmax = myglmnet.max(t(simData))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
