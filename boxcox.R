bcloglik <- function(x, lower=-1, upper=2, freq = 0.05) {
  n <- length(x)
  if (any(x <= 0, na.rm = TRUE)) {
    stop("x must be positive")
  }
  logx <- log(na.omit(c(x)))
  xdot <- exp(mean(logx))
  if (all(class(x) != "ts")) {
    fit <- lm(x ~ 1, data = data.frame(x = x), na.action = na.exclude)
  } else if (frequency(x) > 1) {
    fit <- tslm(x ~ trend + season, data = data.frame(x = x))
  } else {
    fit <- tslm(x ~ trend, data = data.frame(x = x))
  }
  xqr <- fit$qr
  lambda <- seq(lower, upper, by = freq)
  xl <- loglik <- as.vector(lambda)
  m <- length(xl)
  x <- na.omit(c(x))
  for (i in 1L:m)
  {
    if (abs(la <- xl[i]) > (2*freq)) {
      xt <- (x ^ la - 1) / la
    } else {
      xt <- logx * (1 + (la * logx) / 2 * (1 + (la * logx) / 3 * (1 + (la * logx) / 4)))
    }
    loglik[i] <- -n / 2 * log(sum(qr.resid(xqr, xt / xdot ^ (la - 1)) ^ 2))
  }
  return(loglik)
}

slambda = function(x, lower=-1, upper=2, freq = 0.05){
  lambda_seq = seq(lower, upper, by = freq)
  cum_ll = rep(0,length(lambda_seq))
  for (i in 1:dim(x)[2]){
    tmp_mat = x[,i] + matrix(x[,i]==0) * 1e-10
    cum_ll  = cum_ll + bcloglik(tmp_mat)
  }
  return(lambda_seq[which.max(cum_ll)])
}