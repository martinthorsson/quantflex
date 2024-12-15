# Likelihood function
f <- function(p = NULL, y = NULL, pars = NULL, trace = F, i = NULL, n_th = NULL) {

  if (i > length(pars)) stop("i out of bounds")
  # update current parameter
  pars[i] <- p

  e <- 1e-3
  # Get expression for h(y)
  hu <- eval(get.expr("y", n_th, e = e))
  hl <- eval(get.expr("y", n_th, e = -e))

  p <- (pnorm(hu) - pnorm(hl)) / (2 * e)
  if (any(p < 0)) {
    warning("Negative probabilities\n")
    p[p < 0] <- .Machine$double.eps
  }
  l <- -sum(log(p))

  if (is.infinite(l)) l <- .Machine$double.xmax

  if (trace) {
    cat("Current theta:", i - 1, "\n",
        "Likelihood: ", l, "\n",
        "Parameter estimates: ", pars[i], "\n")
    readline("Press [Enter] to continue...")
  }

  return(l)
}

# get the parameters for transformation
get.h.pars <- function(y = NULL, n_th = NULL, trace = F) {
  y <- y[!is.na(y)] # complete case

  # init parameters
  pars <- rep(0, n_th)
  names(pars) <- paste0("th_",seq(0, n_th - 1)) # name the vector

  for (i in 1:n_th) {
    fit <- nlm(f = f, p = pars[i], y = y, pars = pars,
               trace = trace, i = i, gradtol = 1e-8, n_th = n_th)
    pars[i] <- fit$estimate
  }

  if (trace) {
    cat("Parameter estimates: \n")
    for (i in 1:n_th) {
      cat("Theta", i - 1, ":", pars[i], "\n")
    }

  }
  return(pars)
}

# Gets the dynamic h(y) function depending on the number of thetas
get.expr <- function(y = NULL, n_th = NULL, e = 0) {
  s <- seq(2, n_th)
  z <- paste0(y, " + ", e, " - pars[1]")
  sym <- paste0(z, " + ", paste0("pars[", s, "]/(10^", s - 2, ")*(", z,")^", s - 1,
                                 collapse = " + "))
  expr <- parse(text = sym)
  return(expr[[1]])
}
