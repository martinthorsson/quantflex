# Header ----------------------------------
# Purpose: Transformations using Taylor approximations
# Date: 2024-12-13
# Author: Martin Thorsson

#' trans
#' @description
#'  Transforms a vector or a list of vectors to a gaussian distribution using a Taylor-like transformation function.
#' @param l A list with vectors or a vector to transform
#' @param n_theta A integer specifying the number of theta parameters to use in the transformation.
#' If a vector of integers is supplied, there will be one transformation per element.
#'
#' @return Returns a data frame with the transformed vectors in long format.
#' @export
#'
#' @examples
#' x1 <- rgamma(100, 1, 5)
#' x2 <- rgamma(1000, 1, 5)
#' x3 <- rbeta(100, 2, 4)
#' l <- list(x1, x2, x3)
#' trans(l, n_theta = 3)
trans <- function(l, n_theta = 3) {

  if (is.numeric(l)) l <- list(l)

  l_org <- l # save original list

  # logit transform

  l <- lapply(l_org, \(y) {
    y_norm <- (y - min(y)) / ( max(y) - min(y) )
    y_norm[y_norm == 0] <- NA
    logit <- log(y_norm / (1 - y_norm))
    logit[is.infinite(logit)] <- NA
    return(logit)
  })


  # estimate transformation parameters
  pars_list <- lapply(n_theta, \(n_th) {
    lapply(l, \(y) {
      get.h.pars(y, n_th = n_th, trace = F)
    })
  })

  # compute h(y) and bind to data frame
  t_list <- vector("list")
  j <- 1
  for (th in seq_along(n_theta)) {

    for (i in seq_along(l)) {
      y <- l[[i]]
      pars <- pars_list[[th]][[i]]

      h <- eval(get.expr("y", n_theta[th]))
      t_df <- data.frame(distr = paste0("distr_", i),
                         th = paste0("theta_",n_theta[th]),
                         h.y = h,
                         y_org = l_org[[i]],
                         y_scaled = y)

      t_list[[j]] <- t_df
      j <- j + 1
    }
  }
  df <- do.call(rbind, t_list)

  t_list <- vector("list")
  j <- 1

  for (th in seq_along(n_theta)) {

    for (i in seq_along(l)) {
      pars <- pars_list[[th]][[i]]

      y_min <- min(l[[i]], na.rm = T)
      y_max <- max(l[[i]], na.rm = T)
      new_y <- seq(y_min, y_max, length.out = 100)
      e <- 1e-3
      # Get expression for h(y)
      h <- eval(get.expr("new_y", n_theta[th]))
      hu <- eval(get.expr("new_y", n_theta[th], e = e))
      hl <- eval(get.expr("new_y", n_theta[th], e = -e))

      dens <- (pnorm(hu) - pnorm(hl)) / (2 * e)

      t_df <- data.frame(distr = paste0("distr_", i),
                         th = paste0("theta_",n_theta[th]),
                         h.y = h,
                         dens = dens)

      t_list[[j]] <- t_df
      j <- j + 1
    }
  }

  pred_df <- do.call(rbind, t_list)

  return_list         <- vector("list")
  return_list$df      <- df
  return_list$pred_df <- pred_df
  return_list$l_org   <- l_org
  return_list         <- structure(return_list, class = "flxTrans")

  return(return_list)
}



