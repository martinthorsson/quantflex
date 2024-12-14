# Header ----------------------------------
# Purpose: Transformations using Taylor approximations
# Date: 2024-12-13
# Author: Martin Thorsson

#' trans
#' @description
#'  Transforms a vector or a list of vectors to a gaussian distribution using a Taylor-like transformation function.
#' @param l A list with vectors or a vector to transform
#' @param n_theta Number of theta parameters to use in the transformation
#' @param plot Plot original and transformed distribition together with various diagnostic plots.
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
trans <- function(l, n_theta = 3, plot = F) {
  n <- sapply(l, length)

  if (plot) {

    df <- data.frame()

    for (i in seq_along(n)) {
      t_df <- data.frame(n = n[i],
                         y = l[[i]])
      df <- rbind(df, t_df)
    }

    ggplot(df, aes(x = y)) +
      geom_density() +
      geom_histogram(aes(y = after_stat(density)), alpha = 0.3, color = "lightblue") +
      facet_grid(rows = vars(n)) +
      labs(title = "Densities and histograms for y")

  }

  l_org <- l

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
  j <- 1; i <- 3; th <- 3
  for (th in seq_along(n_theta)) {

    for (i in seq_along(n)) {
      y <- l[[i]]
      pars <- pars_list[[th]][[i]]

      h <- eval(get.expr("y", n_theta[th]))
      t_df <- data.frame(n = n[i],
                         th = n_theta[th],
                         h = h,
                         y_org = l_org[[i]],
                         y_scaled = l[[i]])

      t_list[[j]] <- t_df
      j <- j + 1
    }
  }
  df <- do.call(rbind, t_list)

  t_list <- vector("list")
  j <- 1
  for (th in seq_along(n_theta)) {

    for (i in seq_along(n)) {
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

      t_df <- data.frame(n = n[i],
                         th = n_theta[th],
                         h = h,
                         dens = dens)

      t_list[[j]] <- t_df
      j <- j + 1
    }
  }
  pred_df <- do.call(rbind, t_list)

  if (plot) {
    # Density plots för h(y)
    p1 <- ggplot(df, aes(x = h)) +
      # geom_density() +
      geom_histogram(aes(y = after_stat(density)), alpha = 0.3, color = "lightblue",
                     binwidth = 0.5) +
      geom_line(data = pred_df, mapping = aes(x = h, y = dens)) +
      facet_grid(rows = vars(n), cols = vars(th)) +
      labs(title = "Densities and histograms for h(y)")

    # Monotonicity?
    p2 <- ggplot(df, aes(x = y_scaled, y = h)) +
      geom_point() +
      facet_grid(rows = vars(n), cols = vars(th)) |>
      labs(title = "Plot of y vs h(y) to verify monotonicity")

    # QQ plots of transformed y
    p3 <- ggplot(df, aes(sample = h)) +
      geom_qq() +
      geom_qq_line() +
      facet_grid(rows = vars(n), cols = vars(th)) +
      labs(title = "QQ plots h(y)")

    # PP plots for h(y)

    t_df <- df |>
      group_by(n, th) |>
      mutate(ecdf = rank(h) / (length(h) + 1) ,
             tcdf = pnorm(h, 0, 1))

    p4 <- ggplot(t_df, aes(x = tcdf, y = ecdf)) +
      geom_point() +
      facet_grid(rows = vars(n), cols = vars(th)) +
      labs(title = "PP plots h(y) for N(0,1)")

    print(p1); print(p2);print(p3);print(p4)

  }

  return(df)
}



