

utils::globalVariables(c(
  "y", "distr", "h.y", "dens", "th", "y_scaled", "tcdf"
))

#' plot
#' @description
#' Plots a flextrans object
#'
#' @param x flexTrans object
#' @param ... Additional arguments passed to the plot function
#'
#' @export
#'
#' @examples
#' d1 <- rnorm(100)
#' d2 <- rgamma(1000,1,5)
#' l <- list(d1, d2)
#' obj <- trans(l)
#' plot(obj)
plot.flxTrans <- function(x, ...) {

  df <- x$df
  pred_df <- x$pred_df
  l_org <- x$l_org

  t_list <- vector("list")
  for (i in seq_along(l_org)) {
    t_df <- data.frame(distr = paste0("distr_", i),
                       y = l_org[[i]])
    t_list[[i]] <- t_df
  }
  t_df <- do.call(rbind, t_list)

  p0 <- ggplot(t_df, aes(x = y)) +
    geom_density() +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.3, color = "lightblue") +
    facet_grid(rows = vars(distr), scales = "free") +
    labs(title = "Densities and histograms for y")


  p1 <- ggplot(df, aes(x = h.y)) +
    # geom_density() +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.3, color = "lightblue",
                   binwidth = 0.5) +
    geom_line(data = pred_df, mapping = aes(x = h.y, y = dens)) +
    facet_grid(rows = vars(distr), cols = vars(th)) +
    labs(title = "Densities and histograms for h(y)")

  # Monotonicity?
  p2 <- ggplot(df, aes(x = y_scaled, y = h.y)) +
    geom_point() +
    facet_grid(rows = vars(distr), cols = vars(th)) +
    labs(title = "Plot of y vs h(y) to verify monotonicity")

  # QQ plots of transformed y
  p3 <- ggplot(df, aes(sample = h.y)) +
    geom_qq() +
    geom_qq_line() +
    facet_grid(rows = vars(distr), cols = vars(th)) +
    labs(title = "QQ plots h(y)")

  # PP plots for h(y)

  t_df <- df |>
    group_by(distr, th) |>
    mutate(ecdf = rank(h.y) / (length(h.y) + 1) ,
           tcdf = pnorm(h.y, 0, 1))

  p4 <- ggplot(t_df, aes(x = tcdf, y = ecdf)) +
    geom_point() +
    facet_grid(rows = vars(distr), cols = vars(th)) +
    labs(title = "PP plots h(y) for N(0,1)")

  print(p0); print(p1); print(p2);print(p3);print(p4)

}
