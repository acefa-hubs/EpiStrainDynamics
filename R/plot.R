#' Generic Method for plotting
#'
#' S3 generic for plotting
#'
#' @param df Fitted model object with appropriate class
#' @param ... Additional arguments passed to methods
#' @importFrom viridis viridis
#' @importFrom stats setNames
#' @import ggplot2
#' @return Data frame with proportion analysis results
#' @export
plot <- function(df, ...) {
  UseMethod("plot")
}

#' @rdname plot
#' @export
plot.incidence <- function(df) {

  measure_df <- df$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels)), other_levels)
  )

  input_data <- data.frame(time = df$constructed_model$data$time,
                           case_timeseries = df$constructed_model$data$case_timeseries)
  ggplot(measure_df) +
    geom_line(aes(x = time, y = y, color = pathogen)) +
    geom_ribbon(aes(x = time, y = y,
                    ymin = lb_50, ymax = ub_50,
                    fill = pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = time, y = y,
                    ymin = lb_95, ymax = ub_95,
                    fill = pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    geom_point(data = input_data,
               aes(x = time, y = case_timeseries),
               size = 0.8) +
    geom_line(data = input_data,
              aes(x = time, y = case_timeseries),
              linewidth = 0.2) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")
    ) +
    ylab("Modelled influenza cases")

}

#' @rdname plot
#' @export
plot.growth_rate <- function(df) {

  measure_df <- df$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels)), other_levels)
  )

  p <- ggplot(measure_df) +
    geom_line(aes(x = time, y = y, color = pathogen)) +
    geom_ribbon(aes(x = time, y = y,
                    ymin = lb_50, ymax = ub_50,
                    fill = pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = time, y = y,
                    ymin = lb_95, ymax = ub_95,
                    fill = pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Growth rate")

  # Get the y-axis range from the built plot
  y_range <- ggplot_build(p)$layout$panel_params[[1]]$y.range

  # Calculate appropriate breaks
  n_breaks <- 5
  breaks <- pretty(y_range, n = n_breaks)

  # Create labels
  times <- log(2) / abs(breaks)
  signs <- ifelse(breaks >= 0, "", "-")
  labels <- ifelse(abs(breaks) < 1e-6, Inf,
                   paste0(signs, round(times, 1)))

  # Add secondary axis with calculated breaks and labels
  p + scale_y_continuous(
    sec.axis = sec_axis(~ .,
                        breaks = breaks,
                        labels = labels,
                        name = "Doubling(+) / Halving(-) time")
  )

}

#' @rdname plot
#' @export
plot.Rt <- function(df) {

  measure_df <- df$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels)), other_levels)
  )

  ggplot(measure_df) +
    geom_line(aes(x = time, y = y, color = pathogen)) +
    geom_ribbon(aes(x = time, y = y,
                    ymin = lb_50, ymax = ub_50,
                    fill = pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = time, y = y,
                    ymin = lb_95, ymax = ub_95,
                    fill = pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("Effective reproduction number")

}
