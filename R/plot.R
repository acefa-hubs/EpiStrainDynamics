#' Generic Method for plotting metrics calculation outputs
#'
#' S3 generic for plotting
#'
#' @param df Metrics calculation output of class `EpiStrainDynamics.metric`
#'  from either `incidence()`, `growth_rate()`, `Rt()`, or `proportion()`.
#' @param xlab Time label for x axis, defaults to "Time"
#' @param ... Additional arguments passed to plot
#' @importFrom viridis viridis
#' @importFrom stats setNames
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @return ggplot2 plot output
#' @export
#'
#' @srrstats {BS6.1} All calculated epidemiological metrics can be plotted with
#'   a default plot function.
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {TS5.0} default plot methods implemented
#' @srrstats {TS5.1, TS5.2, TS5.3} Time on x axis with units printed
#'
#' @examples
#' \dontrun{
#'   mod <- construct_model(
#'     method = random_walk(),
#'     pathogen_structure = single(
#'       case_timeseries = sarscov2$cases,
#'       time = sarscov2$date))
#'
#'   fit <- fit_model(mod)
#'   gr <- growth_rate(mod)
#'   plot(gr)
#' }
plot <- function(df, xlab = 'Time', ...) {
  validate_class_inherits(df, 'EpiStrainDynamics.metric')
  UseMethod("plot")
}

#' @rdname plot
#' @export
plot.incidence <- function(df, xlab = 'Time', ...) {

  measure_df <- df$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels)), other_levels)
  )

  tsbl <- inc$constructed_model$validated_tsbl
  time_col <- tsibble::index_var(tsbl)
  input_data <- tsbl[, c(time_col, 'case_timeseries')]

  ggplot(measure_df) +
    geom_line(aes(x = .data$time, y = .data$y, color = .data$pathogen)) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_50, ymax = .data$ub_50,
                    fill = .data$pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_95, ymax = .data$ub_95,
                    fill = .data$pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    geom_point(data = input_data,
               aes(x = .data[[time_col]], y = .data$case_timeseries),
               size = 0.8) +
    geom_line(data = input_data,
              aes(x = .data[[time_col]], y = .data$case_timeseries),
              linewidth = 0.2) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")
    ) +
    ylab("Modelled influenza cases") +
    theme(legend.title = element_blank()) +
    xlab(xlab)

}

#' @rdname plot
#' @export
plot.growth_rate <- function(df, xlab = 'Time', ...) {

  measure_df <- df$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels)), other_levels)
  )

  p <- ggplot(measure_df) +
    geom_line(aes(x = .data$time, y = .data$y, color = .data$pathogen)) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_50, ymax = .data$ub_50,
                    fill = .data$pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_95, ymax = .data$ub_95,
                    fill = .data$pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Growth rate") +
    theme(legend.title = element_blank()) +
    xlab(xlab)

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
plot.Rt <- function(df, xlab = 'Time', ...) {

  measure_df <- df$measure

  total_pathogens <- "Total"
  all_levels <- unique(measure_df$pathogen)
  other_levels <- setdiff(all_levels, total_pathogens)

  colors <- c(
    setNames("black", total_pathogens),
    setNames(viridis::viridis(length(other_levels)), other_levels)
  )

  ggplot(measure_df) +
    geom_line(aes(x = .data$time, y = .data$y, color = .data$pathogen)) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_50, ymax = .data$ub_50,
                    fill = .data$pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_95, ymax = .data$ub_95,
                    fill = .data$pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("Effective reproduction number") +
    theme(legend.title = element_blank()) +
    xlab(xlab)

}

#' @rdname plot
#' @export
plot.proportion <- function(df, xlab = 'Time', ...) {

  measure_df <- df$measure

  combos <- unique(measure_df$pathogen)
  colors <- setNames(viridis::viridis(length(combos)), combos)

  ggplot(measure_df) +
    geom_line(aes(x = .data$time, y = .data$y, color = .data$pathogen)) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_50, ymax = .data$ub_50,
                    fill = .data$pathogen),
                alpha = 0.2) +
    geom_ribbon(aes(x = .data$time, y = .data$y,
                    ymin = .data$lb_95, ymax = .data$ub_95,
                    fill = .data$pathogen),
                alpha = 0.2) +
    theme_bw(base_size = 14) +
    scale_colour_manual(
      values = colors,
      aesthetics = c("colour", "fill")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    ylab("Modelled proportion of cases") +
    theme(legend.title = element_blank()) +
    xlab(xlab)

}
