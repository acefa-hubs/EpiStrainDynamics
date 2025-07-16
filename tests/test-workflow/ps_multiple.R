ps_multiple_mod <- construct_model(
  method = p_spline(spline_degree = 3, days_per_knot = 5),
  pathogen_structure = multiple(
    case_timeseries = sarscov2$cases,
    component_pathogen_timeseries = list(
      B.1.177 = sarscov2$B.1.177,
      B.1.1.7 = sarscov2$B.1.1.7,
      B.1.617.2 = sarscov2$B.1.617.2,
      BA.1 = sarscov2$BA.1,
      BA.2 = sarscov2$BA.2,
      BA.2.75 = sarscov2$BA.2.75,
      BA.4 = sarscov2$BA.4,
      BA.5 = sarscov2$BA.5,
      BQ.1 = sarscov2$BQ.1,
      XBB = sarscov2$XBB,
      Other = sarscov2$Other
    ),
    smoothing_structure = 'single',
    observation_noise = 'observation_noise_only'
  ),
  dow_effect = TRUE
)
ps_multiple_fit <- fit_model(
  ps_multiple_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4
)

# PAPER VERSION
# knots <- get_knots(df$t, days_per_knot = 5, spline_degree = 3)
#
# ps_mp_data <- list(
#   num_data = nrow(df),
#   num_path = 11,
#   num_knots = length(knots),
#   knots = knots,
#   spline_degree = 3,
#   Y = df$cases,
#   X = df$t,
#   P = t(df[3:13]),
#   # Might have to fix this line
#   week_effect = 7,
#   DOW = (df$t %% 7) + 1,
#   cov_structure = 0,
#   noise_structure = 0
# )
#
# ps_mp_fit <- sampling(
#   ps_mp_mod,
#   iter = 5000,
#   warmup = 1000,
#   chains = 4,
#   data = ps_mp_data
# )
# # Get incidences
# mod_inc <- ps_incidence(
#   ps_mp_fit,
#   df$t,
#   num_days = nrow(df),
#   time_labels = df$date,
#   num_path = 11,
#   pathogen_names = c(lineages_considered, "Other")
# )
#
# # Get growth rates
# mod_gr <- ps_growth_rate(
#   ps_mp_fit,
#   df$t,
#   num_days = nrow(df),
#   time_labels = df$date,
#   num_path = 11,
#   pathogen_names = c(lineages_considered, "Other")
# )
#
# # Set mask for values not within date ranges
# mod_inc$mask <- 0
# mod_gr$mask <- 0
#
# pathogens <- c(lineages_considered, "Other")
# for (i in 1:length(pathogens)) {
#   index <- df[colnames(df) == pathogens[i]] > 0
#   min_date <- min(df[index, ]$date)
#   max_date <- max(df[index, ]$date)
#   mod_gr[mod_gr$pathogen == pathogens[i] &
#            (mod_gr$time >= min_date & mod_gr$time <= max_date), ]$mask <- 1
#   mod_inc[mod_inc$pathogen == pathogens[i] &
#             (mod_inc$time >= min_date & mod_inc$time <= max_date), ]$mask <- 1
#   print(pathogens[i])
#   #mod_gr<-mod_gr[ (!(mod_gr$pathogen==pathogens[i])) | (mod_gr$time>min_date &
#   # mod_gr$time<max_date) ,]
# }
# mod_gr[mod_gr$pathogen == "Other", ]$mask <- 0
# mod_gr[mod_gr$pathogen == "Other" &
#          (mod_gr$time <= as.Date("2021-01-19") |
#             mod_gr$time >= as.Date("2022-03-05")), ]$mask <- 1
#
# mod_gr$grp <- 0
# mod_gr[mod_gr$pathogen == "Other" &
#          mod_gr$time <= as.Date("2021-01-19"), ]$grp <- 1
#
# # Set colours
# cols <- c(brewer.pal(11, "Paired"), "grey")
# cols <- brewer.pal(12, "Paired")
# cols <- c(cols[1:10], cols[12])
#
# # Plot sub figures
# plt1 <- ggplot(mod_inc[mod_inc$pathogen == "Total", ]) +
#   geom_line(aes(x = time, y = y)) +
#   geom_ribbon(aes(
#     x = time,
#     y = y,
#     ymin = lb_50,
#     ymax = ub_50
#   ), alpha = 0.2) +
#   geom_ribbon(aes(
#     x = time,
#     y = y,
#     ymin = lb_95,
#     ymax = ub_95
#   ), alpha = 0.2) +
#   geom_point(data = df, aes(x = date, y = cases), size = 0.2) +
#   #geom_line(data=df, aes(x=date, y=cases))+
#   ylab("Cases") +
#   coord_cartesian(xlim = c(min_date + 30, max_date - 30),
#                   ylim = c(0, max(df$cases))) +
#   scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
#   scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000, 250000)) +
#   theme_bw(base_size = 14)
#
# plt2 <- ggplot(mod_inc[mod_inc$pathogen != "Total" &
#                          mod_inc$mask == 1, ]) +
#   geom_line(aes(x = time, y = y, color = pathogen)) +
#   geom_ribbon(aes(
#     x = time,
#     y = y,
#     ymin = lb_50,
#     ymax = ub_50,
#     fill = pathogen
#   ), alpha = 0.2) +
#   geom_ribbon(aes(
#     x = time,
#     y = y,
#     ymin = lb_95,
#     ymax = ub_95,
#     fill = pathogen
#   ), alpha = 0.2) +
#   coord_cartesian(xlim = c(min_date + 30, max_date - 30),
#                   ylim = c(0, max(mod_inc$y))) +
#   scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
#   scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000, 250000)) +
#   ylab("Modelled cases") +
#   scale_color_manual(values = cols) +
#   scale_fill_manual(values = cols) +
#   theme_bw(base_size = 14)
#
# plt3 <- ggplot(mod_gr[mod_gr$pathogen != "Total" &
#                         mod_gr$mask == 1, ]) +
#   geom_line(aes(
#     x = time,
#     y = y,
#     color = pathogen,
#     group = interaction(pathogen, grp)
#   )) +
#   geom_ribbon(aes(
#     x = time,
#     y = y,
#     ymin = lb_50,
#     ymax = ub_50,
#     fill = pathogen,
#     group = interaction(pathogen, grp)
#   ),
#   alpha = 0.2) +
#   geom_ribbon(aes(
#     x = time,
#     y = y,
#     ymin = lb_95,
#     ymax = ub_95,
#     fill = pathogen,
#     group = interaction(pathogen, grp)
#   ),
#   alpha = 0.2) +
#   theme_bw(base_size = 14) +
#   coord_cartesian(xlim = c(min_date + 30, max_date - 30)) +
#   scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
#   ylab("Modelled cases") +
#   scale_color_manual("Variant", values = cols) +
#   scale_fill_manual("Variant", values = cols) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   xlab("Date") +
#   scale_y_continuous("Growth rate") +
#   guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
#
# plt1i <- plt1 +
#   coord_cartesian(xlim = c(as.Date("2022-05-01"), max_date),
#                   ylim = c(0, 34000)) +
#   theme(
#     axis.title = element_blank(),
#     axis.text.x = element_blank(),
#     plot.background = element_rect(color = "black", fill = alpha("grey", 0.2))
#   ) +
#   scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000))
#
# plt2i <- plt2 +
#   coord_cartesian(xlim = c(as.Date("2022-05-01"), max_date),
#                   ylim = c(0, 8000)) +
#   theme(
#     axis.title = element_blank(),
#     legend.position = "none",
#     plot.background = element_rect(color = "black", fill = alpha("grey", 0.2))
#   ) +
#   scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000,
#                                 6000, 7000, 8000, 10000, 15000))
#
# plt1 <- plt1 +
#   labs(tag = "A") +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     plot.tag.position = c(0.01, 0.96)
#   ) +
#   annotate(
#     "rect",
#     xmin = as.Date("2022-05-01"),
#     xmax = max_date + 2,
#     ymin = -5000,
#     ymax = 35000,
#     color = "black",
#     fill = 'grey',
#     alpha = 0.2,
#     linetype = "dashed"
#   )
#
# plt2 <- plt2 +
#   labs(tag = "B") +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     legend.position = "none",
#     plot.tag.position = c(0.01, 0.96)
#   ) +
#   annotate(
#     "rect",
#     xmin = as.Date("2022-05-01"),
#     xmax = max_date + 2,
#     ymin = -5000,
#     ymax = 9000,
#     color = "black",
#     fill = "grey",
#     alpha = 0.2,
#     linetype = "dashed"
#   )
#
# plt3 <- plt3 +
#   labs(tag = "C") +
#   theme(legend.position = "bottom",
#         plot.tag.position = c(0.01, 0.96))
#
# plt1 + plt1i + plt2 + plt2i + plt3 + plot_spacer() +
#   plot_layout(ncol = 2,
#               widths = c(3, 1),
#               heights = c(2, 2, 1.5))
#
# patch1 <- plt1 + plt2 + plt3 +
#   plot_layout(nrow = 3, heights = c(2, 2, 1.5))
# patch2 <- plt1i + plt2i +
#   plot_layout(nrow = 2)
#
# cowplot::plot_grid(patch1,
#                    patch2,
#                    rel_widths = c(3, 1),
#                    rel_heights = c(1, 0.8))

