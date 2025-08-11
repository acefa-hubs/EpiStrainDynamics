df <- as.data.frame(rw_subtyped_fit$constructed_model$data)

first_date <- as.Date("2011-08-29")
origin_date <- as.Date("2020-03-01")
date_labels <- seq(first_date, origin_date, length.out = nrow(df))
df$date <- date_labels

rw_subtyped_incidence <- rw_subtyped_incidence |>
  dplyr::group_by(pathogen) |>
  dplyr::mutate(time = dplyr::row_number(),
                date = date_labels)

rw_subtyped_gr <- growth_rate(rw_subtyped_fit) |>
  dplyr::group_by(pathogen) |>
  dplyr::mutate(time = dplyr::row_number(),
                date = date_labels[2:length(date_labels)])

rw_subtyped_rt <- Rt(rw_subtyped_fit,
                     gi_dist = rw_gi_dist) |>
  dplyr::group_by(pathogen) |>
  dplyr::mutate(time = dplyr::row_number(),
                date = date_labels[7:length(date_labels)])

inc_plot <- ggplot(rw_subtyped_incidence) +
  geom_line(aes(x = date, y = y, color = pathogen)) +
  geom_ribbon(aes(x = date, y = y, ymin = lb_50, ymax = ub_50, fill = pathogen), alpha = 0.2) +
  geom_ribbon(aes(x = date, y = y, ymin = lb_95, ymax = ub_95, fill = pathogen), alpha = 0.2) +
  theme_bw(base_size = 14) +
  geom_point(data = df, aes(x = date, y = case_timeseries), size = 0.8) +
  geom_line(data = df, aes(x = date, y = case_timeseries), linewidth = 0.2) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none") +
  ylab("Modelled influenza cases") +
  xlab("Date")
  # coord_cartesian(xlim=c(first_date, origin_date), ylim=c(0,max(mod_inc[mod_inc$time> (max(mod_inc$time)-180),]$ub_95, df[df$notification_date> (max(df$notification_date)-180),]$cases )))+
  # scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") #+
  # geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Influenza cases"),hjust=-0,vjust=1.2, fill = "white", size=5)+
  # theme( axis.text.x=element_blank(),
  #        axis.title.x = element_blank())

gr_plot <- ggplot(rw_subtyped_gr) +
  geom_line(aes(x = date, y = y, color = pathogen)) +
  geom_ribbon(aes(x = date, y = y, ymin = lb_50, ymax = ub_50, fill = pathogen), alpha = 0.2) +
  geom_ribbon(aes(x = date, y = y, ymin = lb_95, ymax = ub_95, fill = pathogen), alpha = 0.2) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(
    "Growth rate",
    sec.axis = sec_axis(~., breaks = c(log(2)/7, log(2)/14, log(2)/21,0,
                                       -log(2)/21, -log(2)/14, -log(2)/7),
                        labels = c( "7", "14","21",
                                    expression(infinity/-infinity),
                                    "-21", "-14", "-7"),
                        name = "Doubling(+) / Halving(-) time (days)")) +
  xlab("Date") #+
# coord_cartesian(ylim=c(-0.12,0.12), xlim=c(first_date, origin_date))+
  # scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")+
  # geom_label(data=data.frame(), aes(x = first_date, y = Inf, label ="Growth rate"),hjust=-0,vjust=1.2, fill = "white", size=5)+
  # theme(legend.position = "none",
  #       axis.text.x=element_blank(),
  #       axis.title.x = element_blank())

rt_plot <- ggplot(rw_subtyped_rt) +
  geom_line(aes(x = date, y = y, color = pathogen)) +
  geom_ribbon(aes(x = date, y = y, ymin = lb_50, ymax = ub_50, fill = pathogen), alpha = 0.2) +
  geom_ribbon(aes(x = date, y = y, ymin = lb_95, ymax = ub_95, fill = pathogen), alpha = 0.2) +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer("Sub-type", palette = "Dark2") +
  scale_color_brewer("Sub-type", palette = "Dark2") +
  theme(legend.position = "bottom",
        axis.title.x = element_blank()) +
  ylab("Effective reproduction number") +
  xlab("Date")
  # coord_cartesian(xlim = c(first_date, origin_date)) +
  # geom_label(data = data.frame(), aes(x = first_date, y = Inf, label = paste(expression(R[t])) ),parse=T,hjust=-0,vjust=1.2, fill = "white", size=5)+
  # scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y")
