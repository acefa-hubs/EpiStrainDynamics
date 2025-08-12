
## Loading the influenza data
df <- read.csv('data-raw/aus_influenza_data.csv')

## Selecting subset of data (2012-2023)
# Set limits on dates to consider
min_date <- as.Date("2012-01-01")
max_date <- as.Date("2020-03-01")

df$week <- lubridate::dmy(df$week)

df <- df[df$week < max_date & df$week >= min_date, ]

cols <- c('ili', 'week', 'inf_A', 'inf_B', 'inf_all',
          'inf_neg', 'inf_H3N2', 'inf_H1N1', 'num_spec')

influenza <- df[order(df$week), cols]

usethis::use_data(influenza, overwrite = TRUE)
