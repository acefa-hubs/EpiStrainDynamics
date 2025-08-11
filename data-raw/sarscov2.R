# Set date limits to consider
min_date <- as.Date("2020-06-01")
max_date <- as.Date("2022-12-31")

# Load the UK COVID data from dashboard
df1a <- read.csv('data-raw/newCasesBySpecimenDate_nation_2020.csv')
df1b <- read.csv('data-raw/newCasesBySpecimenDate_nation_2021.csv')
df1c <- read.csv('data-raw/newCasesBySpecimenDate_nation_2022.csv')
df1 <- rbind(df1a, df1b, df1c)
df1 <- df1[df1$date >= min_date & df1$date <= max_date, ]

# Combine the data for each nation into one time series
df <- data.frame()
for (i in 1:length(unique(df1$date))) {
  tmp_df <- df1[df1$date %in% unique(df1$date)[i], ]

  row_df <- data.frame(date = as.Date(unique(df1$date)[i]),
                       cases = sum(tmp_df$value))
  df <- rbind(df, row_df)
}

# Load the variant data available from Lythgoe et al and format
df1S <- read.csv('data-raw/SampleInfo.csv')
df1S <- df1S[df1S$major_lineage != "Unassingned", ]

lineages_considered <- c(
  "B.1.177",
  "B.1.1.7",
  "B.1.617.2",
  "BA.1",
  "BA.2",
  "BA.2.75",
  "BA.4",
  "BA.5",
  "BQ.1",
  "XBB"
)
df1S$major_lineage <- factor(df1S$major_lineage, levels = lineages_considered)

df1S$collection_date <- as.Date(df1S$collection_date)

uniq_dates <- seq.Date(min(df1S$collection_date),
                       max(df1S$collection_date),
                       by = 1)
uniq_dates <- seq.Date(min_date, max_date, by = 1)

mat <- matrix(
  data = NA,
  nrow = length(uniq_dates),
  ncol = length(lineages_considered) + 1
)
rownames(mat) <- uniq_dates
colnames(mat) <- c(lineages_considered, "Other")

for (i in 1:length(uniq_dates)) {
  print(i)
  df_tmp <- df1S[df1S$collection_date %in% uniq_dates[i], ]
  tab <- table(df_tmp$major_lineage, useNA = "always")
  mat[i, ] <- tab
}

# Combine all data into one
df <- cbind(df, mat)
df$total <- rowSums(mat)
df$t <- as.numeric(df$date)
df$t <- df$t - min(df$t) + 1

# Include data after threshold for number of variants determined per day is
# reached
sarscov2 <- df[df$date >= as.Date("2020-09-23"), ]

usethis::use_data(sarscov2, overwrite = TRUE)
