inf_sampler <- function(df1, sr = 100, seed = 12345) {

  set.seed(seed)

  df_new <- df1

  for (i in 1:nrow(df1)) {
    print(i)

    p1 <- c(
      rep("A", df1$inf_A[i]),
      rep("B", df1$inf_B[i]),
      rep("N", df1$num_spec[i] - df1$inf_all[i])
    )

    if (length(p1) >= sr) {
      p1_sample <- sample(p1, sr)
    } else{
      p1_sample <- p1
    }

    df_new$inf_A[i] <- length(p1_sample[p1_sample == "A"])
    df_new$inf_B[i] <- length(p1_sample[p1_sample == "B"])
    df_new$inf_all[i] <- df_new$inf_A[i] + df_new$inf_B[i]
    df_new$num_spec[i] <- min(sr, length(p1))

    p2 <- c(
      rep("H3N2", df1$inf_H3N2[i]),
      rep("H1N1", df1$inf_H1N1[i]),
      rep("US", df1$inf_A[i] - df1$inf_H3N2[i] - df1$inf_H1N1[i])
    )

    p2_sample <- sample(p2, df_new$inf_A[i])

    df_new$inf_H3N2[i] <- length(p2_sample[p2_sample == "H3N2"])
    df_new$inf_H1N1[i] <- length(p2_sample[p2_sample == "H1N1"])

  }

  df_new

}
