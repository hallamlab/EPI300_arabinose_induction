# ============================================================
# Arabinose induction: control vs treatment (+ara) at each time point
# Welch t-test (n = 3 vs 3) per hour, 0-14 h, with Holm correction
# Input : Arab_Induction_replicates.csv  (raw replicates, 3 rows per strain)
# Output: Arab_Induction_pvalues_0-14h.csv
# ============================================================

library(tidyverse)

# ---- 0. Settings -----------------------------------------------------------
infile   <- "Arab_Induction_replicates.csv"   # put in working dir or give full path
outfile  <- "Arab_Induction_pvalues_0-14h.csv"
max_time <- 14                                # analyse time points up to this hour
alpha    <- 0.05

# Control / treatment pairs (names must match the first column of the CSV)
pairs <- tribble(
  ~Control,                 ~Treatment,
  "EPI300",                 "EPI300_ara",
  "EPI300_pCC1",            "EPI300_pCC1_ara",
  "EPI300_pCC1_ATF1",       "EPI300_pCC1_ATF1_ara",
  "EPI300_pCC1_ATF1_iAAl",  "EPI300_pCC1_ATF1_iAAl_ara"
)

# ---- 1. Read the raw replicate file ---------------------------------------
# Layout of the file:
#   Row 1: header ("h", "Average" ...)
#   Row 2: time points (0,1,...,14,20) in every 2nd column (columns 2,4,6,...)
#   Row 3+: strain name in column 1 on the first of 3 replicate rows;
#           replicate OD600 values in columns 2,4,6,...
# The Average/SD columns (3,5,7,...) are ignored: they are recalculated here.
raw <- read_csv(infile, col_names = FALSE, col_types = cols(.default = "c"),
                show_col_types = FALSE)

n_time   <- (ncol(raw) - 1) %/% 2           # number of time points
val_cols <- seq(2, by = 2, length.out = n_time)

times <- as.numeric(unlist(raw[2, val_cols]))

data_rows <- raw[-(1:2), ]
strain <- data_rows[[1]]
strain[strain == ""] <- NA
for (i in seq_along(strain)) {              # fill strain name down the 3 rows
  if (is.na(strain[i]) && i > 1) strain[i] <- strain[i - 1]
}
rep_id <- ave(seq_along(strain), strain, FUN = seq_along)

# ---- 2. Long format: Strain, Replicate, Time, OD600 ------------------------
df_long <- map_dfr(seq_len(n_time), function(j) {
  tibble(
    Strain    = strain,
    Replicate = rep_id,
    Time      = times[j],
    OD600     = as.numeric(data_rows[[val_cols[j]]])
  )
}) |>
  filter(!is.na(OD600), Time <= max_time)

# ---- 3. Welch t-test at each time point for each pair ----------------------
one_test <- function(control, treatment, time) {
  a <- df_long |> filter(Strain == control,   Time == time) |> pull(OD600)
  b <- df_long |> filter(Strain == treatment, Time == time) |> pull(OD600)
  tt <- t.test(a, b, var.equal = FALSE)     # Welch (unequal variances)
  tibble(
    Control_mean_OD600   = mean(a),
    Control_SD           = sd(a),
    Treatment_mean_OD600 = mean(b),
    Treatment_SD         = sd(b),
    Difference_OD600     = mean(b) - mean(a),
    Difference_percent   = 100 * (mean(b) / mean(a) - 1),
    Raw_p_value          = tt$p.value
  )
}

results <- pairs |>
  crossing(Time_h = sort(unique(df_long$Time))) |>
  mutate(stats = pmap(list(Control, Treatment, Time_h), one_test)) |>
  unnest(stats) |>
  group_by(Control) |>                      # Holm correction within each strain
  mutate(Holm_adjusted_p = p.adjust(Raw_p_value, method = "holm")) |>
  ungroup() |>
  mutate(
    `Significant_raw_p<0.05`  = if_else(Raw_p_value     < alpha, "yes", "no"),
    `Significant_Holm_p<0.05` = if_else(Holm_adjusted_p < alpha, "yes", "no")
  ) |>
  arrange(match(Control, pairs$Control), Time_h) |>
  mutate(
    across(c(Control_mean_OD600, Control_SD, Treatment_mean_OD600,
             Treatment_SD, Difference_OD600), ~ round(.x, 3)),
    Difference_percent = round(Difference_percent, 1),
    across(c(Raw_p_value, Holm_adjusted_p), ~ round(.x, 4))
  )

print(results, n = Inf, width = Inf)

# ---- 4. Wide table of raw p-values (time x strain) -------------------------
p_table <- results |>
  select(Time_h, Control, Raw_p_value) |>
  pivot_wider(names_from = Control, values_from = Raw_p_value)
print(p_table, n = Inf)

# ---- 5. First time point where each strain becomes significant -------------
first_sig <- results |>
  group_by(Control) |>
  summarise(
    First_time_raw_p_below_alpha  = suppressWarnings(min(Time_h[Raw_p_value     < alpha])),
    First_time_Holm_p_below_alpha = suppressWarnings(min(Time_h[Holm_adjusted_p < alpha])),
    .groups = "drop"
  ) |>
  mutate(across(starts_with("First"), ~ if_else(is.infinite(.x), NA_real_, .x)))
print(first_sig)

# ---- 6. Save ----------------------------------------------------------------
write_csv(results, outfile)
message("Saved: ", outfile)
