## =============================================================================
## Tubes vs 96-well plate growth-curve comparison -- FROM RAW TRIPLICATE DATA
## =============================================================================
## This is a revision of tubes_vs_96well_analysis.R. That version only had
## Avg/SD per timepoint, so it had to (a) reconstruct a Welch t-test from
## summary statistics at every timepoint, assuming n = 3, and (b) run an
## "overall" test by pairing the 15 timepoints against each other -- which
## pseudoreplicates, since OD at hour 6 is highly correlated with hour 7.
##
## Now that raw triplicate OD values are available (one row per replicate),
## this script:
##   1. Uses the real per-replicate values directly for the per-timepoint
##      Welch t-tests (no more assuming n = 3 -- it's read straight off the
##      data, and mathematically gives the identical result as before).
##   2. Replaces the pseudoreplicated "overall" test with a proper
##      replicate-level test: collapses each replicate's whole curve to one
##      number (area under the curve), then runs an unpaired Welch t-test on
##      n = 3 vs n = 3 AUC values per strain. This is the statistically
##      correct way to ask "does format affect growth, across replicates?"
##
## Uses base R only (no external packages required, except ggplot2 for the plot).
## =============================================================================

## ---- 0. Setup --------------------------------------------------------------
ALPHA      <- 0.05
input_file <- "Arab_Induction_Tu96w_Repl_Comb.csv"
output_dir <- "."

## ---- 1. Parse the raw triplicate CSV into tidy long format -----------------
## File layout: two stacked blocks ("Tubes" then "96 well plate"), each with:
##   header row: "<Label>,,<hour_0>,<hour_1>,...,,,,,"
##   then repeating groups of 3 rows per strain:
##     row 1 of group: ",<strain name>,<OD values...>"
##     rows 2-3 of group: ",,<OD values...>"           (strain name blank = same strain)
## separated by a blank line.

raw <- readLines(input_file, encoding = "UTF-8", warn = FALSE)
raw <- sub("^\uFEFF", "", raw)      # strip BOM if present
raw <- gsub("\r$", "", raw)         # strip stray carriage returns

parse_triplicate_block <- function(lines, hour_line) {
  header <- strsplit(lines[hour_line], ",", fixed = TRUE)[[1]]
  hours  <- suppressWarnings(as.numeric(header[-1]))
  hours  <- hours[!is.na(hours)]
  n_hours <- length(hours)

  i <- hour_line + 1
  cur_strain <- NA_character_
  rep_counter <- 0
  out <- vector("list", 0)

  while (i <= length(lines) && nzchar(gsub(",", "", lines[i]))) {
    cells <- strsplit(lines[i], ",", fixed = TRUE)[[1]]
    length(cells) <- max(length(cells), n_hours + 2)   # pad dropped trailing commas

    strain_cell <- trimws(cells[2])
    if (nzchar(strain_cell)) {
      cur_strain  <- strain_cell
      rep_counter <- 1
    } else {
      rep_counter <- rep_counter + 1
    }

    vals <- suppressWarnings(as.numeric(cells[3:(2 + n_hours)]))
    out[[length(out) + 1]] <- data.frame(
      strain = cur_strain, replicate = rep_counter,
      hour = hours, OD = vals, stringsAsFactors = FALSE
    )
    i <- i + 1
  }
  do.call(rbind, out)
}

tubes_hour_line <- grep("^Tubes,", raw)[1]
wells_hour_line <- grep("^96 well plate,", raw)[1]

tubes_long <- parse_triplicate_block(raw, tubes_hour_line)
tubes_long$Format <- "Tubes"

wells_long <- parse_triplicate_block(raw, wells_hour_line)
wells_long$Format <- "96-well"

# fix a stray double-underscore typo so strain names match across formats
wells_long$strain <- sub("_ATF1__ara_iAAl$", "_ATF1_ara_iAAl", wells_long$strain)

raw_long <- rbind(tubes_long, wells_long)
raw_long <- raw_long[!is.na(raw_long$OD), ]

cat("Parsed", nrow(raw_long), "raw OD readings across",
    length(unique(raw_long$strain)), "strains,",
    length(unique(raw_long$Format)), "formats,",
    "up to", max(tapply(raw_long$replicate, raw_long$Format, max)), "replicates each.\n\n")

## ---- 2. Per-timepoint Welch t-test, computed directly from raw replicates --
## Only hours present in BOTH formats can be compared.
common_hours <- intersect(unique(tubes_long$hour), unique(wells_long$hour))

pertimepoint <- do.call(rbind, lapply(unique(raw_long$strain), function(s) {
  do.call(rbind, lapply(common_hours, function(h) {
    t_vals <- raw_long$OD[raw_long$strain == s & raw_long$Format == "Tubes"   & raw_long$hour == h]
    w_vals <- raw_long$OD[raw_long$strain == s & raw_long$Format == "96-well" & raw_long$hour == h]
    if (length(t_vals) < 2 || length(w_vals) < 2) return(NULL)
    tt <- t.test(t_vals, w_vals)
    data.frame(
      strain = s, hour = h,
      avg_t = mean(t_vals), sd_t = sd(t_vals),
      avg_w = mean(w_vals), sd_w = sd(w_vals),
      t_stat = unname(tt$statistic), df = unname(tt$parameter), p_raw = tt$p.value
    )
  }))
}))

pertimepoint$p_adj_strain  <- ave(pertimepoint$p_raw, pertimepoint$strain,
                                   FUN = function(x) p.adjust(x, method = "BH"))
pertimepoint$p_adj_overall <- p.adjust(pertimepoint$p_raw, method = "BH")
pertimepoint$significant_raw <- pertimepoint$p_raw < ALPHA
pertimepoint$significant_adj <- pertimepoint$p_adj_overall < ALPHA
pertimepoint <- pertimepoint[order(pertimepoint$strain, pertimepoint$hour), ]

cat("=== Per-timepoint Welch t-test (Tubes vs 96-well), computed from raw n = 3 replicates ===\n")
print(pertimepoint, digits = 3, row.names = FALSE)

write.csv(pertimepoint, file.path(output_dir, "tubes_vs_96well_pertimepoint_ttests_fromraw.csv"),
          row.names = FALSE)

## ---- 3. Replicate-level AUC test (replaces the pseudoreplicated "overall" test) --
## Collapses each replicate's whole curve to one number (trapezoidal AUC),
## then runs an unpaired Welch t-test on n = 3 vs n = 3 AUC values per strain.
auc_trapz <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

auc_by_rep <- do.call(rbind, lapply(split(raw_long, list(raw_long$strain, raw_long$Format, raw_long$replicate),
                                            drop = TRUE), function(d) {
  data.frame(strain = d$strain[1], Format = d$Format[1], replicate = d$replicate[1],
             AUC = auc_trapz(d$hour, d$OD))
}))

replicate_auc_test <- do.call(rbind, lapply(unique(auc_by_rep$strain), function(s) {
  t_auc <- auc_by_rep$AUC[auc_by_rep$strain == s & auc_by_rep$Format == "Tubes"]
  w_auc <- auc_by_rep$AUC[auc_by_rep$strain == s & auc_by_rep$Format == "96-well"]
  tt <- t.test(t_auc, w_auc)   # unpaired Welch, n = 3 vs n = 3
  data.frame(
    strain = s,
    mean_AUC_tubes = mean(t_auc), mean_AUC_96well = mean(w_auc),
    t_stat = unname(tt$statistic), df = unname(tt$parameter), p_raw = tt$p.value
  )
}))
replicate_auc_test$p_BH <- p.adjust(replicate_auc_test$p_raw, method = "BH")
replicate_auc_test <- replicate_auc_test[order(replicate_auc_test$p_BH), ]

cat("\n=== Replicate-level AUC test per strain (n = 3 vs n = 3, unpaired Welch) ===\n")
cat("(Replaces the earlier timepoint-paired test, which pseudoreplicated across\n")
cat(" correlated hours within a growth curve.)\n")
print(replicate_auc_test, digits = 3, row.names = FALSE)

write.csv(replicate_auc_test, file.path(output_dir, "tubes_vs_96well_replicate_AUC_ttests.csv"),
          row.names = FALSE)

## ---- 4. Per-strain summary: how much of the curve differs? -----------------
strain_summary <- do.call(rbind, lapply(split(pertimepoint, pertimepoint$strain), function(d) {
  data.frame(
    strain = d$strain[1],
    n_timepoints = nrow(d),
    n_sig_raw_p05 = sum(d$significant_raw, na.rm = TRUE),
    n_sig_BH_p05  = sum(d$significant_adj, na.rm = TRUE),
    mean_diff_tubes_minus_96w = mean(d$avg_t - d$avg_w, na.rm = TRUE)
  )
}))
strain_summary <- strain_summary[order(-strain_summary$n_sig_BH_p05), ]

cat("\n=== Per-strain summary across all comparable timepoints ===\n")
print(strain_summary, digits = 3, row.names = FALSE)
write.csv(strain_summary, file.path(output_dir, "tubes_vs_96well_strain_summary_fromraw.csv"),
          row.names = FALSE)

## ---- 5. Plot: growth curves with error bars, by strain (ggplot2) ----------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

summary_agg <- aggregate(OD ~ strain + hour + Format, data = raw_long,
                          FUN = function(x) c(avg = mean(x), sd = sd(x)))
summary_agg <- do.call(data.frame, summary_agg)
names(summary_agg)[4:5] <- c("avg", "sd")

strain_order <- c(
  "EPI300", "EPI300_ara", "EPI300_pCC1", "EPI300_pCC1_ara",
  "EPI300_pCC1_ATF1", "EPI300_pCC1_ATF1_ara",
  "EPI300_pCC1_ATF1_iAAl", "EPI300_pCC1_ATF1_ara_iAAl"
)
summary_agg$strain <- factor(summary_agg$strain, levels = strain_order)

p <- ggplot(summary_agg, aes(x = hour, y = avg, color = Format, fill = Format)) +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 0.3, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.0) +
  facet_wrap(~ strain, scales = "free_y") +
  labs(
    title = "Growth curves using triplicates: Tubes vs 96-well plate",
    x = "Time (h)", y = "OD (Avg \u00b1 SD)",
    color = "Format", fill = "Format"
  ) +
  scale_color_manual(values = c("Tubes" = "#1b6ca8", "96-well" = "#d9822b")) +
  scale_fill_manual(values = c("Tubes" = "#1b6ca8", "96-well" = "#d9822b")) +
  theme_bw(base_size = 11)

print(p)

ggsave(file.path(output_dir, "tubes_vs_96well_growthcurves_fromraw.png"),
       plot = p, width = 11, height = 7, dpi = 150)

cat("\nSaved: tubes_vs_96well_pertimepoint_ttests_fromraw.csv\n")
cat("Saved: tubes_vs_96well_replicate_AUC_ttests.csv\n")
cat("Saved: tubes_vs_96well_strain_summary_fromraw.csv\n")
cat("Saved: tubes_vs_96well_growthcurves_fromraw.png\n")

## ---- 6. Overall conclusion --------------------------------------------------
cat("\n=== Bottom line ===\n")
n_strains_sig <- sum(replicate_auc_test$p_BH < ALPHA)
cat(n_strains_sig, "out of", nrow(replicate_auc_test),
    "strains show a statistically significant (BH-adjusted p <", ALPHA,
    ") difference in whole-curve AUC between Tubes and 96-well plate,\n")
cat("using the proper replicate-level test (n = 3 vs n = 3) rather than\n")
cat("pseudoreplicated timepoints.\n")
