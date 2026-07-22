## =============================================================================
## Tubes vs 96-well plate growth-curve comparison
## =============================================================================
## The source file only reports Avg and SD per timepoint (no raw replicate
## values), so a classic two-sample t-test can't be run directly on raw data.
## Instead this script reconstructs a Welch two-sample t-test from the
## reported summary statistics (mean, SD, n) at every strain x timepoint
## combination where both formats were measured.
##
## IMPORTANT ASSUMPTION: n = 3 replicates for BOTH tubes and 96-well cultures
## (set below). Change N_TUBES / N_96W if your actual replicate count differs.
##
## Uses base R only (no external packages required).
## =============================================================================

## ---- 0. Setup --------------------------------------------------------------
N_TUBES <- 3   # number of replicates used to calculate each Tubes Avg/SD
N_96W   <- 3   # number of replicates used to calculate each 96-well Avg/SD
ALPHA   <- 0.05
input_file  <- "Arab_Induction_Tu96w_Comb.csv"
output_dir  <- "."

## ---- 1. Parse the (irregular) CSV into tidy long format --------------------
## The file has two stacked blocks ("Tubes" then "96-well plate"), each with:
##   line 1: block label, then strain name repeated every 2 columns
##   line 2: "h,Avg,SD,Avg,SD,..."
##   data rows: hour, avg1, sd1, avg2, sd2, ...
## separated by a blank line. We parse this manually rather than relying on
## read.csv, since the header spans two rows and columns are unlabeled.

raw <- readLines(input_file, encoding = "UTF-8", warn = FALSE)
raw <- sub("^\uFEFF", "", raw)      # strip BOM if present
raw <- gsub("\r$", "", raw)         # strip stray carriage returns

parse_block <- function(lines, start) {
  header  <- strsplit(lines[start], ",")[[1]]
  fmt     <- header[1]
  strains <- header[-1]
  strains <- strains[nzchar(trimws(strains))]

  i <- start + 2   # skip block label row + "h,Avg,SD,..." row
  out <- vector("list", 0)
  while (i <= length(lines) && nzchar(gsub(",", "", lines[i]))) {
    cells <- strsplit(lines[i], ",")[[1]]
    hour  <- as.numeric(cells[1])
    vals  <- as.numeric(cells[-1])
    for (j in seq_along(strains)) {
      avg <- vals[2 * j - 1]
      sd  <- vals[2 * j]
      out[[length(out) + 1]] <- data.frame(
        format = fmt, strain_raw = strains[j], hour = hour, avg = avg, sd = sd,
        stringsAsFactors = FALSE
      )
    }
    i <- i + 1
  }
  list(df = do.call(rbind, out), next_line = i)
}

block1 <- parse_block(raw, 1)
next_start <- block1$next_line
while (next_start <= length(raw) && !nzchar(gsub(",", "", raw[next_start]))) {
  next_start <- next_start + 1
}
block2 <- parse_block(raw, next_start)

df <- rbind(block1$df, block2$df)
df$strain <- sub("_T$", "", df$strain_raw)
df$strain <- sub("_96w$", "", df$strain)
df$format <- ifelse(df$format == "Tubes", "Tubes", "96-well")
df <- df[, c("format", "strain", "hour", "avg", "sd")]

cat("Parsed", nrow(df), "rows across", length(unique(df$strain)), "strains and",
    length(unique(df$format)), "formats.\n\n")

## ---- 2. Pair Tubes vs 96-well at matching strain x hour --------------------
## Tubes were sampled at hours 0-14 & 20; 96-well at hours 0-19 & 23.
## Only hours present in BOTH formats can be compared -> merge keeps 0-14.
tubes <- df[df$format == "Tubes", ]
w96   <- df[df$format == "96-well", ]
names(tubes)[names(tubes) %in% c("avg", "sd")] <- c("avg_t", "sd_t")
names(w96)[names(w96)   %in% c("avg", "sd")] <- c("avg_w", "sd_w")

paired <- merge(
  tubes[, c("strain", "hour", "avg_t", "sd_t")],
  w96[,   c("strain", "hour", "avg_w", "sd_w")],
  by = c("strain", "hour")
)
paired <- paired[order(paired$strain, paired$hour), ]

cat("Comparable timepoints per strain (present in both formats): hours",
    paste(sort(unique(paired$hour)), collapse = ", "), "\n\n")

## ---- 3. Welch two-sample t-test reconstructed from summary statistics -----
welch_from_summary <- function(m1, s1, n1, m2, s2, n2) {
  se <- sqrt(s1^2 / n1 + s2^2 / n2)
  if (se == 0) return(c(t = NA, df = NA, p = NA))
  t  <- (m1 - m2) / se
  dfree <- se^4 / ( (s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1) )
  p  <- 2 * pt(-abs(t), dfree)
  c(t = t, df = dfree, p = p)
}

stat_mat <- t(mapply(welch_from_summary,
                      paired$avg_t, paired$sd_t, N_TUBES,
                      paired$avg_w, paired$sd_w, N_96W))
results <- cbind(paired, as.data.frame(stat_mat))
names(results)[names(results) %in% c("t", "df", "p")] <- c("t_stat", "df", "p_raw")

results$p_adj_strain <- ave(results$p_raw, results$strain,
                            FUN = function(x) p.adjust(x, method = "BH"))
results$p_adj_overall <- p.adjust(results$p_raw, method = "BH")
results$significant_raw <- results$p_raw < ALPHA
results$significant_adj <- results$p_adj_overall < ALPHA
results <- results[order(results$strain, results$hour), ]

cat("=== Per-timepoint Welch t-test (Tubes vs 96-well), n =", N_TUBES,
    "vs", N_96W, "assumed ===\n")
print(results, digits = 3, row.names = FALSE)

write.csv(results, file.path(output_dir, "tubes_vs_96well_pertimepoint_ttests.csv"),
          row.names = FALSE)

## ---- 4. Per-strain summary: how much of the curve differs? ----------------
strain_summary <- do.call(rbind, lapply(split(results, results$strain), function(d) {
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
write.csv(strain_summary, file.path(output_dir, "tubes_vs_96well_strain_summary.csv"),
          row.names = FALSE)

## ---- 5. Overall (whole-curve) paired t-test per strain ---------------------
## Complementary check: treats each comparable hour as a paired observation
## and asks whether, overall, one format reads systematically higher than the
## other across the growth curve for that strain (paired t-test on the Avg
## values themselves; doesn't require the N_TUBES/N_96W assumption).
overall_paired <- do.call(rbind, lapply(split(paired, paired$strain), function(d) {
  tt <- t.test(d$avg_t, d$avg_w, paired = TRUE)
  data.frame(
    strain = d$strain[1],
    paired_t = unname(tt$statistic),
    paired_df = unname(tt$parameter),
    paired_p = tt$p.value,
    mean_diff = mean(d$avg_t - d$avg_w)
  )
}))
overall_paired$paired_p_BH <- p.adjust(overall_paired$paired_p, method = "BH")
overall_paired <- overall_paired[order(overall_paired$paired_p_BH), ]

cat("\n=== Overall paired t-test per strain (Tubes vs 96-well across the whole curve) ===\n")
print(overall_paired, digits = 3, row.names = FALSE)
write.csv(overall_paired, file.path(output_dir, "tubes_vs_96well_overall_paired.csv"),
          row.names = FALSE)

## ---- 6. Plot: growth curves with error bars, by strain (ggplot2) ----------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

plot_df <- rbind(
  data.frame(strain = tubes$strain, hour = tubes$hour, avg = tubes$avg_t,
             sd = tubes$sd_t, format = "Tubes"),
  data.frame(strain = w96$strain,   hour = w96$hour,   avg = w96$avg_w,
             sd = w96$sd_w,   format = "96-well")
)

# Facet order: manually specified (left-to-right, top-to-bottom in facet_wrap).
strain_order <- c(
  "EPI300", "EPI300_ara", "EPI300_pCC1", "EPI300_pCC1_ara",
  "EPI300_pCC1_ATF1", "EPI300_pCC1_ATF1_ara",
  "EPI300_pCC1_ATF1_iAAl", "EPI300_pCC1_ATF1_ara_iAAl"
)
plot_df$strain <- factor(plot_df$strain, levels = strain_order)

p <- ggplot(plot_df, aes(x = hour, y = avg, color = format, fill = format)) +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd), width = 0.3, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.0) +
  facet_wrap(~ strain, scales = "free_y") +
  labs(
    title = "Growth curves using Avg: Tubes vs 96-well plate",
    x = "Time (h)", y = "OD (Avg \u00b1 SD)",
    color = "Format", fill = "Format"
  ) +
  scale_color_manual(values = c("Tubes" = "#1b6ca8", "96-well" = "#d9822b")) +
  scale_fill_manual(values = c("Tubes" = "#1b6ca8", "96-well" = "#d9822b")) +
  theme_bw(base_size = 11)

print(p)   # display the figure in the plot window

ggsave(file.path(output_dir, "tubes_vs_96well_Avg_curves.png"),
       plot = p, width = 11, height = 7, dpi = 150)

cat("\nSaved: tubes_vs_96well_pertimepoint_ttests.csv\n")
cat("Saved: tubes_vs_96well_strain_summary.csv\n")
cat("Saved: tubes_vs_96well_overall_paired.csv\n")
cat("Saved: tubes_vs_96well_growthcurves.png\n")

## ---- 7. Overall conclusion --------------------------------------------------
cat("\n=== Bottom line ===\n")
n_strains_sig <- sum(overall_paired$paired_p_BH < ALPHA)
cat(n_strains_sig, "out of", nrow(overall_paired),
    "strains show a statistically significant (BH-adjusted p <", ALPHA,
    ") difference between Tubes and 96-well plate across the growth curve.\n")
