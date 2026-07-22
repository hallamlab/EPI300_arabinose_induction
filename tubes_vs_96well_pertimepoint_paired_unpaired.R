## =============================================================================
## Per-timepoint Tubes vs 96-well t-tests (PAIRED and UNPAIRED), from raw
## triplicate OD data.
## =============================================================================
## For every strain x hour combination present in BOTH formats, this runs:
##   - an UNPAIRED (Welch) two-sample t-test: treats the 3 tube replicates and
##     3 well replicates as independent samples. Always valid.
##   - a PAIRED t-test: pairs tube replicate i with well replicate i (i = 1,2,3).
##     ONLY valid if replicate i really is the "same" biological unit measured
##     in both formats (e.g. same starter culture split into a tube and a
##     well on the same day). If your replicates are just independent
##     cultures with no such 1-to-1 correspondence, ignore the paired columns
##     and use the unpaired test.
##
## Requires: raw triplicate CSV already parsed the same way as in
## tubes_vs_96well_analysis_triplicates.R (parsing code repeated here so this
## script can be run standalone).
## =============================================================================

## ---- 0. Setup --------------------------------------------------------------
ALPHA      <- 0.05
input_file <- "Arab_Induction_Tu96w_Repl_Comb.csv"
output_dir <- "."

## ---- 1. Parse the raw triplicate CSV (same logic as before) ----------------
raw <- readLines(input_file, encoding = "UTF-8", warn = FALSE)
raw <- sub("^\uFEFF", "", raw)
raw <- gsub("\r$", "", raw)

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
    length(cells) <- max(length(cells), n_hours + 2)

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
wells_long$strain <- sub("_ATF1__ara_iAAl$", "_ATF1_ara_iAAl", wells_long$strain)

raw_long <- rbind(tubes_long, wells_long)
raw_long <- raw_long[!is.na(raw_long$OD), ]

## ---- 2. Per-timepoint paired + unpaired t-tests -----------------------------
common_hours <- intersect(unique(tubes_long$hour), unique(wells_long$hour))
strains <- unique(raw_long$strain)

pertimepoint <- do.call(rbind, lapply(strains, function(s) {
  do.call(rbind, lapply(common_hours, function(h) {

    t_df <- raw_long[raw_long$strain == s & raw_long$Format == "Tubes"   & raw_long$hour == h, ]
    w_df <- raw_long[raw_long$strain == s & raw_long$Format == "96-well" & raw_long$hour == h, ]
    t_df <- t_df[order(t_df$replicate), ]
    w_df <- w_df[order(w_df$replicate), ]

    if (nrow(t_df) < 2 || nrow(w_df) < 2) return(NULL)

    ## --- Unpaired (Welch) two-sample t-test: always computed ---
    tt_unpaired <- t.test(t_df$OD, w_df$OD)  # var.equal = FALSE (Welch) by default

    ## --- Paired t-test: only if replicate counts match, pairing rep i <-> rep i ---
    paired_ok <- nrow(t_df) == nrow(w_df) && all(t_df$replicate == w_df$replicate)
    if (paired_ok) {
      tt_paired <- t.test(t_df$OD, w_df$OD, paired = TRUE)
      paired_t  <- unname(tt_paired$statistic)
      paired_df <- unname(tt_paired$parameter)
      paired_p  <- tt_paired$p.value
    } else {
      paired_t <- paired_df <- paired_p <- NA
    }

    data.frame(
      strain = s, hour = h,
      avg_t = mean(t_df$OD), sd_t = sd(t_df$OD),
      avg_w = mean(w_df$OD), sd_w = sd(w_df$OD),
      unpaired_t = unname(tt_unpaired$statistic),
      unpaired_df = unname(tt_unpaired$parameter),
      unpaired_p = tt_unpaired$p.value,
      paired_t = paired_t,
      paired_df = paired_df,
      paired_p = paired_p
    )
  }))
}))

## BH correction, done separately within each test type (and across all
## strain x timepoint rows, matching the earlier script's convention)
pertimepoint$unpaired_p_BH <- p.adjust(pertimepoint$unpaired_p, method = "BH")
pertimepoint$paired_p_BH   <- p.adjust(pertimepoint$paired_p,   method = "BH")
pertimepoint <- pertimepoint[order(pertimepoint$strain, pertimepoint$hour), ]

cat("=== Per-timepoint Tubes vs 96-well t-tests, from raw triplicates ===\n")
cat("(paired columns assume tube replicate i corresponds to well replicate i;\n")
cat(" set to NA if that assumption doesn't hold for your experimental design)\n\n")
print(pertimepoint, digits = 3, row.names = FALSE)

write.csv(pertimepoint,
          file.path(output_dir, "tubes_vs_96well_pertimepoint_paired_unpaired.csv"),
          row.names = FALSE)

cat("\nSaved: tubes_vs_96well_pertimepoint_paired_unpaired.csv\n")

## ---- 3. Quick summary: do paired and unpaired tests agree? -----------------
agree_summary <- data.frame(
  n_timepoints = nrow(pertimepoint),
  n_sig_unpaired_p05 = sum(pertimepoint$unpaired_p < ALPHA, na.rm = TRUE),
  n_sig_paired_p05   = sum(pertimepoint$paired_p   < ALPHA, na.rm = TRUE),
  n_sig_unpaired_BH  = sum(pertimepoint$unpaired_p_BH < ALPHA, na.rm = TRUE),
  n_sig_paired_BH    = sum(pertimepoint$paired_p_BH   < ALPHA, na.rm = TRUE)
)
cat("\n=== Agreement summary across all strain x timepoint rows ===\n")
print(agree_summary, row.names = FALSE)
