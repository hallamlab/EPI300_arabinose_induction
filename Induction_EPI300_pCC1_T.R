# ============================================================
# Plot: Induction, isoamyl alcohol (OD600 vs Time)
# Strains: EPI300_pCC1, EPI300_pCC1_ara
# ============================================================

library(tidyverse)

# ---- 1. Read the raw CSV --------------------------------------------------
# The file has a 2-row header:
#   Row 1: strain names (merged across Avg/SD column pairs)
#   Row 2: "Avg" / "SD" labels
# We read both header rows manually, then build clean column names.

raw <- read_csv("Arab_Induction.csv", col_names = FALSE, skip = 2,
                show_col_types = FALSE)

header1 <- read_csv("Arab_Induction.csv", col_names = FALSE, n_max = 1,
                    show_col_types = FALSE) |> unlist(use.names = FALSE)
header2 <- read_csv("Arab_Induction.csv", col_names = FALSE, skip = 1, n_max = 1,
                    show_col_types = FALSE) |> unlist(use.names = FALSE)

# Forward-fill the strain names across their Avg/SD column pair
strain_names <- header1
for (i in seq_along(strain_names)) {
  if (is.na(strain_names[i]) && i > 1) strain_names[i] <- strain_names[i - 1]
}

# Drop fully empty trailing columns
keep <- !is.na(header2) | seq_along(header2) == 1
raw         <- raw[, keep]
strain_names <- strain_names[keep]
header2      <- header2[keep]

# Build final column names: "Time", "<strain>_Avg", "<strain>_SD"
col_names <- c("Time",
               paste0(strain_names[-1], "_", header2[-1]))
colnames(raw) <- col_names

# ---- 2. Reshape to long format (Avg + SD, for plotting with error bars) --
# Only keep the two strains of interest
strains <- c("EPI300_pCC1", "EPI300_pCC1_ara")

df_long <- raw |>
  pivot_longer(
    cols = -Time,
    names_to = c("Strain", ".value"),
    names_pattern = "(.+)_(Avg|SD)"
  ) |>
  rename(OD600 = Avg) |>
  filter(Strain %in% strains) |>
  mutate(Strain = factor(Strain, levels = strains))

# ---- 3. Plot ---------------------------------------------------------------
# Visually distinct shapes (mix of filled/open so they're tellable apart
# even without color): square, open triangle
shape_values <- c(15, 2)

p <- ggplot(df_long, aes(x = Time, y = OD600, color = Strain, shape = Strain)) +
  geom_errorbar(aes(ymin = OD600 - SD, ymax = OD600 + SD),
                width = 0.3, linewidth = 0.4, alpha = 0.7) +
  geom_point(size = 2.4) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(0, 22, by = 2), limits = c(0, 22)) +
  scale_y_continuous(breaks = seq(0, 1.6, by = 0.2), limits = c(0, 1.6)) +
  labs(
    title = "Induction: EPI300_pCC1, Tubes",
    x = "Time (h)",
    y = "OD 600",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.key = element_blank()
  )

print(p)

ggsave("induction_EPI300_pCC1_T.png", p, width = 9, height = 5, dpi = 300)
