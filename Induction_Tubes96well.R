# ============================================================
# Plot: Induction, Tubes vs 96-well plate (OD600 vs Time, with SD error bars)
# ============================================================

library(tidyverse)
library(RColorBrewer)

file <- "Arab_Induction_Tu96w_Comb.csv"

# ---- 1. Split the file into its stacked blocks ----------------------------
# The file has two stacked tables, each with a 2-row header
# ("Tubes"/"96-well plate" + strain names, then "Avg"/"SD"),
# separated by a blank (all-comma) row.

lines <- read_lines(file)
lines[1] <- str_remove(lines[1], "^\ufeff")   # strip BOM if present

is_blank  <- str_detect(lines, "^,*$")
block_id  <- cumsum(is_blank)
blocks    <- split(lines[!is_blank], block_id[!is_blank])

# ---- 2. Parser for one block (same logic as the single-table scripts) -----
parse_block <- function(block_lines) {

  header1 <- read_csv(I(block_lines[1]), col_names = FALSE,
                       show_col_types = FALSE) |> unlist(use.names = FALSE)
  header2 <- read_csv(I(block_lines[2]), col_names = FALSE,
                       show_col_types = FALSE) |> unlist(use.names = FALSE)
  raw     <- read_csv(I(paste(block_lines[-c(1, 2)], collapse = "\n")),
                       col_names = FALSE, show_col_types = FALSE)

  method <- header1[1]   # "Tubes" or "96-well plate"

  # Forward-fill the strain names across their Avg/SD column pair
  strain_names <- header1
  for (i in seq_along(strain_names)) {
    if (is.na(strain_names[i]) && i > 1) strain_names[i] <- strain_names[i - 1]
  }

  # Drop fully empty trailing columns
  keep <- !is.na(header2) | seq_along(header2) == 1
  raw          <- raw[, keep]
  strain_names <- strain_names[keep]
  header2      <- header2[keep]

  col_names <- c("Time", paste0(strain_names[-1], "_", header2[-1]))
  colnames(raw) <- col_names

  raw |>
    pivot_longer(
      cols = -Time,
      names_to = c("Strain", ".value"),
      names_pattern = "(.*)_(Avg|SD)"
    ) |>
    rename(OD600 = Avg) |>
    mutate(
      Method = method,
      # Strip the method-specific suffix (_T or _96w) so the same strain
      # carries one label across both panels
      Strain = str_remove(Strain, "_(T|96w)$")
    )
}

df_long <- map_dfr(blocks, parse_block)

# Preserve strain order (as they first appear) and method order (Tubes, then 96-well)
strains <- unique(df_long$Strain)
methods <- unique(df_long$Method)
df_long <- df_long |>
  mutate(
    Strain = factor(Strain, levels = strains),
    Method = factor(Method, levels = methods)
  )

# ---- 3. Plot ---------------------------------------------------------------
# Shape = strain (8 distinct symbols), Color = method (Tubes vs 96-well plate)
method_colors <- c("#1f78b4", "#e31a1c")   # Tubes = blue, 96-well plate = red
shape_values  <- c(16, 15, 17, 18, 3, 4, 1, 2)

p <- ggplot(df_long, aes(x = Time, y = OD600, color = Method, shape = Strain)) +
  geom_errorbar(aes(ymin = OD600 - SD, ymax = OD600 + SD), width = 0.3,
                linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 2.0) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(0, 22, by = 2)) +
  labs(
    title = "Induction, isoamyl alcohol, tube vs 96-well plate",
    x = "Time (h)",
    y = "OD 600",
    color = NULL,
    shape = NULL
  ) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.key = element_blank()
  )

print(p)

ggsave("induction_tubes96well_combined.png", p, width = 10, height = 5.5, dpi = 300)

# ============================================================
# Plot: Induction, Tubes vs 96-well plate (OD600 vs Time, with SD error bars)
# ============================================================

library(tidyverse)
library(RColorBrewer)

file <- "Arab_Induction_Tu96w_Comb.csv"

# ---- 1. Split the file into its stacked blocks ----------------------------
# The file has two stacked tables, each with a 2-row header
# ("Tubes"/"96-well plate" + strain names, then "Avg"/"SD"),
# separated by a blank (all-comma) row.

lines <- read_lines(file)
lines[1] <- str_remove(lines[1], "^\ufeff")   # strip BOM if present

is_blank  <- str_detect(lines, "^,*$")
block_id  <- cumsum(is_blank)
blocks    <- split(lines[!is_blank], block_id[!is_blank])

# ---- 2. Parser for one block (same logic as the single-table scripts) -----
parse_block <- function(block_lines) {

  header1 <- read_csv(I(block_lines[1]), col_names = FALSE,
                       show_col_types = FALSE) |> unlist(use.names = FALSE)
  header2 <- read_csv(I(block_lines[2]), col_names = FALSE,
                       show_col_types = FALSE) |> unlist(use.names = FALSE)
  raw     <- read_csv(I(paste(block_lines[-c(1, 2)], collapse = "\n")),
                       col_names = FALSE, show_col_types = FALSE)

  method <- header1[1]   # "Tubes" or "96-well plate"

  # Forward-fill the strain names across their Avg/SD column pair
  strain_names <- header1
  for (i in seq_along(strain_names)) {
    if (is.na(strain_names[i]) && i > 1) strain_names[i] <- strain_names[i - 1]
  }

  # Drop fully empty trailing columns
  keep <- !is.na(header2) | seq_along(header2) == 1
  raw          <- raw[, keep]
  strain_names <- strain_names[keep]
  header2      <- header2[keep]

  col_names <- c("Time", paste0(strain_names[-1], "_", header2[-1]))
  colnames(raw) <- col_names

  raw |>
    pivot_longer(
      cols = -Time,
      names_to = c("Strain", ".value"),
      names_pattern = "(.*)_(Avg|SD)"
    ) |>
    rename(OD600 = Avg) |>
    mutate(
      Method = method,
      # Strip the method-specific suffix (_T or _96w) so the same strain
      # carries one label across both panels
      Strain = str_remove(Strain, "_(T|96w)$")
    )
}

df_long <- map_dfr(blocks, parse_block)

# Preserve strain order (as they first appear) and method order (Tubes, then 96-well)
strains <- unique(df_long$Strain)
methods <- unique(df_long$Method)
df_long <- df_long |>
  mutate(
    Strain = factor(Strain, levels = strains),
    Method = factor(Method, levels = methods)
  )

# ---- 3. Plot ---------------------------------------------------------------
# Shape = strain (8 distinct symbols), Color = method (Tubes vs 96-well plate)
method_colors <- c("#1f78b4", "#e31a1c")   # Tubes = blue, 96-well plate = red
shape_values  <- c(16, 15, 17, 18, 3, 4, 1, 2)

p <- ggplot(df_long, aes(x = Time, y = OD600, color = Method, shape = Strain)) +
  geom_errorbar(aes(ymin = OD600 - SD, ymax = OD600 + SD), width = 0.3,
                linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 2.0) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = seq(0, 22, by = 2)) +
  labs(
    title = "Induction, isoamyl alcohol, tube vs 96-well plate",
    x = "Time (h)",
    y = "OD 600",
    color = NULL,
    shape = NULL
  ) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.key = element_blank()
  )

print(p)

ggsave("induction_tubes96well_combined.png", p, width = 10, height = 5.5, dpi = 300)

