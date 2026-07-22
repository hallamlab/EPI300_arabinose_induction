library(tidyverse)

# File
csv_file <- "Arab_Induction_Tu96w_Repl_Comb.csv"

# Read everything exactly as-is
raw <- read.csv(csv_file,
                header = FALSE,
                stringsAsFactors = FALSE,
                check.names = FALSE)

#---------------------------------------------------------
# Function to convert one block to long format
#---------------------------------------------------------

extract_block <- function(df, start_row, end_row, format_name){
  
  block <- df[start_row:end_row, ]
  
  # First row contains time points
  hours <- as.numeric(unlist(block[1, -(1:2)]))
  
  # Keep only columns that have valid time values
  keep <- !is.na(hours)
  
  hours <- hours[keep]
  
  block <- block[-1, ]
  
  strain <- block[[2]]
  
  # fill down strain names
  for(i in seq_along(strain)){
    if(i > 1 && (is.na(strain[i]) || strain[i] == "")){
      strain[i] <- strain[i-1]
    }
  }
  
  values <- block[, which(keep) + 2]
  
  colnames(values) <- hours
  
  dat <- bind_cols(
    strain = strain,
    values
  )
  
  dat <- dat %>%
    filter(!is.na(strain),
           strain != "") %>%
    group_by(strain) %>%
    mutate(replicate = row_number()) %>%
    ungroup()
  
  dat %>%
    pivot_longer(
      cols = -c(strain, replicate),
      names_to = "hour",
      values_to = "OD"
    ) %>%
    mutate(
      hour = as.numeric(hour),
      OD = as.numeric(OD),
      Format = format_name
    ) %>%
    filter(!is.na(OD))
}

#---------------------------------------------------------
# Tubes block
#---------------------------------------------------------

tubes <- extract_block(
  raw,
  start_row = 2,
  end_row = 26,
  format_name = "Tubes"
)

#---------------------------------------------------------
# 96-well block
#---------------------------------------------------------

wells <- extract_block(
  raw,
  start_row = 29,
  end_row = 53,
  format_name = "96-well"
)

#---------------------------------------------------------
# Combine
#---------------------------------------------------------

growth <- bind_rows(tubes, wells)

#---------------------------------------------------------
# Summary statistics
#---------------------------------------------------------

summary_df <- growth %>%
  group_by(strain, hour, Format) %>%
  summarise(
    mean_OD = mean(OD),
    sd_OD = sd(OD),
    .groups = "drop"
  )

#---------------------------------------------------------
# Plot
#---------------------------------------------------------

p <- ggplot(summary_df,
            aes(hour,
                mean_OD,
                colour = Format,
                group = Format)) +
  
  geom_line(linewidth = 0.8) +
  
  geom_point(size = 1) +
  
  geom_errorbar(
    aes(
      ymin = mean_OD - sd_OD,
      ymax = mean_OD + sd_OD
    ),
    width = 0.3
  ) +
  
  geom_jitter(
    data = growth,
    aes(hour, OD, colour = Format),
    width = 0.15,
    alpha = 0.4,
    size = 1,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~strain, scales = "free_y") +
  
  labs(
    title = "Growth curves using triplicates: Tubes vs 96-well plate",
    x = "Time (h)",
    y = "OD600"
  ) +
  
  theme_bw() +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(p)

ggsave(
  "Growth_Curves_Tu96w_Repl_Comb.png",
  p,
  width = 12,
  height = 8,
  dpi = 300
)
