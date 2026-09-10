library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data
metadata <- read_tsv("biome_map.tsv")
sasa     <- read_tsv("../mapping_monomers_BSA/sasa_results.tsv")

# Extract protein ID from SASA filename (everything before the 3rd underscore-separated block's suffix)
sasa <- sasa |>
  mutate(protein = str_extract(protein, "^IMGVR_UViG_[^_]+_[^_]+")) |>
  group_by(protein) |>
  summarise(sasa = mean(sasa), .groups = "drop")

metadata <- metadata |>
  mutate(protein = str_extract(protein, "^IMGVR_UViG_[^|]+")) |>
  distinct(protein, .keep_all = TRUE)

# Merge and normalise
merged <- sasa |>
  inner_join(metadata, by = "protein") |>
  mutate(sasa_norm = sasa * 0.5402)

# Plot
ggplot(merged, aes(x = ecosystem_subtype, y = sasa_norm, fill = ecosystem_subtype)) +
  geom_boxplot() +
  geom_point()+
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(
    x     = "Ecosystem",
    y     = "Normalised SASA (Å² / 0.5402)",
    title = "SASA by ecosystem type"
  ) +
  theme_classic() +
  theme(legend.position = "none")
