library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

df_ocean <- read.table("../ocean_switch_df_for_heatmap_653_NEW.tsv", header=TRUE, sep="\t")
df_lake  <- read.table("../lake_switch_df_for_heatmap_653_NEW.tsv",  header=TRUE, sep="\t")


# =========================
# Amino-acid groupings (same as before)
# =========================
aa_groups <- list(
  "Positive charged" = list(aas = c("R", "H", "K"), color = "#EF4806"),
  "Negative charged" = list(aas = c("D", "E"),       color = "#244ABB"),
  "Polar uncharged"  = list(aas = c("S", "T", "N", "Q"), color = "#1EAE83"),
  "Hydrophobic"      = list(aas = c("A", "V", "I", "L", "M"), color = "#F3B713"),
  "Aromatic"         = list(aas = c("F", "Y", "W"),  color = "#BF22A2"),
  "Other"            = list(aas = c("C", "U", "G", "P"), color = "#4B4444")
)

aa_to_group <- setNames(
  rep(names(aa_groups), sapply(aa_groups, function(g) length(g$aas))),
  unlist(lapply(aa_groups, `[[`, "aas"))
)

# =========================
# Classify each row as within- or between-group
# =========================
classify_within_between <- function(df, biome_label, exclude_self = TRUE) {
  df <- df %>%
    mutate(
      other_group = aa_to_group[other_AA],
      feat_group  = aa_to_group[feat_AA],
      comparison  = ifelse(other_group == feat_group, "Within-group", "Between-group"),
      biome = biome_label
    )
  if (exclude_self) {
    df <- df %>% filter(other_AA != feat_AA)
  }
  df
}

df_lake_class  <- classify_within_between(df_lake,  "Lake")
df_ocean_class <- classify_within_between(df_ocean, "Ocean")

df_combined <- bind_rows(df_lake_class, df_ocean_class) %>%
  mutate(
    biome      = factor(biome, levels = c("Ocean", "Lake")),
    comparison = factor(comparison, levels = c("Within-group", "Between-group"))
  )

# =========================
# Violin plot: x = biome, facet = comparison type
# Wilcoxon test: Lake vs Ocean, within each comparison type
# =========================
ggplot(df_combined, aes(x = biome, y = proportion, fill = biome)) +
  geom_violin(trim = FALSE, alpha = 1, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6, fill = "white") +
  facet_wrap(~comparison) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Ocean", "Lake")),
    label = "p.format",
    tip.length = 0.01
  ) +
  scale_fill_manual(values = c("Ocean" = "#2b8cbf", "Lake" = "#9dcc7e")) +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = "Likelihood of substitution")+
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")
  )

# =========================
# Explicit p-values / summary stats
# =========================
df_combined %>%
  group_by(comparison) %>%
  summarise(
    n_lake   = sum(biome == "Lake"),
    n_ocean  = sum(biome == "Ocean"),
    median_lake  = median(proportion[biome == "Lake"]),
    median_ocean = median(proportion[biome == "Ocean"]),
    wilcox_p = wilcox.test(
      proportion[biome == "Lake"],
      proportion[biome == "Ocean"]
    )$p.value,
    .groups = "drop"
  )
