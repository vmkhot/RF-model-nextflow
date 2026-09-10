library(dplyr)
library(ggplot2)
library(ggpubr)

# df_big <- read.csv("mapped_positions_with_biomes_2.tsv",sep='\t')

# ── Tien et al. 2013 max SASA values (Å²) for RSA normalization ──────────────
tien_max_sasa <- c(
  ALA = 129.0, ARG = 274.0, ASN = 195.0, ASP = 193.0, CYS = 167.0,
  GLN = 225.0, GLU = 223.0, GLY = 209.9, HIS = 224.0, ILE = 197.0,
  LEU = 201.0, LYS = 236.0, MET = 224.0, PHE = 240.0, PRO = 159.0,
  SER = 155.0, THR = 172.0, TRP = 285.0, TYR = 263.0, VAL = 174.0
)


# ── 1. Filter: remove rows with empty ecosystem_subtype (non reference) ──────────────────────
df <- df_big %>%
  filter(!is.na(ecosystem_subtype) & ecosystem_subtype != "")

# binning feat ranks
max_rank <- max(df$imp_order, na.rm = TRUE)
bars=20
breaks   <- c(0, seq(bars, ceiling(max_rank / bars) * bars, by = bars))
labels   <- paste0(breaks[-length(breaks)] + 1, "-", breaks[-1])

df <- df %>%
  mutate(
    max_sasa = tien_max_sasa[A_AA],          # look up max SASA by amino acid
    RSA      = surface_area / max_sasa,       # relative surface area (0–1)
    rank_bin = cut(imp_order, breaks = breaks, labels = labels,
                   include.lowest = TRUE, right = TRUE),
    burial_class = if_else(RSA < 0.2, "Buried (RSA < 20%)", "Surface-accessible"),
    exposure_class = case_when(
      RSA < 0.20 ~ "Buried",
      RSA < 0.40 ~ "Partially buried",
      RSA < 0.60 ~ "Partially exposed",
      TRUE       ~ "Exposed"),
    exposure_class = factor(exposure_class,
                            levels = c("Buried", "Partially buried",
                                       "Partially exposed", "Exposed")),
    important = if_else(imp_order <= 200, "important", "not important"))

# df <- df %>%
#   group_by(A_AA) %>%
#   mutate(
#     max_sasa_observed = max(surface_area, na.rm = TRUE),
#     RSA_observed      = surface_area / max_sasa_observed
#   ) %>%
#   ungroup()
# 
# df %>% 
#   filter(RSA > 1) %>% 
#   arrange(desc(RSA)) %>% 
#   select(A_AA, A_pos, surface_area, max_sasa, RSA) %>% 
#   head(50)
# 
# df %>%
#   filter(RSA > 1) %>%
#   summarise(min_pos = min(A_pos), max_pos = max(A_pos), 
#             n = n(), mean_pos = mean(A_pos))
# ── 1.5. Filter: keep only reference ──────────────────────
df_ref <- df_big %>% filter(protein == "EC6098_reference")
# binning feat ranks (reference)
max_rank <- max(df_ref$imp_order, na.rm = TRUE)
bars=20
breaks   <- c(0, seq(bars, ceiling(max_rank / bars) * bars, by = bars))
labels   <- paste0(breaks[-length(breaks)] + 1, "-", breaks[-1])

df_ref <- df_ref %>%
  mutate(
    max_sasa = tien_max_sasa[A_AA],          # look up max SASA by amino acid
    RSA      = surface_area / max_sasa,       # relative surface area (0–1)
    rank_bin = cut(imp_order, breaks = breaks, labels = labels,
                   include.lowest = TRUE, right = TRUE),
    burial_class = if_else(RSA < 0.2, "Buried (RSA < 20%)", "Surface-accessible"),
    exposure_class = case_when(
      RSA < 0.20 ~ "Buried",
      RSA < 0.40 ~ "Partially buried",
      RSA < 0.60 ~ "Partially exposed",
      TRUE       ~ "Exposed"
    ),
    exposure_class = factor(exposure_class,
                            levels = c("Buried", "Partially buried",
                                       "Partially exposed", "Exposed")),
   important = if_else(imp_order <= 200, "important", "not important"))

# df_ref <- df_ref %>%
#   group_by(A_AA) %>%
#   mutate(
#     max_sasa_observed = max(surface_area, na.rm = TRUE),
#     RSA_observed      = surface_area / max_sasa_observed
#   ) %>%
#   ungroup()


# ──────────────────────────────────────────────────────────────────────────────────────────
# PLOTTING
# ── scatter plot position or gini_imp v. RSA or surface area ────────────────────────────

p_scatter <- ggplot(df_ref, aes(x = feature_importance_vals, y = RSA, colour = feature_importance_vals )) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_colour_gradientn(colours = c("#CFCFEA","#BE77BF","#AD1F93"),name = "Gini\nImportance") +
  # scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Gini Importance (log scale)", y = "Relative Surface Area") +
  theme_classic(base_size = 13)

p_scatter
# ggsave("gini_v_RSA_scatter.svg")


df %>%
  filter(!is.na(RSA) & !is.na(ecosystem_subtype) & ecosystem_subtype != "") %>%
  ggplot(aes(x = ecosystem_subtype, y = RSA, fill = ecosystem_subtype)) +
  geom_violin(trim = TRUE) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  facet_wrap(~ A_AA, scales = "free_y") +
  labs(x = "Ecosystem", y = "Relative Surface Area (RSA)", title = "RSA by ecosystem per amino acid") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
# ──────────────────────────────────────────────────────────────────────────────────────────
# ── 4-class exposure + cumulative lines per class ────────────────────────────

# Proportions per bin (stacked bars — all classes, sum to 1)
df_bar <- df %>%
  filter(!is.na(rank_bin)) %>%
  count(rank_bin, exposure_class) %>%
  group_by(rank_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Cumulative % of each class's own residues across bins
df_cum <- df %>%
  filter(!is.na(rank_bin)) %>%
  count(rank_bin, exposure_class) %>%
  group_by(exposure_class) %>%
  mutate(cum_pct = cumsum(n) / sum(n)) %>%
  ungroup()

class_colors <- c(
  "Buried"            = "#2166ac",
  "Partially buried"  = "#867682",
  "Partially exposed" = "#e6ab02",
  "Exposed"           = "#b33605"
)

p_combined <- ggplot() +
  geom_col(data = df_bar,
           aes(x = rank_bin, y = prop, fill = exposure_class),
           position = "stack", color = "white", linewidth = 0.3,alpha=0.7) +
  geom_line(data = df_cum,
            aes(x = as.numeric(rank_bin), y = cum_pct,
                color = exposure_class, group = exposure_class),
            linewidth = 1.1) +
  geom_point(data = df_cum,
             aes(x = as.numeric(rank_bin), y = cum_pct,
                 color = exposure_class),
             size = 1.8) +
  scale_fill_manual(values  = class_colors) +
  scale_color_manual(values = class_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  guides(
    fill  = guide_legend(title = "Exposure class"),
    color = guide_legend(title = "Cumulative %")
  ) +
  labs(
    # title = "Residue Exposure Class per Feature Rank Bin\nwith Per-Class Cumulative Enrichment",
    x     = "Feature Rank Bin",
    y     = "Proportion / Cumulative %"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, size = 9),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13)
  )

p_combined

# ──────────────────────────────────────────────────────────────────────────────────────────
# ── Violin + boxplot: RSA by feat_rank bin ────────────────────────────────────

p_violin <- ggplot(df_ref %>% filter(!is.na(rank_bin)),
                   aes(x = rank_bin, y = RSA, fill = rank_bin)) +
  geom_violin(trim = FALSE, alpha = 0.9, color = NA) +
  geom_boxplot(width = 0.08, fill = "white", color = "grey30",
               outlier.shape = NA, alpha = 0.8) +
  geom_hline(yintercept = 0.20, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0.40, linetype = "dashed",
             color = "orange", linewidth = 0.8) +
  annotate("text", x = 0.6, y = 0.21, label = "buried",
           color = "red", size = 3, hjust = 0) +
  annotate("text", x = 0.6, y = 0.41, label = "partially buried",
           color = "orange", size = 3, hjust = 0) +
  scale_fill_viridis_d(option = "mako", guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "RSA Distribution by Feature Rank Bin",
    x     = "Feature Rank Bin",
    y     = "RSA"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    plot.title   = element_text(face = "bold", size = 13)
  )

p_violin
# ──boxplot: RSA by feat_rank bin ────────────────────────────────────

p_box <- ggplot(df %>% filter(!is.na(rank_bin)),
                aes(x = rank_bin, y = surface_area, fill = rank_bin)) +
  geom_boxplot(color = "grey30", outlier.size = 0.6,
               outlier.alpha = 0.4, alpha = 0.8) +
  # geom_hline(yintercept = 0.20, linetype = "dashed",
             # color = "red", linewidth = 0.8) +
  # geom_hline(yintercept = 0.40, linetype = "dashed",
             # color = "orange", linewidth = 0.8) +
  # annotate("text", x = 0.6, y = 0.21, label = "buried",
  #          color = "red", size = 3, hjust = 0) +
  # annotate("text", x = 0.6, y = 0.41, label = "partially buried",
  #          color = "orange", size = 3, hjust = 0) +
  scale_fill_viridis_d(option = "mako", guide = "none") +
  # scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "Feature Rank Bin",
    y     = "Surface Area (Angstrom squared)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    plot.title  = element_text(face = "bold", size = 13)
  )
p_box

