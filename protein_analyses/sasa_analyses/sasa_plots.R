library(dplyr)
library(ggplot2)
library(ggpubr)

df_big <- read.csv("mapped_positions_with_biomes.tsv",sep='\t')

# ── Tien et al. 2013 max SASA values (Å²) for RSA normalization ──────────────
tien_max_sasa <- c(
  ALA = 129.0, ARG = 274.0, ASN = 195.0, ASP = 193.0, CYS = 167.0,
  GLN = 225.0, GLU = 223.0, GLY = 104.0, HIS = 224.0, ILE = 197.0,
  LEU = 201.0, LYS = 236.0, MET = 224.0, PHE = 240.0, PRO = 159.0,
  SER = 155.0, THR = 172.0, TRP = 285.0, TYR = 263.0, VAL = 174.0
)

# ── 1. Filter: remove rows with empty ecosystem_subtype ──────────────────────
df <- df_big %>%
  filter(!is.na(ecosystem_subtype) & ecosystem_subtype != "")

# ── 2. Compute RSA and burial classification ──────────────────────────────────
df <- df %>%
  mutate(
    max_sasa = tien_max_sasa[A_AA],          # look up max SASA by amino acid
    RSA      = surface_area / max_sasa,       # relative surface area (0–1)
    buried   = RSA < 0.20                     # TRUE = buried (<20% RSA)
  )


# ── 3. Violin plot: surface_area by ecosystem_subtype ────────────────────────
p1 <- ggplot(df, aes(x = ecosystem_subtype, y = surface_area, fill = ecosystem_subtype)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.08, fill = "white", color = "grey30",
               outlier.shape = NA, alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Solvent-Accessible Surface Area by Ecosystem Subtype",
    x        = "Ecosystem Subtype",
    y        = "Surface Area (Ų)",
    fill     = "Ecosystem Subtype"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 30, hjust = 1),
    plot.title       = element_text(face = "bold", size = 14)
  )

p1
# ggsave("violin_surface_area.png", p1, width = 9, height = 5.5, dpi = 150)
# cat("Saved: violin_surface_area.png\n")

# ── 4. Bar plot: proportion buried/surface in feat_rank bins of 10 ────────────

# Create bins: 1–10, 11–20, … (cut labels like "1-10", "11-20")
max_rank <- max(df$feat_rank, na.rm = TRUE)
bars=20
breaks   <- c(0, seq(bars, ceiling(max_rank / bars) * bars, by = bars))
labels   <- paste0(breaks[-length(breaks)] + 1, "-", breaks[-1])

df <- df %>%
  mutate(
    rank_bin = cut(feat_rank, breaks = breaks, labels = labels,
                   include.lowest = TRUE, right = TRUE),
    burial_class = if_else(buried, "Buried (RSA < 20%)", "Surface-accessible"),
    important = if_else(feat_rank <= 200, "important", "not important")
  )

# Summarise proportions per bin × burial_class
df_prop <- df %>%
  filter(!is.na(rank_bin)) %>%
  count(rank_bin, burial_class) %>%
  group_by(rank_bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(df_prop, aes(x = rank_bin, y = prop, fill = burial_class)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("Buried (RSA < 20%)" = "#4E79A7",
               "Surface-accessible"  = "#F28E2B")
  ) +
  labs(
    title = "Proportion of Buried vs. Surface-Accessible Residues\nby Feature Rank Bin",
    x     = "Feature Rank Bin",
    y     = "Proportion of Residues",
    fill  = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top",
    plot.title   = element_text(face = "bold", size = 13)
  )

p2
# ggsave("barplot_burial_by_rank.png", p2, width = 10, height = 5.5, dpi = 150)
# cat("Saved: barplot_burial_by_rank.png\n")

# ── Cumulative enrichment of surface-accessible residues by feat_rank ─────────

df_bar <- df %>%
  filter(!is.na(rank_bin)) %>%
  count(rank_bin, burial_class) %>%
  group_by(rank_bin) %>%
  mutate(prop = n / sum(n)) %>%
  filter(burial_class == "Surface-accessible") %>%
  ungroup()

# cumulative % of surface residues by bin order
df_cum <- df %>%
  filter(!is.na(rank_bin), !buried) %>%
  count(rank_bin) %>%
  mutate(cum_pct = cumsum(n) / sum(n)) %>%
  ungroup()

p_combined <- ggplot() +
  geom_col(data = df_bar,
           aes(x = rank_bin, y = prop, fill = "Surface-accessible"),
           color = "white", linewidth = 0.3) +
  geom_line(data = df_cum,
            aes(x = as.numeric(rank_bin), y = cum_pct, color = "Cumulative %"),
            linewidth = 1.2) +
  geom_point(data = df_cum,
             aes(x = as.numeric(rank_bin), y = cum_pct, color = "Cumulative %"),
             size = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("Surface-accessible" = "#F28E2B")) +
  scale_color_manual(values = c("Cumulative %" = "#4E79A7")) +
  labs(
    title = "Surface-Accessible Residues per Feature Rank Bin\nwith Cumulative Enrichment",
    x     = "Feature Rank Bin",
    y     = "Proportion of Residues",
    fill  = NULL, color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = 13)
  )
p_combined


# ── 5. scatter plot: surface area by gini_imp ─────────────────────────────
df_ref <- df_big %>% filter(protein == "EC6098_reference")

df_ref <- df_ref %>%
  mutate(
    max_sasa = tien_max_sasa[A_AA],          # look up max SASA by amino acid
    RSA      = surface_area / max_sasa,       # relative surface area (0–1)
    buried   = RSA < 0.20,# TRUE = buried (<20% RSA)
    important = if_else(feat_rank <= 100, "important", "not important")
    
  )

p3 <- ggplot(df_ref, aes(x = position, y = RSA, colour = important)) +
  geom_point(aes(color = important))+
  facet_grid(~ ecosystem_subtype)
  
p3

# ── 4-class exposure + cumulative lines per class ────────────────────────────
df <- df %>%
  mutate(
    rank_bin = cut(feat_rank, breaks = breaks, labels = labels,
                   include.lowest = TRUE, right = TRUE))

df <- df %>%
  mutate(exposure_class = case_when(
    RSA < 0.20 ~ "Buried",
    RSA < 0.40 ~ "Partially buried",
    RSA < 0.60 ~ "Partially exposed",
    TRUE       ~ "Exposed"
  ),
  exposure_class = factor(exposure_class,
                          levels = c("Buried", "Partially buried",
                                     "Partially exposed", "Exposed")))

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
  "Buried"            = "#4E79A7",
  "Partially buried"  = "#76B7B2",
  "Partially exposed" = "#F28E2B",
  "Exposed"           = "#E15759"
)

p_combined <- ggplot() +
  geom_col(data = df_bar,
           aes(x = rank_bin, y = prop, fill = exposure_class),
           position = "stack", color = "white", linewidth = 0.3) +
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
    title = "Residue Exposure Class per Feature Rank Bin\nwith Per-Class Cumulative Enrichment",
    x     = "Feature Rank Bin",
    y     = "Proportion / Cumulative %"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 13)
  )

p_combined


# ── Continuous RSA: 2D hexbin density with mean RSA per bin ──────────────────

# Option 1: hexbin — shows density + mean RSA as fill
library(hexbin)
p_hex <- ggplot(df, aes(x = gini_imp, y = RSA)) +
  geom_hex(bins = 40) +
  geom_smooth(method = "gam", formula = y ~ s(x),
              color = "red", linewidth = 1, se = TRUE) +
  scale_fill_viridis_c(option = "magma", name = "Residue\ncount") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "RSA vs Feature Rank",
       x = "Feature Rank", y = "RSA") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
p_hex


# Option 2: rolling mean ± SD ribbon — smoothed RSA trend across rank
library(zoo)

df_roll <- df %>%
  arrange(feat_rank) %>%
  mutate(
    roll_mean = rollmean(RSA, k = 50, fill = NA, align = "center"),
    roll_sd   = rollapply(RSA, width = 50, FUN = sd, fill = NA, align = "center")
  )

p_roll <- ggplot(df_roll, aes(x = feat_rank)) +
  geom_point(aes(y = RSA), alpha = 0.15, size = 0.6, color = "grey60") +
  geom_ribbon(aes(ymin = roll_mean - roll_sd,
                  ymax = roll_mean + roll_sd),
              fill = "#4E79A7", alpha = 0.3) +
  geom_line(aes(y = roll_mean), color = "#4E79A7", linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "RSA vs Feature Rank (rolling mean ± 1SD, k=50)",
       x = "Feature Rank", y = "RSA") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
p_roll


# violins for important vs non-important
p_imp <- ggplot(df, aes(x=rank_bin, y=surface_area))+
  geom_violin()+
  geom_point()

  # stat_compare_means(method = "t.test", label = "p.signif")
p_imp


# ── Violin + boxplot: RSA by feat_rank bin ────────────────────────────────────

p_violin <- ggplot(df %>% filter(!is.na(rank_bin)),
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

p_box <- ggplot(df %>% filter(!is.na(rank_bin)),
                aes(x = rank_bin, y = RSA, fill = rank_bin)) +
  geom_boxplot(color = "grey30", outlier.size = 0.6,
               outlier.alpha = 0.4, alpha = 0.8) +
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title  = element_text(face = "bold", size = 13)
  )
p_box
