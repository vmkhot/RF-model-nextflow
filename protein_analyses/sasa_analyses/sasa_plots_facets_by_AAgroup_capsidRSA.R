library(tidyverse)
library(ggpubr)

# ── Tien et al. 2013 max SASA values (Å²) ────────────────────────────────────
tien_max_sasa <- c(
  ALA = 129.0, ARG = 274.0, ASN = 195.0, ASP = 193.0, CYS = 167.0,
  GLN = 225.0, GLU = 223.0, GLY = 209.9, HIS = 224.0, ILE = 197.0,
  LEU = 201.0, LYS = 236.0, MET = 224.0, PHE = 240.0, PRO = 159.0,
  SER = 155.0, THR = 172.0, TRP = 285.0, TYR = 263.0, VAL = 174.0
)

# Load data
metadata <- read_tsv("biome_map.tsv")
sasa     <- read_tsv("sasa_capsid_estimated_all_proteins.tsv")

# Extract protein ID
sasa <- sasa |>
  mutate(protein = str_extract(protein, "^IMGVR_UViG_[^_]+_[^_]+"))

metadata <- metadata |>
  mutate(protein = str_extract(protein, "^IMGVR_UViG_[^|]+")) |>
  distinct(protein, .keep_all = TRUE)

# Three-letter to one-letter AA lookup (to match tien_max_sasa keys)
one_to_three <- c(
  A="ALA", R="ARG", N="ASN", D="ASP", C="CYS",
  Q="GLN", E="GLU", G="GLY", H="HIS", I="ILE",
  L="LEU", K="LYS", M="MET", F="PHE", P="PRO",
  S="SER", T="THR", W="TRP", Y="TYR", V="VAL"
)

# ── Amino acid groupings ─────────────────────────────────────────────────────
aa_groups <- list(
  "Positive charge" = list(aas = c("R", "H", "K"),           color = "#ef4806"),
  "Negative charge" = list(aas = c("D", "E"),                 color = "#234abb"),
  "Polar uncharged" = list(aas = c("S", "T", "N", "Q"),       color = "#1eae83"),
  "Hydrophobic"     = list(aas = c("A", "V", "I", "L", "M"),  color = "#f4b713"),
  "Aromatic"        = list(aas = c("F", "Y", "W"),            color = "#bf21a2"),
  "Cysteine"        = list(aas = c("C"),                      color = "#7a5195"),
  "Glycine"         = list(aas = c("G"),                      color = "#003f5c"),
  "Proline"         = list(aas = c("P"),                      color = "#665191")
)

# Build a lookup: AA (one-letter) -> group name
aa_to_group <- aa_groups |>
  imap(~ set_names(rep(.y, length(.x$aas)), .x$aas)) |>
  reduce(c)

group_order <- names(aa_groups)

# Merge, RSA-normalise, assign AA group, drop AAs not in Tien table (e.g. X, U)
merged <- sasa |>
  inner_join(metadata, by = "protein") |>
  mutate(
    aa_three = one_to_three[AA],                        # one-letter → three-letter
    max_sasa = tien_max_sasa[aa_three],                 # look up Tien max
    rsa      = sasa_capsid_est / max_sasa,              # RSA: 0 = buried, 1 = fully exposed
    aa_group = aa_to_group[AA]                          # map to grouping
  ) |>
  filter(!is.na(max_sasa), !is.na(aa_group)) |>          # drop non-standard AAs
  mutate(aa_group = factor(aa_group, levels = group_order))

counts <- merged |>
  group_by(ecosystem_subtype, aa_group) |>
  summarise(n = n(), .groups = "drop")

comparisons <- list(c("Lake", "Oceanic"))

# ── Strip colour per facet (matches aa_groups colours) ──────────────────────
strip_colors <- map_chr(aa_groups, "color")[group_order]

# Plot: boxplot per ecosystem, faceted by amino acid grouping
p <- ggplot(merged, aes(x = ecosystem_subtype, y = rsa, fill = ecosystem_subtype)) +
  geom_boxplot(alpha = 0.8) +
  scale_fill_manual(
    values = c(
      "Lake"    = "#9dcc7e",
      "Oceanic" = "#2B8CBF")) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "t.test",
    label       = "p.signif",
    hide.ns     = TRUE,
    size        = 7,          # controls the ** symbol size
    label.y.npc = "top",
    label.x.npc = "center",
    vjust       = 1.5,
    bracket.size = 0.8
  ) +
  geom_text(
    data = counts,
    aes(x = ecosystem_subtype, y = -0.1, label = n),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~ aa_group, ncol = 4) +
  labs(
    x = "Ecosystem",
    y = "Relative Solvent Accessibility (RSA)",
  ) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "#960089") +
  theme_classic(base_size = 15) +
  theme(
    legend.position   = "none",
    axis.text.x       = element_text(size = 15, colour = "black"),
    axis.text.y       = element_text(size = 15, colour = "black"),
    strip.text        = element_text(face = "bold", colour = "white"),
    strip.background  = element_rect(fill = "grey", color = "black", linewidth = 0.3)
  ) +
  coord_cartesian(clip = "off")   # important so text below axis is visible

p

# ── Optional: colour each facet strip by its group colour ──────────────────
# Requires ggh4x. Uncomment to use per-group strip colours instead of grey.
# library(ggh4x)
# p + facet_wrap2(~ aa_group, ncol = 4, strip = strip_themed(
#   background_x = elem_list_rect(fill = strip_colors)
# ))

# ggsave("rsa_per_aa_group_boxplot.pdf", width = 16, height = 10)