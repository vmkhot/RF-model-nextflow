library(dplyr)
library(tidyr)
library(ggplot2)

# df_int <- read.csv("mapped_residue_interactions_allfeats_2.tsv", sep = "\t", na.strings = "") %>%
#   mutate(
#     biochemical_interaction = replace_na(biochemical_interaction, "non-interacting"),
#     interaction             = replace_na(interaction, "non-interacting"),
#     A_msa_col               = as.integer(A_msa_col))

aa_groups <- list(
  "Positive charged" = list(aas = c("R", "H", "K"), color = "#EF4806"),
  "Negative charged" = list(aas = c("D", "E"),       color = "#244ABB"),
  "Polar uncharged"  = list(aas = c("S", "T", "N", "Q"), color = "#1EAE83"),
  "Hydrophobic"      = list(aas = c("A", "V", "I", "L", "M"), color = "#F3B713"),
  "Aromatic"         = list(aas = c("F", "Y", "W"),  color = "#BF22A2"),
  "Other"            = list(aas = c("C", "U", "G", "P"), color = "#4B4444")
)

# Flatten into ordered AA vector + matching color vector
aa_order  <- unlist(lapply(aa_groups, function(x) x$aas), use.names = FALSE)
aa_colors <- unlist(lapply(aa_groups, function(x) rep(x$color, length(x$aas))), use.names = FALSE)
names(aa_colors) <- aa_order

# Group boundary positions (for divider lines), based on cumulative group sizes
group_sizes <- sapply(aa_groups, function(x) length(x$aas))
group_boundaries <- cumsum(group_sizes)[-length(group_sizes)] + 0.5

# Build a symmetric AA-AA pair count table for a given ecosystem_subtype,
# normalized by total residues (sum of protein lengths) in that biome
make_pair_counts <- function(df, subtype, aa_order) {
  df_sub <- df %>%
    filter(ecosystem_subtype == subtype) %>%
    filter(!is.na(AA_A), !is.na(AA_B))
  
  protein_info <- df_sub %>% distinct(protein, length)
  n_proteins  <- nrow(protein_info)
  total_res   <- sum(protein_info$length)
  
  df_counts <- df_sub %>%
    mutate(AA1 = pmin(AA_A, AA_B),
           AA2 = pmax(AA_A, AA_B)) %>%
    count(AA1, AA2, name = "n")
  
  present_levels <- intersect(aa_order, unique(c(df_counts$AA1, df_counts$AA2)))
  
  mirrored <- bind_rows(
    df_counts,
    df_counts %>% filter(AA1 != AA2) %>% rename(AA1 = AA2, AA2 = AA1)
  )
  
  mirrored %>%
    complete(AA1 = present_levels, AA2 = present_levels, fill = list(n = 0)) %>%
    mutate(
      AA1 = factor(AA1, levels = aa_order),
      AA2 = factor(AA2, levels = aa_order),
      n_norm = n / total_res,
      n_proteins = n_proteins,
      total_res = total_res
    )
}

lake_counts   <- make_pair_counts(df_int, "Lake",    aa_order)
ocean_counts  <- make_pair_counts(df_int, "Oceanic", aa_order)

plot_heatmap <- function(counts, title, aa_colors, group_boundaries) {
  n_prot <- unique(counts$n_proteins)
  tot_res <- unique(counts$total_res)
  present <- levels(droplevels(counts$AA1))
  cols_present <- aa_colors[present]
  
  ggplot(counts, aes(x = AA1, y = AA2, fill = n_norm)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = sprintf("%.5f", n_norm)), size = 2.5) +
    { if (length(group_boundaries) > 0)
      list(geom_vline(xintercept = group_boundaries, color = "black", linewidth = 0.6),
           geom_hline(yintercept = group_boundaries, color = "black", linewidth = 0.6))
    } +
    scale_fill_gradient(low = "white", high = "darkred") +
    theme_minimal(base_size = 12) +
    labs(title = title,
         subtitle = paste0("n proteins = ", n_prot, ", total residues = ", tot_res),
         x = "Amino Acid", y = "Amino Acid",
         fill = "Bonds per\nresidue") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = cols_present, face = "bold"),
      axis.text.y = element_text(colour = cols_present, face = "bold"),
      panel.grid = element_blank()
    )
}

p_lake  <- plot_heatmap(lake_counts,  "Lake proteins: normalized AA–AA bond counts", aa_colors, group_boundaries)
p_ocean <- plot_heatmap(ocean_counts, "Oceanic proteins: normalized AA–AA bond counts", aa_colors, group_boundaries)

p_lake |p_ocean
