library(ComplexHeatmap)
library(tidyverse)
library(circlize)

df_ocean <- read.table("../ocean_switch_df_for_heatmap_653_NEW.tsv", header=TRUE, sep="\t")
df_lake  <- read.table("../lake_switch_df_for_heatmap_653_NEW.tsv",  header=TRUE, sep="\t")

aa_counts_ocean <- read.csv("../ocean_aa_counts_in100_features.csv")
aa_counts_lake <- read.csv("../lake_aa_counts_in100_features.csv")

lake_col_counts <- setNames(
  aa_counts_lake$position,
  aa_counts_lake$aa
)

ocean_col_counts <- setNames(
  aa_counts_ocean$position,
  aa_counts_ocean$aa
)

# =========================
# Amino-acid groupings
# =========================
aa_groups <- list(
  "Positive charged" = list(aas = c("R", "H", "K"), color = "#EF4806"), #B22222
  "Negative charged" = list(aas = c("D", "E"),       color = "#244ABB"), #3a5fcd
  "Polar uncharged"  = list(aas = c("S", "T", "N", "Q"), color = "#1EAE83"),  #228b22
  "Hydrophobic"      = list(aas = c("A", "V", "I", "L", "M"), color = "#F3B713"),  #ffc125
  "Aromatic"         = list(aas = c("F", "Y", "W"),  color = "#BF22A2"), #cd00cd
  "Other"    = list(aas = c("C", "U", "G", "P"), color = "#4B4444")  #ff6347
)



# Ordered AA vector and colour lookup
aa_order  <- unlist(lapply(aa_groups, `[[`, "aas"))
aa_colors <- setNames(
  rep(unlist(lapply(aa_groups, `[[`, "color")), sapply(aa_groups, function(g) length(g$aas))),
  aa_order
)

# Helper: pivot + reorder rows/cols to aa_order
make_matrix <- function(df, value_col) {
  mat <- df %>%
    select(other_AA, feat_AA, all_of(value_col)) %>%
    replace_na(list(seqs=0,proportion=0,feat_counts=0)) %>%
    pivot_wider(names_from = feat_AA, values_from = all_of(value_col)) %>%
    column_to_rownames("other_AA") %>%
    as.matrix()
  
  row_order <- intersect(aa_order, rownames(mat))
  col_order <- intersect(aa_order, colnames(mat))
  mat[row_order, col_order]
}

prop_matrix_ocean <- make_matrix(df_ocean, "proportion")
prop_matrix_lake  <- make_matrix(df_lake,  "proportion")

number = "seqs"
# number = "feat_counts"
annot_matrix_ocean <- make_matrix(df_ocean, number)
annot_matrix_lake  <- make_matrix(df_lake,  number)

# =========================
# Row + column annotations
# =========================
make_annotations <- function(mat, col_counts) {
  
  row_labels <- rownames(mat)
  col_labels <- colnames(mat)
  
  row_ann <- rowAnnotation(
    Group = anno_simple(
      row_labels,
      col     = aa_colors[row_labels],
      width   = unit(6, "mm"),
      pch     = row_labels,
      pt_gp   = gpar(fontsize = 7, col = "white", fontface = "bold"),
      pt_size = unit(0, "mm"),
      gp      = gpar(alpha = 0.8)    # alpha goes here
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  # align counts to matrix columns
  col_vals <- col_counts[colnames(mat)]
  col_vals[is.na(col_vals)] <- 0
  
  col_ann <- HeatmapAnnotation(
    # 2) Barplot of counts
    Count = anno_barplot(
      col_vals,
      gp = gpar(
        fill = "#6e6e6e",
        col  = NA
      ),
      border = FALSE,
      height = unit(18, "mm")
    ),
    # 1) AA color strip
    AA = anno_simple(
      colnames(mat),
      col = aa_colors[colnames(mat)],
      height = unit(7, "mm"),
      width   = unit(6, "mm"),
      gp  = gpar(alpha = 0.8)    # alpha goes here
    ),
    
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 12),
    show_legend = FALSE
  )
  
  list(row = row_ann, col = col_ann)
}

ann_ocean <- make_annotations(prop_matrix_ocean,ocean_col_counts)
ann_lake  <- make_annotations(prop_matrix_lake, lake_col_counts)

# Shared colour scale across both matrices
col_range <- range(c(prop_matrix_ocean, prop_matrix_lake), na.rm = TRUE)
col_fun   <- colorRamp2(
  seq(col_range[1], col_range[2], length.out = 9),
  RColorBrewer::brewer.pal(9, "YlGn")
)

# =========================
# Build heatmaps
# =========================
make_heatmap <- function(prop_mat, annot_mat, ann, title) {
  Heatmap(prop_mat,
          name              = "Proportion",
          col               = col_fun,
          column_title      = title,
          column_title_gp = gpar(fontsize=16),
          cluster_rows      = FALSE,
          cluster_columns   = FALSE,
          row_names_side    = "left",
          column_names_side = "top",
          row_names_gp = gpar(fontsize = 16),
          column_names_gp = gpar(fontsize = 16),
          column_names_rot = 0,
          top_annotation    = ann$col,
          left_annotation   = ann$row,
          heatmap_legend_param = list(
            title="Proportion", fontsize=16),
          cell_fun = function(j, i, x, y, width, height, fill) {
            val <- annot_mat[i, j]
            if (!is.na(val) && val != 0) {
              grid.text(sprintf("%.0f", val), x, y,
                        gp = gpar(fontsize = 10))
            }
          }
  )
}

ht_ocean <- make_heatmap(prop_matrix_ocean, annot_matrix_ocean, ann_ocean, "Ocean-important features")
ht_lake  <- make_heatmap(prop_matrix_lake,  annot_matrix_lake,  ann_lake,  "Lake-important features")
ht_ocean
# =========================
# Legend for AA groups
# =========================
group_legend <- Legend(
  labels    = names(aa_groups),
  title_gp = gpar(fontsize=16, fontface="bold"),
  legend_gp = gpar(fill = sapply(aa_groups, `[[`, "color"),fontsize=16),
  labels_gp = gpar(fontsize = 16),
  title     = "AA Group",
  grid_width = unit(0.5, "cm"),
  grid_height = unit(0.5, "cm"),
  type = "points")

proportion_legend <- Legend(
  col_fun = "YLGN", at=c(0,0.2,0.4,0.6,0.8,1),
  legend_width = unit(0.5, "cm"),
  legend_height =  unit(0.5, "cm"),
  title = "Proportion",
  title_gp = gpar(fontsize=16, fontface="bold"),
  labels_gp = gpar(fontsize = 16))

# Draw side by side with shared legend
draw(ht_ocean + ht_lake,
     annotation_legend_list = list(group_legend,proportion_legend),
     merge_legends = TRUE)
c
