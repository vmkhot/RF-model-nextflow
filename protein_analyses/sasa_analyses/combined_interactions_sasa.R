library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
df <- read.csv("combined_interactions_sasa.tsv",sep='\t')

df <- df %>% 
  filter(ecosystem_subtype == 0) %>%
  filter(interaction != "non-interacting")

a <- group_by(df,interaction, ecosystem_subtype) 

facet_var <- "ecosystem_subtype"

facet_maxes <- df %>%
  group_by(.data[[facet_var]]) %>%
  summarise(facet_max = max(sasa_monomer, na.rm = TRUE), .groups = "drop")

stats <- df %>% 
    group_by(ecosystem_subtype) %>%
    wilcox_test(sasa_monomer ~ interaction, paired = FALSE) %>%
    add_significance(
      p.col     = "p",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols   = c("p<0.001", "p<0.01", "p<0.05", "ns")) %>%
    # add_xy_position(data = df, x = x_var, dodge = 0.8, fun = "max") %>% # for full violins
    add_xy_position(data = df, x = "interaction", dodge = 0.25, fun = "max") %>%  # for full violins
    filter(is.finite(y.position)) %>%
    left_join(facet_maxes, by = facet_var) %>%
    mutate(y.position = pmin(y.position, facet_max * 1.1)) %>%  # cap to per-facet max + 30%
    select(-facet_max)




p <- ggplot(df,aes(x=interaction,y=sasa_monomer,fill=ecosystem_subtype)) +
  geom_violin()+
  facet_wrap(. ~ ecosystem_subtype)+
  stat_pvalue_manual(stats, label = "p", tip.length = 0.05,
                     hide.ns = TRUE, bracket.nudge.y = 0.15,size = 5)
p

p1 <- ggplot(df, aes(x=feature_importance_vals_A,y=sasa_monomer, fill=interaction)) +
  geom_point(aes(fill=interaction,colour=interaction,alpha=0.3))
p1
