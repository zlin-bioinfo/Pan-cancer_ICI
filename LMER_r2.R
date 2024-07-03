rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','tibble','qs','MetBrewer','forcats','lmerTest','grid','rstatix','ggpubr','RColorBrewer','ggrepel')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Change at r2 level
meta_int <- read.csv('./tables/meta_int.csv') 
freq_mat <- meta_int |> 
  distinct(celltype_r2, sample, .keep_all = T) |> 
  select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  pivot_longer(cols = c('Pre', 'Post'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(celltype_r2) |> 
  mutate(freq_r2_comp_scale = scale(freq_r2_comp)) |> 
  ungroup()
# Linear Mixed-Effects Regression Models
uni_lmer <- function(freq_mat){
  subtypes <- unique(freq_mat$celltype_r2)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    freq_mat$timepoint <- ifelse(freq_mat$time_point == 'Pre', 0, 1)
    freq_mat$int_cat <- ifelse(freq_mat$int_cat == '<21d', 0, 1)
    freq_mat <- freq_mat |> filter(celltype_r2 == subtype)
    variables_to_use <- c('freq_r2_comp_scale', 'dataset'[length(unique(freq_mat$dataset)) > 1], 'modality'[length(freq_mat$modality)], 'int_cat'[length(freq_mat$int_cat)])
    variables_to_use <- c(variables_to_use, "(1 | patient)")
    variables_to_use <- variables_to_use[!is.na(variables_to_use)]
    formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
    # formula <- as.formula("timepoint ~ freq_r2_comp_scale + dataset + modality + int_cat + (1 | patient)")
    model <- lmer(formula, data = freq_mat, REML = FALSE)
    model_summary <- summary(model)
    coef_table <- coef(summary(model))
    confint_table <- confint(model, level = 0.95)
    subtype_results <- coef_table['freq_r2_comp_scale', , drop = FALSE]
    subtype_confint <- confint_table['freq_r2_comp_scale', , drop = FALSE]
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = subtype_results[1],
      StdError = subtype_results[2],
      tValue = subtype_results[4],
      pValue = coef_table['freq_r2_comp_scale', 'Pr(>|t|)'],
      CI_lower = subtype_confint[1],
      CI_upper = subtype_confint[2]
    )
    return(combined_results)
  })
  df <- do.call(rbind, uni_models) |> data.frame()
  df$fdr <- p.adjust(df$pValue, method = 'fdr', n = nrow(df)) 
  # celltype_keep <- filter(df, pValue < 0.05) |> pull(Celltypes) |> unique()
  # df <- filter(df, Celltypes %in% celltype_keep) |> 
  #   arrange(desc(Estimate)) 
  df <- df |> 
    arrange(desc(Estimate)) 
  df$fdr_cat <- '<0.05'
  df$fdr_cat[df$fdr>0.05] <- '>0.05'
  return(df)
}

# Overall
res_lmer <- uni_lmer(freq_mat)
write.csv(res_lmer, file = './tables/lmer_all.csv', row.names = F)
res_lmer <- filter(res_lmer, pValue < 0.05)

# Making plot
pdf('./figures/Change/uni_lmer_r2.pdf', height = 6, width = 5)
p <- ggplot(res_lmer, aes(x= factor(Celltypes, levels = rev(res_lmer$Celltypes)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(aes(color = fdr_cat), size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.001), -log10(0.005), -log10(0.01)), labels = c('<0.001', '<0.005', '<0.01'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-0.3, 0.3)) +
  xlab("") +
  scale_color_manual(values = c("black", "grey")) +
  coord_flip() + theme_minimal() +
  labs(size = "P value", color = "FDR") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 10, colour = "black"))
grid.draw(p)
y_text = 0.12
x1 = 0.3
x2 = 0.75
grid.segments(
  x = unit(x1 + 0.05, "npc"),
  x1 = unit(x1 - 0.05, "npc"),
  y = unit(y_text+0.02, "npc"),
  y1 = unit(y_text+0.02, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.05, "npc"),
  x1 = unit(x2 + 0.05, "npc"),
  y = unit(y_text+0.02, "npc"),
  y1 = unit(y_text+0.02, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("On", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
dev.off()
