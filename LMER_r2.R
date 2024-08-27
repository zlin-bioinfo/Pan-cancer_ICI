rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','tibble','qs','MetBrewer','forcats','lmerTest','grid','rstatix','ggpubr','RColorBrewer','ggrepel','effsize')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Change at r2 level
meta_int <- read.csv('./PanCancer_ICI/tables/meta_int.csv') 
freq_mat <- meta_int |> 
  distinct(celltype_r2, sample, .keep_all = T) |> 
  select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
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
write.csv(res_lmer, file = './PanCancer_ICI/tables/lmer_all.csv', row.names = F)
res_lmer <- filter(res_lmer, pValue < 0.05)

# Making plot
pdf('./PanCancer_ICI/figures/Change/uni_lmer_r2.pdf', height = 6, width = 5)
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

order <- c("Overall", "SKCM_Becker", "BCC_Yost", "SCC_Yost","BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "HNSC_Franken", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "CRC_Chen", "NSCLC_Liu")
subtypes <- res_lmer |> 
  filter(fdr<0.05) |> 
  select(Celltypes) |> 
  pull()
list_res <- lapply(subtypes, function(subtype){
  print(subtype)
  res_lmer |> filter(Celltypes == subtype)
  freq_mat_sub <- meta_int |> 
    distinct(celltype_r2, sample, .keep_all = T) |> 
    select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
    pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
    pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
    filter(celltype_r2 == subtype) |> 
    group_by(dataset) |> 
    mutate(freq_r2_comp_scale = scale(freq_r2_comp),
           n_sample = n()/2) |> 
    filter(n_sample > 2) |> 
    ungroup()
  uni_ds <- lapply(unique(freq_mat_sub$dataset), function(ds){
    mat <- freq_mat_sub[freq_mat_sub$dataset == ds,]
    # print(ds)
    mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
    # variables_to_use <- c('freq_r2_comp_scale', 'response'[length(unique(mat$response)) > 1])
    # variables_to_use <- c(variables_to_use, "(1 | patient)")
    # formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
    formula <- as.formula("timepoint~freq_r2_comp_scale + (1 | patient)")
    model <- lmer(formula, data = mat, REML = FALSE)
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
      CI_upper = subtype_confint[2],
      dataset = ds
    )
    return(combined_results)
  })
  res_df <- do.call(rbind, uni_ds)
  res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
  res_subtype_all$fdr <- 'Overall'
  names(res_subtype_all)[[8]] <- 'dataset'
  res_df <- rbind(res_subtype_all, res_df) 
})
names(list_res) <- subtypes
i <- 6
subtype <- subtypes[[i]]; subtype
res <- list_res[[i]]
min(res$CI_lower)
max(res$CI_upper)
pdf(paste0('./PanCancer_ICI/figures/Change/', subtype, '.pdf'), height = 4, width = 5)
p <- ggplot(res, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-0.6, 0.6)) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.2, "npc"),
  y1 = unit(0.2, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.2, "npc"),
  y1 = unit(0.2, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.3, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 8))
grid.text("On", x = unit(0.75, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 8))
dev.off()

# boxplot
df_add <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% res_lmer$Celltypes[!grepl('ISG15', res_lmer$Celltypes)]) |> 
  mutate(celltype_r2 = factor(celltype_r2, levels = res_lmer$Celltypes)) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(direction = ifelse(On > Pre*1.1, 'Up', ifelse(On < Pre*0.9, 'Down', 'not'))) |>
  group_by(celltype_r2) |> 
  mutate(pt_count = n()) |> 
  group_by(direction, .add = T) |> 
  mutate(dir_count = n()) |> 
  distinct(celltype_r2, direction, pt_count, dir_count) |> 
  filter(!direction == 'not') |> 
  pivot_wider(values_from = 'dir_count', names_from = 'direction')
# P1
stat.test <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% res_lmer$Celltypes[!grepl('ISG15', res_lmer$Celltypes)]) |> 
  mutate(celltype_r2 = factor(celltype_r2, levels = res_lmer$Celltypes)) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(celltype_r2) |> 
  t_test(freq_r2_comp~time_point, paired = T) |> 
  add_significance("p")
stat.test <- merge(stat.test, df_add[df_add$celltype_r2 %in% res_lmer$Celltypes ,], by = 'celltype_r2') 
stat.test <- stat.test |> 
  mutate(celltype_r2 = factor(celltype_r2, levels = res_lmer$Celltypes)) |> 
  arrange(celltype_r2) |> 
  mutate(strip_label = paste0('Up: ', Up, '/', pt_count, '\n Down: ', Down, '/', pt_count, '\n p=', p)) 

meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% res_lmer$Celltypes[!grepl('ISG15|B_LMO2', res_lmer$Celltypes)]) |> 
  mutate(celltype_r2 = factor(celltype_r2, levels = res_lmer$Celltypes)) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  merge(stat.test, by = 'celltype_r2') |> 
  ggplot(aes(x = factor(time_point, levels = c('Pre', 'On')), y = freq_r2_comp)) +
  geom_violin(aes(fill = time_point), alpha = 0.6) +
  geom_boxplot(width = 0.3, alpha = 0.2) +
  geom_line(aes(group = patient), linetype = "dashed", size = 0.2, alpha = 0.3) +
  geom_point(aes(color = time_point), size = 0.1) +
  scale_color_manual(values = c('#154999', '#CF0034')) +
  scale_fill_manual(values = c('#154999', '#CF0034')) +
  facet_wrap(.~celltype_r2 + strip_label, scales = 'free', ncol = 5) +
  xlab("") + ylab("Relative Frequency") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 8, colour = 'black'),
        strip.text = element_text(size = 8, margin = margin(2,1,2,1))) +
  guides(fill="none", color = "none",
         color = guide_legend(override.aes = list(size=6))) 
ggsave('./PanCancer_ICI/figures/Change/boxplot.pdf', width = 8, height = 7)

# ISG15 specific
subtypes <- unique(meta_int$celltype_r2[str_detect(meta_int$celltype_r2,'ISG15')])
subtypes <- factor(subtypes, levels = c('CD4_Th_ISG15', 'CD4_Treg_ISG15', 'CD8_Tex_ISG15', 'cDC2_ISG15', 'Macro_ISG15', 'B_ISG15'))
subtypes <- sort(subtypes)
df_add <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% subtypes) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(direction = ifelse(On > Pre*1.1, 'Up', ifelse(On < Pre*0.9, 'Down', 'not'))) |>
  group_by(celltype_r2) |> 
  mutate(pt_count = n()) |> 
  group_by(direction, .add = T) |> 
  mutate(dir_count = n()) |> 
  distinct(celltype_r2, direction, pt_count, dir_count) |> 
  filter(!direction == 'not') |> 
  pivot_wider(values_from = 'dir_count', names_from = 'direction')
stat.test <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% subtypes) |>  
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(celltype_r2) |> 
  t_test(freq_r2_comp~time_point, paired = T) |> 
  add_significance("p")
stat.test <- merge(stat.test, df_add[df_add$celltype_r2 %in% subtypes ,], by = 'celltype_r2') 
stat.test <- stat.test |> 
  mutate(strip_label = paste0('Up: ', Up, '/', pt_count, '\n Down: ', Down, '/', pt_count, '\n p=', p)) 

meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% subtypes) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  merge(stat.test, by = 'celltype_r2') |> 
  ggplot(aes(x = factor(time_point, levels = c('Pre', 'On')), y = freq_r2_comp)) +
  geom_violin(aes(fill = time_point), alpha = 0.6) +
  geom_boxplot(width = 0.3, alpha = 0.2) +
  geom_line(aes(group = patient), linetype = "dashed", size = 0.2, alpha = 0.3) +
  geom_point(aes(color = time_point), size = 0.1) +
  scale_color_manual(values = c('#154999', '#CF0034')) +
  scale_fill_manual(values = c('#154999', '#CF0034')) +
  facet_wrap(.~factor(celltype_r2, levels = subtypes) + strip_label, scales = 'free', ncol = 6) +
  xlab("") + ylab("Relative Frequency") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 8, colour = 'black'),
        strip.text = element_text(size = 8, margin = margin(2,1,2,1))) +
  guides(fill="none", color = "none",
         color = guide_legend(override.aes = list(size=6))) 
ggsave('./PanCancer_ICI/figures/Change/boxplot_ISG.pdf', width = 10, height = 3)

freq_mat <- meta_int |> 
  distinct(celltype_r2, sample, .keep_all = T) |> 
  select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(celltype_r2) |> 
  mutate(freq_r2_comp_scale = scale(freq_r2_comp)) |> 
  ungroup()
uni_lmer <- function(freq_mat){
  subtypes <- unique(freq_mat$celltype_r2)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    freq_mat$timepoint <- ifelse(freq_mat$time_point == 'Pre', 0, 1)
    freq_mat$int_cat <- ifelse(freq_mat$int_cat == '<21d', 0, 1)
    formula <- as.formula("timepoint ~ freq_r2_comp_scale + dataset + int_cat + (1 | patient)")
    freq_mat <- freq_mat |> filter(celltype_r2 == subtype)
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
  df$group <- unique(freq_mat$modality)
  df$fdr_cat <- '<0.05'
  df$fdr_cat[df$fdr>0.05] <- '>0.05'
  # celltype_keep <- filter(df, pValue < 0.05) |> pull(Celltypes) |> unique()
  # df <- filter(df, Celltypes %in% celltype_keep)
  return(df)
}

# Modality
list_res <- split(freq_mat, freq_mat$modality) |> lapply(uni_lmer) 
pdf('./PanCancer_ICI/figures/Change/scatter_modality.pdf', height = 5, width = 6.5)
com_celltypes <- intersect(list_res[[1]]$Celltypes, list_res[[2]]$Celltypes)
list_res[[1]] <- list_res[[1]][match(com_celltypes, list_res[[1]]$Celltypes),]
list_res[[2]] <- list_res[[2]][match(com_celltypes, list_res[[2]]$Celltypes),]
dual <- list_res[[1]] |> select(Celltypes, Estimate, CI_lower, CI_upper, pValue) |> rename_with(~ paste0(., "_", names(list_res)[1]), .cols = c("Estimate", "CI_lower", "CI_upper", "pValue"))
mono <- list_res[[2]] |> select(Celltypes, Estimate, CI_lower, CI_upper, pValue) |> rename_with(~ paste0(., "_", names(list_res)[2]), .cols = c("Estimate", "CI_lower", "CI_upper", "pValue"))
df <- merge(dual, mono)
# write.csv(df, file = './PanCancer_ICI/tables/lmer_modality.csv', row.names = F)
df |> mutate(celltype_label = ifelse(abs(Estimate_Mono - Estimate_Dual) > 0.2, Celltypes,
                                     ifelse((pValue_Mono < 0.05 & pValue_Dual>= 0.05)|(pValue_Mono >= 0.05 & pValue_Dual < 0.05), Celltypes,
                                            ifelse(pValue_Mono < 0.05 & pValue_Dual<0.05, Celltypes, ''))),
             Significance = ifelse(pValue_Mono < 0.05 & pValue_Dual<0.05, '<0.05 in Both', 
                                   ifelse(pValue_Mono < 0.05 & pValue_Dual>= 0.05, '<0.05 in Mono', 
                                          ifelse(pValue_Mono >= 0.05 & pValue_Dual < 0.05, '<0.05 in Dual', 'Non-signif')))) |>
  mutate(Significance = factor(Significance, levels = c('<0.05 in Mono', '<0.05 in Dual', '<0.05 in Both', 'Non-signif'))) |> 
  ggplot(aes(x = Estimate_Mono, y = Estimate_Dual)) + 
  geom_point(aes(color = Significance), size = 2) + 
  scale_color_manual(values = c('#E41A1C', '#377EB8', '#4DAF4A', '#999999'), name = 'P value') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
  theme_classic() + 
  xlim(c(-0.4, 0.4)) +
  ylim(c(-0.8, 0.8)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "lightgray", size = 0.05, linetype = "dashed"),
        panel.grid.minor = element_line(color = "lightgray", size = 0.05, linetype = "dashed")) +
  geom_text_repel(aes(label = celltype_label), box.padding = 0.7, size = 2.5, max.overlaps = 50) + 
  xlab("Coefficient\n(Mono)") + 
  ylab("Coefficient\n (Dual)") + ggtitle('Modality') + annotate("text", x = -0.35, y = -0.45, label = "y = x", colour = '#999999')
y_text = 0.033
x1 = 0.15
x2 = 0.7
grid.segments(
  x = unit(x1 + 0.05, "npc"),
  x1 = unit(x1 - 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.05, "npc"),
  x1 = unit(x2 + 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("On", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))

x_text = 0.05
y1 = 0.17
y2 = 0.87
grid.segments(
  x = unit(x_text + 0.02, "npc"),
  x1 = unit(x_text + 0.02, "npc"),
  y = unit(y1 + 0.05, "npc"),
  y1 = unit(y1 - 0.05, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x_text + 0.02, "npc"),
  x1 = unit(x_text + 0.02, "npc"),
  y = unit(y2 - 0.05, "npc"),
  y1 = unit(y2 + 0.05, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", y = unit(y1, "npc"), x = unit(x_text, "npc"), gp = gpar(fontsize = 8), rot = 90)
grid.text("On", y = unit(y2, "npc"), x = unit(x_text, "npc"), gp = gpar(fontsize = 8), rot = 90)
dev.off()

res <- do.call(rbind, list_res) |> data.frame()
subtypes <- res |> 
  filter(pValue<0.05) |> 
  pull(Celltypes) |> 
  unique()
res <- res |> filter(Celltypes %in% subtypes) |> arrange(desc(Estimate)) |> mutate(group = factor(group, levels = c('Mono', 'Dual')))
order_row <- res |> 
  select(Celltypes, Estimate, group) |> 
  pivot_wider(values_from = 'Estimate', names_from = 'group') |> 
  mutate(diff = Dual - Mono) |> 
  arrange(desc(diff), desc(Dual)) |> 
  pull(Celltypes)
res$CI_lower[res$CI_lower < -1] <- -1
res$CI_upper[res$CI_upper > 1] <- 1
dotCOLS = c("#FDBF6F","#CAB2D6")
barCOLS = c("#FF7F00","#6A3D9A")
pdf('./PanCancer_ICI/figures/Change/uni_lmer_modality.pdf', height = 5, width = 5)
ggplot(res, aes(x= factor(Celltypes, levels = rev(order_row)), y=Estimate, ymin=CI_lower, ymax=CI_upper,col=group, fill=group, size = -log10(fdr))) +
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05)), labels = c('<0.001', '<0.01', '<0.05'), range = c(0,4)) +
  scale_fill_manual(values=barCOLS) +
  scale_color_manual(values=dotCOLS) +
  scale_y_continuous(name= "Effect Size", limits = c(-1, 1)) + xlab("") +
  coord_flip() +
  theme_minimal() +
  labs(fill = "Modality", color = "Modality", size = "FDR") +
  guides(size = guide_legend(override.aes = list(fill = "black"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
y_text = 0.14
x1 = 0.3
x2 = 0.7
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



pt_df <- read.csv('./PanCancer_ICI/tables/meta_patient.csv')
pt_df$tx_ds <- paste0(pt_df$treatment, ' (', pt_df$dataset, ')')
ct_order <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
              'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
              'CD4_Treg_Early', 'CD4_Treg_ISG15', 'CD4_Treg_TNFRSF9', 
              'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
              'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
              'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13',
              'CD8_Tex_ISG15', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 
              'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo',
              'B_Naive', 'B_ISG15', 'B_HSP', 'B_MT2A', 'ACB_EGR1', 'ACB_NR4A2', 'ACB_CCR7', 'B_Memory', 'B_AtM',
              'GCB_Pre', 'GCB_SUGCT', 'GCB_LMO2', 'GCB_Prolif', 'Plasmablast', 'Plasma_cell',
              'Mast','pDC','cDC1', 
              'cDC2_CD1C', 'cDC2_IL1B','cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'MigrDC', 'MoDC', 
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro_ISG15', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2',
              'Endo_lymphatic','Endo_artery','Endo_capillary','Endo_tip','Endo_vein',
              'EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
subtypes <- df |> 
  filter(pValue_Dual < 0.05 | pValue_Mono < 0.05) |> 
  pull(Celltypes) |> 
  factor(levels = ct_order) |> 
  sort() |> 
  as.character()
meta_int$dataset[meta_int$dataset %in% c('BCC_Yost', 'SCC_Yost')] <- 'BCC/SCC_Yost'
meta_int$tx_ds <- paste0(meta_int$treatment, ' (', meta_int$dataset, ')')
check_and_add_columns <- function(data, columns) {
  for (column in columns) {
    if (!column %in% names(data)) {
      data[[column]] <- 0
    }
  }
  return(data)
}
cohorts <- c("aPD1 (SKCM_Becker)", "aPD1 (BCC/SCC_Yost)", 
             "aPD1 (BRCA_Bassez1)", "aPD1(pre-Chemo) (BRCA_Bassez2)", "aPD1 (TNBC_Shiao)",
             "aPD1 (HNSC_Luoma)", "aPD1 (CRC_Li)", "aPD1 (CRC_Chen)", "aPD1+Chemo (CRC_Chen)", "aPD1+Chemo (NSCLC_Liu)", 
             "aPDL1+Chemo (TNBC_Zhang)", 
             "aPDL1 (HNSC_Franken)", "aPD1+CTLA4 (HNSC_IMCISION)", "aPDL1+CTLA4 (HNSC_Franken)")
res_list <- lapply(cohorts, function(cohort){
  print(cohort)
  df <- filter(meta_int, tx_ds == cohort) |> 
    select(freq_r2_comp, time_point, patient, celltype_r2) |> 
    distinct(celltype_r2, patient, time_point, .keep_all = T)
  dep_comp <- lapply(subtypes, function(subtype){
    print(subtype)
    df_sub <- df |> 
      filter(celltype_r2 == subtype) |>
      pivot_wider(names_from = time_point, values_from = freq_r2_comp, values_fill = 0)
    if (nrow(df_sub) >= 3){
      # Check for columns 'Pre' and 'On'
      df_sub <- check_and_add_columns(df_sub, c("Pre", "On"))
      t <- t.test(df_sub$On, df_sub$Pre, paired = T, alternative = 'two.sided')
      cohen_d <- cohen.d(df_sub$On, df_sub$Pre, paired = TRUE, hedges.correction = TRUE)
      comp <- c(subtype, t$statistic, t$p.value, cohen_d$estimate)
      return(comp)
    } else {
      return(c(subtype, NA, NA, NA))
    }
  })
  results <- do.call(rbind, dep_comp) |> data.frame()
  rownames(results) <- subtypes
  colnames(results) <- c('celltype','t_score', 'pvalue','effectsize')
  results$fdr <- p.adjust(results$pvalue, method = 'fdr', n = nrow(results))
  return(results)
})
names(res_list) <- cohorts
mat_es <- res_list |> lapply(function(res){return(as.numeric(res$effectsize))}) 
mat_es <- do.call(rbind, mat_es)
colnames(mat_es) <- subtypes
mat_sig <- res_list |> lapply(function(res){return(res$pvalue)}) 
mat_sig <- do.call(rbind, mat_sig)
colnames(mat_sig) <- subtypes
pdf('./PanCancer_ICI/figures/Change/ht_tx_ds.pdf', height = 5, width = 10)
Heatmap(mat_es, name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), na_col = 'lightgray',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "right", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8"), legend_direction = "horizontal", legend_side = 'bottom'),
        column_names_gp = gpar(fontsize = 10),
        row_split = factor(c(rep('Mono',12),rep('Dual',2)), levels = c('Mono','Dual')),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6.3, "mm"), 
        height = nrow(mat_es)*unit(6.3, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()


list_res <- split(freq_mat, freq_mat$response) |> lapply(uni_lmer) 
names(list_res) <- unique(freq_mat$response)
com_celltypes <- intersect(list_res[[1]]$Celltypes, list_res[[2]]$Celltypes)
list_res[[1]] <- list_res[[1]][match(com_celltypes, list_res[[1]]$Celltypes),]
list_res[[2]] <- list_res[[2]][match(com_celltypes, list_res[[2]]$Celltypes),]
nr <- list_res[[1]] |> select(Celltypes, Estimate, CI_lower, CI_upper, pValue) |> rename_with(~ paste0(., "_", names(list_res)[1]), .cols = c("Estimate", "CI_lower", "CI_upper", "pValue"))
re <- list_res[[2]] |> select(Celltypes, Estimate, CI_lower, CI_upper, pValue) |> rename_with(~ paste0(., "_", names(list_res)[2]), .cols = c("Estimate", "CI_lower", "CI_upper", "pValue"))
df <- merge(nr, re)
write.csv(df, file = './PanCancer_ICI/tables/lmer_response.csv', row.names = F)
pdf('./PanCancer_ICI/figures/Change/scatter_response.pdf', height = 5, width = 6.5)
df |> mutate(celltype_label = ifelse(abs(Estimate_RE - Estimate_NR) > 0.2, Celltypes,
                                     ifelse((pValue_RE < 0.05 & pValue_NR>= 0.05)|(pValue_RE >= 0.05 & pValue_NR < 0.05), Celltypes,
                                            ifelse(pValue_RE < 0.05 & pValue_NR<0.05, Celltypes, ''))),
             Significance = ifelse(pValue_RE < 0.05 & pValue_NR<0.05, '<0.05 in Both', 
                                   ifelse(pValue_RE < 0.05 & pValue_NR>= 0.05, '<0.05 in RE', 
                                          ifelse(pValue_RE >= 0.05 & pValue_NR < 0.05, '<0.05 in NR', 'Non-signif')))) |>
  mutate(Significance = factor(Significance, levels = c('<0.05 in RE', '<0.05 in NR', '<0.05 in Both', 'Non-signif'))) |> 
  ggplot(aes(x = Estimate_RE, y = Estimate_NR)) + 
  geom_point(aes(color = Significance), size = 2) + 
  scale_color_manual(values = c('#E41A1C', '#377EB8', '#4DAF4A', '#999999'), name = 'P value') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.1) +
  theme_classic() + 
  xlim(c(-0.3, 0.5)) +
  ylim(c(-0.3, 0.3)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "lightgray", size = 0.05, linetype = "dashed"),
        panel.grid.minor = element_line(color = "lightgray", size = 0.05, linetype = "dashed")) +
  geom_text_repel(aes(label = celltype_label), box.padding = 0.7, size = 2.5, max.overlaps = 50) + 
  xlab("Coefficient\n(RE)") +
  ylab("Coefficient\n (NR)") + ggtitle('Response') + annotate("text", x = 0.3, y = 0.25, label = "y = x", colour = '#999999')
y_text = 0.033
x1 = 0.15
x2 = 0.7
grid.segments(
  x = unit(x1 + 0.05, "npc"),
  x1 = unit(x1 - 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.05, "npc"),
  x1 = unit(x2 + 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Down", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("Up", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))

x_text = 0.05
y1 = 0.17
y2 = 0.87
grid.segments(
  x = unit(x_text + 0.02, "npc"),
  x1 = unit(x_text + 0.02, "npc"),
  y = unit(y1 + 0.05, "npc"),
  y1 = unit(y1 - 0.05, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x_text + 0.02, "npc"),
  x1 = unit(x_text + 0.02, "npc"),
  y = unit(y2 - 0.05, "npc"),
  y1 = unit(y2 + 0.05, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Down", y = unit(y1, "npc"), x = unit(x_text, "npc"), gp = gpar(fontsize = 8), rot = 90)
grid.text("Up", y = unit(y2, "npc"), x = unit(x_text, "npc"), gp = gpar(fontsize = 8), rot = 90)
dev.off()

pt_df <- read.csv('./PanCancer_ICI/tables/meta_patient.csv')
pt_df$re_ds <- paste0(pt_df$response, ' (', pt_df$dataset, ')')
ct_order <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
              'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
              'CD4_Treg_Early', 'CD4_Treg_ISG15', 'CD4_Treg_TNFRSF9', 
              'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
              'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
              'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13',
              'CD8_Tex_ISG15', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 
              'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo',
              'B_Naive', 'B_ISG15', 'B_HSP', 'B_MT2A', 'ACB_EGR1', 'ACB_NR4A2', 'ACB_CCR7', 'B_Memory', 'B_AtM',
              'GCB_Pre', 'GCB_SUGCT', 'GCB_LMO2', 'GCB_Prolif', 'Plasmablast', 'Plasma_cell',
              'Mast','pDC','cDC1', 
              'cDC2_CD1C', 'cDC2_IL1B','cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'MigrDC', 'MoDC', 
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro_ISG15', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2',
              'Endo_lymphatic','Endo_artery','Endo_capillary','Endo_tip','Endo_vein',
              'EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
subtypes <- df |> 
  filter(pValue_RE < 0.05 | pValue_NR < 0.05) |> 
  pull(Celltypes) |> 
  factor(levels = ct_order) |> 
  sort() |> 
  as.character()
meta_int$dataset[meta_int$dataset %in% c('BCC_Yost','SCC_Yost')] <- 'BCC/SCC_Yost'
meta_int$re_ds <- paste0(meta_int$response, ' (', meta_int$dataset, ')')
check_and_add_columns <- function(data, columns) {
  for (column in columns) {
    if (!column %in% names(data)) {
      data[[column]] <- 0
    }
  }
  return(data)
}
cohorts <- pt_df |> janitor::tabyl(re_ds) |> filter(n>2) |> pull(re_ds)
# cohorts <- c("NR (SKCM_Becker)", "NR (BCC/SCC_Yost)",
#              "NR (BRCA_Bassez1)", "NR (BRCA_Bassez2)","NR (TNBC_Shiao)", "NR (TNBC_Zhang)",
#              "NR (HNSC_Franken)", "NR (HNSC_IMCISION)", "NR (HNSC_Luoma)",
#              "RE (BCC/SCC_Yost)", "RE (BRCA_Bassez1)", "RE (BRCA_Bassez2)","RE (TNBC_Shiao)",
#              "RE (HNSC_Franken)", "RE (HNSC_IMCISION)", "RE (CRC_Li)", "RE (NSCLC_Liu)")
res_list <- lapply(cohorts, function(cohort){
  print(cohort)
  df <- filter(meta_int, re_ds == cohort) |> 
    select(freq_r2_comp, time_point, patient, celltype_r2) |> 
    distinct(celltype_r2, patient, time_point, .keep_all = T)
  dep_comp <- lapply(subtypes, function(subtype){
    print(subtype)
    df_sub <- df |> 
      filter(celltype_r2 == subtype) |>
      pivot_wider(names_from = time_point, values_from = freq_r2_comp, values_fill = 0)
    if (nrow(df_sub) >= 3){
      # Check for columns 'Pre' and 'On'
      df_sub <- check_and_add_columns(df_sub, c("Pre", "On"))
      t <- t.test(df_sub$On, df_sub$Pre, paired = T, alternative = 'two.sided')
      cohen_d <- cohen.d(df_sub$On, df_sub$Pre, paired = TRUE, hedges.correction = TRUE)
      comp <- c(subtype, t$statistic, t$p.value, cohen_d$estimate)
      return(comp)
    } else {
      return(c(subtype, NA, NA, NA))
    }
  })
  results <- do.call(rbind, dep_comp) |> data.frame()
  rownames(results) <- subtypes
  colnames(results) <- c('celltype','t_score', 'pvalue','effectsize')
  results$fdr <- p.adjust(results$pvalue, method = 'fdr', n = nrow(results))
  return(results)
})
names(res_list) <- cohorts
mat_es <- res_list |> lapply(function(res){return(as.numeric(res$effectsize))}) 
mat_es <- do.call(rbind, mat_es)
colnames(mat_es) <- subtypes
mat_sig <- res_list |> lapply(function(res){return(res$pvalue)}) 
mat_sig <- do.call(rbind, mat_sig)
colnames(mat_sig) <- subtypes
pdf('./PanCancer_ICI/figures/Change/ht_re_ds.pdf', height = 6, width = 10)
Heatmap(mat_es, name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), na_col = 'lightgray',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "right", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8"), legend_direction = "horizontal", legend_side = 'bottom'),
        column_names_gp = gpar(fontsize = 10),
        row_split = factor(c(rep('Non-responder',10),rep('Responder',9)), levels = c('Non-responder','Responder')),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6.5, "mm"), 
        height = nrow(mat_es)*unit(6.5, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()

# Response
uni_lmer <- function(freq_mat){
  subtypes <- unique(freq_mat$celltype_r2)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    freq_mat$timepoint <- ifelse(freq_mat$time_point == 'Pre', 0, 1)
    freq_mat$int_cat <- ifelse(freq_mat$int_cat == '<21d', 0, 1)
    formula <- as.formula("timepoint ~ freq_r2_comp_scale + dataset + int_cat + (1 | patient)")
    freq_mat <- freq_mat |> filter(celltype_r2 == subtype)
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
  df$group <- unique(freq_mat$response)
  df$fdr_cat <- '<0.05'
  df$fdr_cat[df$fdr>0.05] <- '>0.05'
  # celltype_keep <- filter(df, pValue < 0.05) |> pull(Celltypes) |> unique()
  # df <- filter(df, Celltypes %in% celltype_keep)
  return(df)
}
list_res <- split(freq_mat, freq_mat$response) |> lapply(uni_lmer) 
res <- do.call(rbind, list_res) |> data.frame()
subtypes <- res |> 
  filter(pValue<0.05) |> 
  pull(Celltypes) |> 
  unique()
res <- res |> filter(Celltypes %in% subtypes) |> arrange(desc(Estimate)) |> mutate(group = factor(group, levels = c('NR','RE')))
order_row <- res |> 
  select(Celltypes, Estimate, group) |> 
  pivot_wider(values_from = 'Estimate', names_from = 'group') |> 
  mutate(diff = RE - NR) |> 
  arrange(desc(diff), desc(NR)) |> 
  pull(Celltypes)
res$CI_upper[res$CI_upper >0.6] <- 0.6
dotCOLS = c("#A6CEE3","#FB9A99")
barCOLS = c("#1F78B4","#E31A1C")
pdf('./PanCancer_ICI/figures/Change/uni_lmer_response.pdf', height = 6, width = 5)
ggplot(res, aes(x= factor(Celltypes, levels = rev(order_row)), y=Estimate, ymin=CI_lower, ymax=CI_upper,col=group, fill=group, size = -log10(fdr))) +
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05)), labels = c('<0.001', '<0.01', '<0.05'), range = c(0,4)) +
  scale_fill_manual(values=barCOLS) +
  scale_color_manual(values=dotCOLS) +
  scale_y_continuous(name= "Effect Size", limits = c(-0.4, 0.6)) + xlab("") +
  coord_flip() +
  theme_minimal() +
  labs(fill = "Response", color = "Response", size = "FDR") +
  guides(size = guide_legend(override.aes = list(fill = "black"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
y_text = 0.11
x1 = 0.3
x2 = 0.7
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



