library(ggbiplot)
library(GGally)
library(tidyverse)

# Load data
rna <- read_tsv("path_to>rnaseq_all_CDS_RPKM_no_filter_orf_scaled.tab")
rp1 <- read_tsv("path_to>rp_all_CDS_RPKM_no_filter_cult1.tab")
rp2 <- read_tsv("path_to>rp_all_CDS_RPKM_no_filter_cult2.tab")
ms <- read_csv("path_to>ms_df_long_annotated.csv")
protein_properties <- read_csv("path_to>protein_properties_syn6803_20180315.csv")

################################################################################################
# RNA data workup
################################################################################################
# Scale RPKM to 1 M
rna <- rna %>% group_by(Sample) %>% mutate(RPKM = RPKM *10^6 / sum(RPKM)) 

# Filter out observations with >30 mapped reads
rna <- filter(rna, Reads > 30)

# Add variables for time point and replicate
rna <- rna %>% separate(col = Sample, into = c("Sample", "Replicate"), sep = 1) %>%
  mutate(Time_point = recode(
    Sample, "A" = "BRISE", "B" = "ARISE", "C" = "MDAY", "D" = "BSET", "E" = "ASET")) %>%
  select(-Sample, -Replicate, Replicate ) %>% 
  mutate(Time_point = factor(Time_point, levels = c("BRISE", "ARISE", "MDAY", "BSET", "ASET")))

# Remove unnecessary variables
rna <- select(rna, Name, Time_point:Replicate, RPKM)

# log2 transform RPKM
rna <- rna %>% mutate(log2rna = log2(RPKM))

# Remove log2RPKM < 0
rna <- rna %>% filter(log2rna >= 0) 

# Remove genes that have <2 replicates at any timepoint
tp <- c("BRISE", "ARISE", "MDAY", "BSET", "ASET")
rna <- rna %>% group_by(Name) %>% filter(all(tp %in% Time_point))
rna <- rna %>% group_by(Name, Time_point) %>% mutate(n_repl = length(Replicate))
rna <- rna %>% group_by(Name) %>% filter(all(n_repl > 1)) %>% select(-n_repl) %>% ungroup
rna %>% group_by(Replicate, Time_point) %>% summarise(n = n())  

# Calculate mean log2RPKM and standard deviation
rna <- rna %>% group_by(Time_point, Name) %>% mutate(mean_log2rna = mean(log2rna),
                                                     sd_log2rna = sd(log2rna)) %>%
  ungroup

# Make rna data tibble wider for quality check analysis
rna_w <- rna %>% select(-RPKM, -mean_log2rna:-sd_log2rna) %>% arrange(Time_point) %>%
  pivot_wider(names_from = Time_point:Replicate, values_from = log2rna) %>% ungroup

# Quality check plots ##########################################################################
# # Correlation plots
# ggpairs(select(rna_w, -Name),
#         lower = list(continuous = wrap("points", alpha = 0.1, size=0.5),
#                      combo = wrap("dot", alpha = 0.4, size=2)),
#         upper = list(continuous = wrap("cor", size = 3)))

# pca
pca_matrix <- rna_w %>% column_to_rownames("Name") %>% t
pca <- prcomp(pca_matrix, center = T, scale = F)
plot(pca)
summary(pca)
ggbiplot(pca, choices = c(1,2), var.axes = F, labels = rownames(pca_matrix))

# Plot SD vs mean abundance
ggplot(rna, aes(x = mean_log2rna, y = sd_log2rna)) + 
  geom_point(alpha = 0.3) + facet_grid(Time_point~.) +
  scale_x_continuous(breaks = seq(0,100,1)) +
  theme_bw()

# Check distribution of standard deviations
ggplot(rna, aes(group = as.character(Replicate), y = sd_log2rna)) +
  geom_boxplot() + facet_grid(.~Time_point) +
  scale_y_continuous(breaks = seq(0,100,0.2)) + theme_bw()

# Check distribution of RPKM
ggplot(rna, aes(group = as.character(Replicate), y = log2rna)) +
  geom_boxplot() + facet_grid(.~Time_point) + theme_bw()

# Check normality
ggplot(rna, aes(sample = log2rna)) + 
  stat_qq() + stat_qq_line() + facet_grid(Time_point~Replicate) + theme_bw()

################################################################################################
# RP data workup
################################################################################################
# Rename RP samples in cultivation 2 to A3, B3, C3
rp2 <- rp2 %>% mutate(Sample = recode(Sample, "A1" = "A3", "B1" = "B3", "E1" = "E3"))

# Combine data from cultivation 1 and 2 (CCE1 and CCE2)
rp <- bind_rows(rp1, rp2)
rm(rp1, rp2)

# Scale RPKM to 1 M
rp <- rp %>% group_by(Sample) %>% mutate(RPKM = RPKM * 10^6 / sum(RPKM)) 

# Filter out observations with >60 mapped reads
rp <- filter(rp, Reads > 60)

# Remove A2 sample because coverage is low 67% (others ~95%)
rp <- filter(rp, Sample != "A2")

# Add variables for time point and replicate
rp <- rp %>% separate(col = Sample, into = c("Sample", "Replicate"), sep = 1) %>%
  mutate(Time_point = recode(
    Sample, "A" = "BRISE", "B" = "ARISE", "C" = "MDAY", "D" = "BSET", "E" = "ASET")) %>%
  mutate(Time_point = factor(Time_point, levels = c("BRISE", "ARISE", "MDAY", "BSET", "ASET"))) %>%
  select(-Sample)

# Remove unnecessary variables
rp <- select(rp, Name, Time_point, Replicate, RPKM)

# log2 transform RPKM
rp <- rp %>% mutate(log2rib = log2(RPKM))

# Remove log2RPKM < 0
rp <- rp %>% filter(log2rib >= 0) 

# Remove genes that have <2 replicates at any timepoint
rp %>% group_by(Replicate, Time_point) %>% summarise(n = n())
tp <- c("BRISE", "ARISE", "MDAY", "BSET", "ASET")
rp <- rp %>% group_by(Name) %>% filter(all(tp %in% Time_point))
rp <- rp %>% group_by(Name, Time_point) %>% mutate(n_repl = length(Replicate))
rp <- rp %>% group_by(Name) %>% filter(all(n_repl > 1)) %>% select(-n_repl) %>% ungroup
rp %>% group_by(Replicate, Time_point) %>% summarise(n = n())

# Calculate mean log2RPKM and standard deviation
rp <- rp %>% group_by(Time_point, Name) %>% mutate(mean_log2rib = mean(log2rib),
                                                     sd_log2rib = sd(log2rib)) %>%
  ungroup

# Make rna data wider for quality check analysis
rp_w <- rp %>% select(-RPKM, -mean_log2rib:-sd_log2rib) %>% arrange(Time_point) %>%
  pivot_wider(names_from = Time_point:Replicate, values_from = log2rib) %>% ungroup
rp_w <- rp_w %>% filter(!is.na(rp_w$ARISE_2))

# Quality check plots ##########################################################################
# Correlation plots (Supplemental figure 1)
# c <- ggpairs(select(rp_w, -Name),
#              lower = list(continuous = wrap("points", alpha = 0.1, size=0.5),
#                           combo = wrap("dot", alpha = 0.4, size=2)),
#              upper = list(continuous = wrap("cor", size = 3)))
# ggsave("path",
#      c, height= 25/2.54, width=25/2.54)

# pca
pca_matrix <- rp_w %>% column_to_rownames("Name") %>% t
pca <- prcomp(pca_matrix, center = T, scale = F)
plot(pca)
summary(pca)
ggbiplot(pca, choices = c(1,2), var.axes = F, labels = rownames(pca_matrix))

# Plot SD vs mean abundance
ggplot(rp, aes(x = mean_log2rib, y = sd_log2rib)) + 
  geom_point(alpha = 0.3) + facet_grid(Time_point~.) +
  scale_x_continuous(breaks = seq(0,100,1)) +
  theme_bw()

# Check distribution of standard deviations
ggplot(rp, aes(group = as.character(Replicate), y = sd_log2rib)) +
  geom_boxplot() + facet_grid(.~Time_point) +
  scale_y_continuous(breaks = seq(0,100,0.2)) + theme_bw()

# Check distribution of RPKM
ggplot(rp, aes(group = as.character(Replicate), y = log2rib)) +
  geom_boxplot() + facet_grid(.~Time_point) + theme_bw()

# Check normality
ggplot(rp, aes(sample = log2rib)) + 
  stat_qq() + stat_qq_line() + facet_grid(Time_point~Replicate) + theme_bw()

################################################################################################
# MS data workup
################################################################################################
# Structure data table
ms <- ms %>% select(protein:R2, -n_peptides, -time) %>%
  gather(key = "replicate", value = "intensity", R1, R2) %>% 
  mutate(replicate = recode(replicate, "R1" = 1, "R2" = 2)) %>%
  rename(Name = protein, Time_point = condition) %>%
  mutate(Time_point = factor(Time_point, levels = c("BRISE", "ARISE", "MDAY", "BSET", "ASET")))

# Remove genes that have <2 replicates at any timepoint
ms <- ms %>% group_by(Name) %>% filter(!any(is.na(intensity))) %>% ungroup

# Switch name on arise1 and mday2 (mislabelled)
ms <- ms %>% unite("run", sep = "_", Time_point:replicate, remove = T)
ms <- ms %>% mutate(run = recode(run,
  "MDAY_2" = "ARISE_1",
  "ARISE_1" = "MDAY_2"
))
ms <- ms %>% separate(run, c("Time_point", "replicate"), remove = T)
ms <- ms %>% mutate(Time_point = factor(Time_point, levels = c("BRISE", "ARISE", "MDAY", "BSET", "ASET")))

# Scale total intensity to 10^12
ms %>% group_by(Time_point, replicate) %>% summarise(sum = sum(intensity, na.rm = T))
ms <- ms %>% group_by(Time_point, replicate) %>%
  mutate(intensity = intensity * 10^12 / sum(intensity, na.rm = T)) %>% ungroup

# log2 transform intensities
ms <- ms %>% group_by(Time_point, replicate) %>% mutate(log2int = log2(intensity))

# Calculate mean log2 intensity standard deviation
ms <- ms %>% group_by(Time_point, Name) %>% mutate(mean_log2int = mean(log2int),
                                                   sd_log2int = sd(log2int)) %>%
  ungroup

# Remove genes that has at least one time point with sd(log2int) < 1.
ms %>% group_by(Time_point, replicate) %>% summarise(n = n())
ms <- ms %>% filter(sd_log2int < 1)
ms %>% group_by(Time_point, replicate) %>% summarise(n = n())

tp <- c("BRISE", "ARISE", "MDAY", "BSET", "ASET")
ms <- ms %>% group_by(Name) %>% filter(all(tp %in% Time_point))
ms <- ms %>% group_by(Name, Time_point) %>% mutate(n_repl = length(replicate))
ms <- ms %>% group_by(Name) %>% filter(all(n_repl > 1)) %>% select(-n_repl)
ms %>% group_by(Time_point, replicate) %>% summarise(n = n())

# Make ms data wider for quality check analysis
ms_w <- ms %>% select(-intensity, -mean_log2int:-sd_log2int) %>% arrange(Time_point) %>%
  pivot_wider(names_from = Time_point:replicate, values_from = log2int) %>% ungroup

### QC Plots ###################################################################################
# # Correlation plots
# ggpairs(select(ms_w, -Name),
#         lower = list(continuous = wrap("points", alpha = 0.1, size=0.5),
#                      combo = wrap("dot", alpha = 0.4, size=2)),
#         upper = list(continuous = wrap("cor", size = 3)))

# pca
pca_matrix <- ms_w %>% column_to_rownames("Name") %>% t
pca <- prcomp(pca_matrix, center = T, scale = F)
plot(pca)
summary(pca)
ggbiplot(pca, choices = c(1,2), var.axes = F, labels = rownames(pca_matrix))

# Plot SD vs mean abundance
ggplot(ms, aes(x = mean_log2int, y = sd_log2int)) + 
  geom_point(alpha = 0.3) + facet_grid(Time_point~.) +
  geom_vline(xintercept = 22.5, color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0,100,1)) +
  theme_bw()

# Check distribution of standard deviations
ggplot(ms, aes(group = as.character(replicate), y = sd_log2int)) +
  geom_boxplot() + facet_grid(.~Time_point) +
  scale_y_continuous(breaks = seq(0,100,0.2)) + theme_bw()

# Check distribution of intensities
ggplot(ms, aes(group = as.character(replicate), y = log2int)) +
  geom_boxplot() + facet_grid(.~Time_point) + theme_bw()

# Check normality
ggplot(ms, aes(sample = log2int)) + 
  stat_qq() + stat_qq_line() + facet_grid(Time_point~replicate) + theme_bw()

shapiro_test_w <- function(data) {
  result <- shapiro.test(data)
  return(result[[1]])
}
shapiro_test_p <- function(data) {
  result <- shapiro.test(data)
  return(result[[2]])
}
ms %>% group_by(Time_point, replicate) %>%
  summarise(W = shapiro_test_w(log2int), p = shapiro_test_p(log2int))

################################################################################################
# Identify shared genes in all three omics data sets (RNAseq, RP and MS) after work up
################################################################################################
shared_genes <- inner_join(distinct(select(rna, Name)), distinct(select(rp, Name)))
shared_genes <- inner_join(shared_genes, distinct(select(ms, Name)))

################################################################################################
# Clustering analysis of mRNA expression
################################################################################################
# Filter genes that are shared accross all omics data sets
rna <- rna %>% filter(Name %in% shared_genes$Name)

# Test if at least one time point has different abundance
# Welch's anova
anova <- function(x, y) {
  d = tibble(Time_point = x, expr = y)
  stat = oneway.test(expr~Time_point, d, var.equal = F)
  return(as.numeric(stat[[3]]))   # return p-value
}
rna <- rna %>% group_by(Name) %>% mutate(p_val = anova(Time_point, log2rna))

# Distribution of p-values
ggplot(rna, aes(x = p_val)) +
  geom_histogram(bins= 40) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw()

# Adjust p-values (Benjamini-Hochberg)
p_adj <- rna %>% select(Name, p_val) %>% ungroup() %>% distinct() %>%
  mutate(q = p.adjust(p_val, method = "fdr"))
rna <- select(rna, -p_val) %>% inner_join(p_adj)
rm(p_adj)

# Center log2rna values and calulate mean values at each time point
rna <- rna %>% group_by(Name) %>% mutate(clog2rna = log2rna - mean(log2rna)) %>% ungroup
rna <- rna %>% group_by(Name, Time_point) %>% mutate(mean_clog2rna = mean(clog2rna)) %>% ungroup

# Reduce data table to only the average value at each time point
rna <- rna %>% select(Name:Time_point, mean_log2rna, mean_clog2rna, sd_log2rna:q) %>% distinct

# Filter genes which have > 2 fold difference in RPKM between at least 2 different time points
rna_diff <- rna %>% group_by(Name) %>% filter( max(mean_log2rna) - min(mean_log2rna) >= 1 )
  
# Calculate fraction (%) >2-fold genes (of shared genes)
perc_changing_genes <- 100 * nrow(rna_diff) / nrow(rna)

# Filter sign. changing genes (controlled FDR = 0.10, q < 0.1)
rna_diff <- rna_diff %>% filter(q < 0.10)
rna_nodiff <- rna %>% filter( !(Name %in% rna_diff$Name) )

# Calculate fraction (%) sign. changing genes (of shared genes)
perc_sgn_changing_genes <- 100 * nrow(rna_diff) / nrow(rna)

# Cluster genes according to diurnal mRNA pattern
rna_diff_w <- rna_diff %>% select(Name, Time_point, mean_clog2rna) %>%
  pivot_wider(names_from = Time_point, values_from = mean_clog2rna) %>%
  column_to_rownames("Name")

library(amap)
rna_clusters <- hcluster(rna_diff_w, method = "pearson", link = "ward", nbproc = 3)

# Determine average silhouette width
dist <- Dist(rna_diff_w, method = "pearson")
library(cluster)
av_sil_width <- sapply(2:40, function(x) {
  summary(
    silhouette(as.numeric(cutree(rna_clusters, k=x)), dist = dist)
  )$si.summary[4]
})
plot(x = 2:40, y = av_sil_width, xlab = "Number of clusters",
     ylab = "Average silhouette width", ylim = c(0,max(av_sil_width)))

# Cut dendogram at different cluster number cutoffs.
clusters_rna <- tibble("Name" = names(cutree(rna_clusters, k = 2)),
                       "c2" = cutree(rna_clusters, k = 2),
                       "c3" = cutree(rna_clusters, k = 3),
                       "c4" = cutree(rna_clusters, k = 4),
                       "c5" = cutree(rna_clusters, k = 5))

# Assign non-changing genes to "cluster 0"
clusters_rna_nodiff <- tibble("Name" = rna_nodiff$Name,
                              "c2" = 0, "c3" = 0, "c4" = 0, "c5" = 0) %>% distinct()
clusters_rna <- bind_rows(clusters_rna, clusters_rna_nodiff)

# Merge cluster annotations with expression levels
rna_clust <- inner_join(rna, clusters_rna)

# Convert time point labels to actual time point 
rna_clust <- rna_clust %>% mutate(Time = recode(
  Time_point, "BRISE" = 23, "ARISE" = 1, "MDAY" = 6, "BSET" = 11, "ASET" = 13))

# Add a copy of the ARISE timepoint (Time = 25) for plotting
ARISE_copy <- rna_clust %>% filter(Time_point == "ARISE") %>% mutate(Time = 25)
rna_clust <- bind_rows(rna_clust, ARISE_copy)
rm(ARISE_copy)

################################################################################################
# Assign mRNA cluster annotations to RP data set
################################################################################################
# Reduce RP data set to shared genes across all omics data sets
rp <- rp %>% filter(Name %in% shared_genes$Name)

# Center log2rib values and calulate mean values at each time point
rp <- rp %>% group_by(Name) %>% mutate(clog2rib = log2rib - mean(unique(mean_log2rib))) %>% ungroup
rp <- rp %>% group_by(Name, Time_point) %>% mutate(mean_clog2rib = mean(clog2rib)) %>% ungroup

# Reduce data table to only the average value at each time point
rp <- rp %>% select(Name:Time_point, mean_log2rib, mean_clog2rib, sd_log2rib) %>% distinct()

# Merge mRNA cluster annotation with RP expression levels
rp_clust <- inner_join(rp, clusters_rna)

# Convert time points labels to actual time point 
rp_clust <- rp_clust %>% mutate(Time = recode(
  Time_point, "BRISE" = 23, "ARISE" = 1, "MDAY" = 6, "BSET" = 11, "ASET" = 13))

# Add a copy of the ARISE timepoint (Time = 25) for plotting
ARISE_copy <- rp_clust %>% filter(Time_point == "ARISE") %>% mutate(Time = 25)
rp_clust <- bind_rows(rp_clust, ARISE_copy)
rm(ARISE_copy)

################################################################################################
# Assess differential protein abundance
################################################################################################
# Reduce MS data set to shared genes across all omics data sets
ms <- ms %>% filter(Name %in% shared_genes$Name)

# Test if at least one time point has different abundance
# Welch's anova
anova <- function(x, y) {
  d = tibble(Time_point = x, expr = y)
  stat = oneway.test(expr~Time_point, d, var.equal = F)
  return(as.numeric(stat[[3]]))   # return p-value
}
ms <- ms %>% group_by(Name) %>% mutate(p_val = anova(Time_point, log2int))

# Distribution of p-values
ggplot(ms, aes(x = p_val)) +
  geom_histogram(bins= 40) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw()

# Adjust p-values (Benjamini-Hochberg)
p_adj <- ms %>% select(Name, p_val) %>% ungroup() %>% distinct() %>%
  mutate(q = p.adjust(p_val, method = "fdr"))
ms <- select(ms, -p_val) %>% inner_join(p_adj)
rm(p_adj)

# Center log2int values and calulate mean values at each time point
ms <- ms %>% group_by(Name) %>% mutate(clog2int = log2int - mean(log2int)) %>% ungroup
ms <- ms %>% group_by(Name, Time_point) %>% mutate(mean_clog2int = mean(clog2int)) %>% ungroup

# Calculate total sum of squares (sst), between group SS (ssb) and within group SS (error, sse)
# Only diff mRNA genes. FOR INDIVIDUAL GENES.
ms0 <- ms %>% filter(Name %in% rna_diff$Name)
ms_ssx <- ms0 %>% group_by(Name, Time_point) %>% mutate(error = clog2int-mean_clog2int) %>%
  group_by(Name) %>% mutate(sst = sum(clog2int^2),
                            ssb = sum(mean_clog2int^2),
                            sse = sum(error^2)) %>%
  mutate(sum_eb = ssb + sse, .before = sst) %>% mutate(ssb_sse = ssb/sse) 
ms_ssx_red <- ms_ssx %>% select(Name, sst:ssb_sse) %>% distinct()
ggplot(ms_ssx_red, aes(y = log2(ssb_sse))) + geom_boxplot() +
  scale_y_continuous(breaks = seq(-6, 7, 1))
perc_error_genes <- nrow(filter(ms_ssx_red, ssb_sse < 1)) / nrow(ms_ssx_red)

# Reduce data table to only the average value at each time point
ms <- ms %>% select(Name:Time_point, mean_log2int, mean_clog2int, sd_log2int:q) %>%
  group_by(Name, Time_point) %>% distinct()

# Make supplementary table 1. Differential protein abundance analysis result.
ms_anot <- left_join(ms, protein_properties, by = c("Name" = "GeneID")) %>% group_by(Name) %>%
  mutate(max_log2fc = max(mean_log2int) - min(mean_log2int), .after = sd_log2int) %>% ungroup %>%
  select(-Time_point:-sd_log2int, p_val:q) %>% distinct() %>%
  relocate(Protein, Gene.names, Pathway, Process, Catalytic.activity,
           psortB_localization, psortB_score, .after = q)
# write_csv(ms_anot, "path")

################################################################################################
# Assign mRNA cluster annotations to MS data set
################################################################################################
# Merge mRNA cluster annotation with MS expression levels
ms_clust <- inner_join(ms, clusters_rna)

# Convert time points labels to actual time point 
ms_clust <- ms_clust %>% mutate(Time = recode(
  Time_point, "BRISE" = 23, "ARISE" = 1, "MDAY" = 6, "BSET" = 11, "ASET" = 13))

# Add a copy of the ARISE timepoint (Time = 25) for plotting
ARISE_copy <- ms_clust %>% filter(Time_point == "ARISE") %>% mutate(Time = 25)
ms_clust <- bind_rows(ms_clust, ARISE_copy)
rm(ARISE_copy)

################################################################################################
# Analyze oscillation amplitudes
################################################################################################
# Join RNA seq, RP and MS data (only cyclic genes)
amp <- inner_join(select(rna_diff, Name:Time_point, mean_clog2rna, mean_log2rna),
                  select(rp, Name:Time_point, mean_clog2rib, mean_log2rib),
                  by = c("Name", "Time_point"))
amp <- inner_join(amp,
                  select(ms, Name:Time_point, mean_clog2int, mean_log2int))

# Calc. the relative mRNA, RP and protein oscillation amplitude for each gene (maximum log2FC)
amp <- amp %>% group_by(Name) %>% mutate(mRNA_amp = max(mean_clog2rna) - min(mean_clog2rna),
                                         rib_amp = max(mean_clog2rib) - min(mean_clog2rib),
                                         protein_amp = max(mean_clog2int) - min(mean_clog2int))

# Calculate daily average mRNA, ribosomes, protein (rescale protein first to make comparable)
amp_red <- amp %>% group_by(Name) %>% mutate(day_avg_rna = mean(mean_log2rna),
                                             day_avg_rib = mean(mean_log2rib),
                                             day_avg_prot = mean(log2(2^mean_log2int/10^6)))

# Reduce to only gene-specific values
amp_red <- amp_red %>% select(-Time_point:-mean_log2int) %>% distinct
summary(amp_red)
# Plot protein amplitude vs. daily synthesis to daily protein abundace ratio in plot section...
# Plot median relative amplitudes in plot section...

################################################################################################
# Plot RNA expression in different clusters 
################################################################################################
p <- ggplot(rna_clust, aes(x = Time, y = mean_clog2rna, group = Name))
p <- p + annotate(
  "rect", xmin=12, xmax=24,
  ymin=-3.3,
  ymax=3.3,
  fill="black", alpha=0.1)

p <- p + geom_line(alpha = 0.5, size = 0.1, color = "#216b21ff") 
p <- p + scale_x_continuous(breaks = seq(0, 24, by = 6), minor_breaks = F,
                            expand = expansion(mult = c(0.05, 0.02)))
p <- p + scale_y_continuous(limits = c(-3.3, 3.3),
                            breaks = seq(-3, 3, 1), minor_breaks = F,
                            expand = expansion(mult = c(0)))
p <- p + theme_bw()
p <- p + theme(legend.position = "none",
               line = element_line(size = 0.2),
               panel.border = element_rect(size = 0.2, fill = NA),
               strip.background = element_blank(),
               strip.text.y = element_blank(),
               axis.text = element_text(size = 7),
               axis.title = element_blank(),
               panel.grid.major.y = element_blank())

# # Two clusters
# c2 <- p + facet_grid(c2~.)
# c2 <- c2 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

# # Three clusters
# c3 <- p + facet_grid(c3~.)
# c3 <- c3 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

# Four clusters. Fig1A (mRNA patterns).
c4 <- p + facet_grid(c4~.)
# ggsave(file ="path",
#        plot = c4, height=10, width=4.3, units = "cm")

# # Five clusters
# c5 <- p + facet_grid(c5~.)
# c5 <- c5 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

################################################################################################
# Plot RP expression in different mRNA clusters
################################################################################################
p <- ggplot(rp_clust, aes(x = Time, y = mean_clog2rib, group = Name, 
                               color = as.factor(c4)))
p <- p + annotate(
  "rect", xmin=12, xmax=24,
  ymin=-3.3,
  ymax=3.3,
  fill="black", alpha=0.1)

p <- p + geom_line(alpha = 0.5, size = 0.1, color = "#55341A") 
p <- p + scale_x_continuous(breaks = seq(0, 24, by = 6), minor_breaks = F,
                            expand = expansion(mult = c(0.05, 0.02)))
p <- p + scale_y_continuous(limits = c(-3.3, 3.3),
                            breaks = seq(-3, 3, 1), minor_breaks = F,
                            expand = expansion(mult = c(0)))
p <- p + theme_bw()
p <- p + theme(legend.position = "none",
               line = element_line(size = 0.2),
               panel.border = element_rect(size = 0.2, fill = NA),
               strip.background = element_blank(),
               strip.text.y = element_blank(),
               axis.text = element_text(size = 7),
               axis.title = element_blank(),
               panel.grid.major.y = element_blank())

# # Two clusters
# c2 <- p + facet_grid(c2~.)
# c2 <- c2 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

# # Three clusters
# c3 <- p + facet_grid(c3~.)
# c3 <- c3 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

# Four clusters. Fig1A (ribosome patterns).
c4 <- p + facet_grid(c4~.)
# ggsave(file ="path",
#        plot = c4, height=10, width=4.3, units = "cm")

# # Five clusters
# c5 <- p + facet_grid(c5~.)
# c5 <- c5 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

################################################################################################
# Plot MS expression in different mRNA clusters
################################################################################################
p <- ggplot(ms_clust, aes(x = Time, y = mean_clog2int, group = Name,
                              color =  as.factor(c4)))
p <- p + annotate(
  "rect", xmin=12, xmax=24,
  ymin=-3.3,
  ymax=3.3,
  fill="black", alpha=0.1)

p <- p + geom_line(alpha = 0.5, size = 0.1, color = "#AE2708") 
p <- p + scale_x_continuous(breaks = seq(0, 24, by = 6), minor_breaks = F,
                            expand = expansion(mult = c(0.05, 0.02)))
p <- p + scale_y_continuous(limits = c(-3.3, 3.3),
                            breaks = seq(-3, 3, 1), minor_breaks = F,
                            expand = expansion(mult = c(0)))
p <- p + theme_bw()
p <- p + theme(legend.position = "none",
               line = element_line(size = 0.2),
               panel.border = element_rect(size = 0.2, fill = NA),
               strip.background = element_blank(),
               strip.text.y = element_blank(),
               axis.text = element_text(size = 7),
               axis.title = element_blank(),
               panel.grid.major.y = element_blank())
# # Two clusters
# c2 <- p + facet_grid(c2~.)
# c2 <- c2 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

# # Three clusters
# c3 <- p + facet_grid(c3~.)
# c3 <- c3 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

# Four clusters. Fig1A (protein patterns).
c4 <- p + facet_grid(c4~.)
# ggsave(file ="path",
#        plot = c4, height=10, width=4.3, units = "cm")

# # Five clusters
# c5 <- p + facet_grid(c5~.)
# c5 <- c5 + theme(strip.background = element_blank(),
#                  strip.text.y = element_text(angle = 0, size = 20),
#                  axis.text = element_text(size = 20),
#                  axis.title = element_text(size = 20))

################################################################################################
# Plot correlation of log2FC relative to mean (= 0) of rna vs. rp and rp vs. ms (Fig1A)
################################################################################################
rrm <- inner_join(select(rna, Name:Time_point, mean_clog2rna),
                  select(rp, Name:Time_point, mean_clog2rib),
                  by = c("Name", "Time_point"))
rrm <- inner_join(rrm,
                  select(ms, Name:Time_point, mean_clog2int))

# Remove non-cyclic genes
rrm <- filter(rrm, Name %in% rna_diff$Name)

# Scatter plots and pearson correaltion
c <- ggplot(rrm, aes(x = mean_clog2rna, y = mean_clog2rib)) +
  geom_point(color = "gray20", size =  0.5, alpha = 0.5, shape = 16) +
  geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.2, size = 0.5) +
  scale_x_continuous(limits = c(-3.3, 3.3),
                     breaks = seq(-3, 3, 1), minor_breaks = F,
                     expand = expansion(mult = c(0))) +
  scale_y_continuous(limits = c(-3.3, 3.3),
                     breaks = seq(-3, 3, 1), minor_breaks = F,
                     expand = expansion(mult = c(0))) + 
  theme_bw() +
  theme(line = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2, fill = NA),
        axis.text = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid = element_blank())

# ggsave(file ="path",
#        plot = c, height=4, width=4, units = "cm")

c <- ggplot(rrm, aes(x = mean_clog2rib, y = mean_clog2int)) +
  geom_point(color = "gray20", size =  0.5, alpha = 0.5, shape = 16) +
  geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.2, size = 0.5) +
  scale_x_continuous(limits = c(-3.3, 3.3),
                     breaks = seq(-3, 3, 1), minor_breaks = F,
                     expand = expansion(mult = c(0))) +
  scale_y_continuous(limits = c(-3.3, 3.3),
                     breaks = seq(-3, 3, 1), minor_breaks = F,
                     expand = expansion(mult = c(0))) + 
  theme_bw() +
  theme(line = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2, fill = NA),
        axis.text = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid = element_blank())

# ggsave(file ="path",
#        plot = c, height=4, width=4, units = "cm")

# Calculate correlation coefficient (Pearson's)
cor_rna_rp <- cor(rrm$mean_clog2rna, rrm$mean_clog2rib)
cor_rp_ms <- cor.test(rrm$mean_clog2rib, rrm$mean_clog2int)

################################################################################################
# Plot median relative amplitude of rna, ribosomes and protein (Fig1B)
################################################################################################
amp_med <- tibble(measure = c("rna","rib","prot"),
                  value = c(
                    median(amp_red$mRNA_amp),
                    median(amp_red$rib_amp),
                    median(amp_red$protein_amp)
                  )) %>% mutate(measure = factor(measure, levels = c("rna", "rib", "prot")))

a <- ggplot(amp_med, aes(x = measure, y = value, fill = measure)) +
  geom_col(position="dodge", width = 0.5) +
  scale_y_continuous(limits = c(0, 1.7),
                     breaks = seq(-0, 1.7, 0.2), minor_breaks = F,
                     expand = expansion(mult = c(0))) +
  scale_fill_manual(values = c("#216b21ff","#55341A", "#AE2708")) +
  theme_bw() +
  theme(legend.position = "none",
        line = element_line(size = 0.5),
        panel.border = element_rect(size = 0.2, fill = NA),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank())

# ggsave(file ="path",
#        plot = a, height=7.8, width=3.8, units = "cm")

################################################################################################
# Plot protein amplitude vs. daily synthesis to daily protein abundace ratio (Fig1C)
################################################################################################
cor.test(amp_red$protein_amp, amp_red$day_avg_rib - amp_red$day_avg_prot)

k <- ggplot(amp_red, aes(x = day_avg_rib - day_avg_prot, y = protein_amp)) +
  geom_point(color = "gray20", size =  0.5, alpha = 0.5, shape = 16) +
  scale_x_continuous(limits = c(-5, 7),
                     breaks = seq(-5, 7, 1), minor_breaks = F,
                     expand = expansion(mult = c(0))) +
  scale_y_continuous(limits = c(0, 4.3),
                     breaks = seq(0, 4, 1), minor_breaks = F,
                     expand = expansion(mult = c(0))) + 
  theme_bw() +
  theme(line = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2, fill = NA),
        axis.text = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid = element_blank())

# ggsave(file ="path",
#        plot = k, height=4, width=4, units = "cm")