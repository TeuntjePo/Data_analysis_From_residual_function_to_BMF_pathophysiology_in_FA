#### Teuntje Poortvliet
#### Single nucleotide variant and indel analysis

### Preparation for analysis

## Downloading and loading the developer version of Mutationalpatterns and all dependencies
install.packages("devtools")
install.packages("tidyr")
install.packages('dplyr')
install.packages('reshape2')
install.packages('ggalluvial')
install.packages("gridExtra")
devtools::load_all('/X/X/X/MutationalPatterns')

## Installing and loading the reference genome
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Loading required packages
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)

## Loading VCF files (AHH1-G2 = FANCC WT, 2B3-3E5 = FANCC K0 & 2B3-20F3-24A7 + APH = FANCC KO + APH)
input_folder <- "X/AHH1_WT_&_FANCC_KO/"

vcf_files <- list.files(input_folder, pattern = ".vcf", full.names = TRUE)
vcf_files <- vcf_files[c(3, 1, 2)]
sample_names <- c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH")
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type='all')

### Analysis SNVs and indels using MutationalPatterns

## SNV analysis

# SBS contributions
snv <- get_mut_type(vcfs, type = "snv")
mut_type_occurrences <- mut_type_occurrences(snv, ref_genome)

plot_spectrum <- plot_spectrum(mut_type_occurrences, by = sample_names, CT = TRUE, error_bars =  "none")

# 96 mutation profile 
mut_mat <- mut_matrix(vcf_list = snv, ref_genome = ref_genome)

profile_96_contribution <- plot_96_profile(mut_mat, ymax = 0.05) 
                          + scale_y_continuous(breaks = c(0, 0.05), labels = c("0", "0.05"))

# Heatmap relative contributions with extensive context
mut_mat_ext_context <- mut_matrix(snv, ref_genome, extension = 2)
heatmap_mut_mat_ext_context <- plot_profile_heatmap(mut_mat_ext_context, by = sample_names)

## Signature SBS analysis

# SBS signature contribution
signatures = get_known_signatures()
refit_snv <- fit_to_signatures(mut_mat, signatures)

signature_contribution_refit_plot <- plot_contribution(refit_snv$contribution,
                                                       coord_flip = FALSE,
                                                       mode = "absolute")

# SBS signature contribution with strict refit
strict_refit_snv <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
fit_res_strict <- strict_refit_snv$fit_res

signature_contribution_strict_refit_plot <- plot_contribution(fit_res_strict$contribution,
                                                              coord_flip = FALSE,
                                                              mode = "absolute")

# Bootstrapping signatures
contribution_bootstrap <- fit_to_signatures_bootstrapped(mut_mat,
                                                        signatures,
                                                        method = "strict")

plot_contribution_bootstrap <- plot_bootstrapped_contribution(contribution_bootstrap,
                               mode = "relative",
                               plot_type = "dotplot")

# Selection most different contributions
plot_contribution_selected_signatures <- plot_bootstrapped_contribution(contribution_bootstrap[, c(5, 10, 18, 21, 22, 28, 39, 41, 42, 43)],
                                                                        mode = "relative",
                                                                        plot_type = "dotplot")

# SBS8 signature 
SBS8 <- as.matrix(signatures[, 11])
names <- rownames(mut_mat)
row.names(SBS8) <- names

plot_96_profile_SBS8 <- plot_96_profile(SBS8, ymax = 0.1) + theme(strip.text.y = element_blank())

# Strict refit reconstruction plot
signature_refit_reconstruction_plot <- plot_original_vs_reconstructed(mut_mat, refit_snv$reconstructed, 
                                                                      y_intercept = 0.95) +
                                                                      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

# Contribution heatmap signatures
signature_contribution_heatmap <- plot_contribution_heatmap(refit_snv$contribution, 
                                                            cluster_samples = TRUE)

# Similarity between mutational profiles and signatures
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, signatures)

Cosine_heatmap_signatures <- plot_cosine_heatmap(cos_sim_samples_signatures, 
                                                 cluster_rows = TRUE, cluster_cols = TRUE)
## Indel analysis

# Indel spectrum 
indels <- get_mut_type(vcfs, type = "indel")
indel <- get_indel_context(indels, ref_genome)
indel_count <- count_indel_contexts(indel)
indel_counts <- indel_count[, c(3, 2, 1)]

indel_spectrum <- plot_main_indel_contexts(indel_counts) + theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
plot_indel_spectrum <- indel_spectrum + scale_y_continuous(limits = c(0, 50))

## Analysis mutation accumulation

# Mutation accumulation per day for SNVs
Mutation_accumulation_3E5_58 <- length(snv$`2B3-3E5`)/58
Mutation_accumulation_24A7 <- length(snv$`2B3-20F3-24A7 + APH`)/(58+43)
Mutation_accumulation_AHH1_G2 <- length(snv$`AHH1-G2`)/45

Mutation_accumulation_df <- data.frame(sample_type=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"),
                                       mutations_per_day=c(7.177778, 5.913793, 5.871287))
Mutation_accumulation_df$sample_type <- factor(Mutation_accumulation_df$sample_type, levels = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

plot_mutation_accumulation <- ggplot(data=Mutation_accumulation_df, aes(x=sample_type, y=mutations_per_day, fill=sample_type)) +
                              geom_bar(stat="identity")+ theme_minimal()

plot_mutation_accumulation_per_day <- plot_mutation_accumulation + scale_fill_manual(values=c("mediumseagreen", "dodgerblue2", "coral2")) +
                                      ggtitle("Mutations per day") +
                                      ylab("Number of mutations per day") +
                                      theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") +
                                      labs(fill = "Sample") +
                                      scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH")) +
                                      scale_y_continuous(breaks = c(0, 2, 4, 6, 8), labels = c("0", "2", "4", "6", "8"), limits = c(0,8))

# Mutation accumulation per day for SNVs, indels and DSBs together
Mutation_accumulation_3E5_58_all <- length(vcfs$`2B3-3E5`)/58
Mutation_accumulation_24A7_all <- length(vcfs$`2B3-20F3-24A7 + APH`)/(58+43)
Mutation_accumulation_AHH1_G2_all <- length(vcfs$`AHH1-G2`)/(45)

Mutation_accumulation_df_all <- data.frame(sample_type=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"),
                                           mutations_per_day=c(8.088889, 6.517241, 6.831683))
Mutation_accumulation_df_all$sample_type <- factor(Mutation_accumulation_df_all$sample_type, levels = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))


plot_mutation_accumulation_all <- ggplot(data=Mutation_accumulation_df_all, aes(x=sample_type, y=mutations_per_day, fill=sample_type)) +
                                  geom_bar(stat="identity")+ theme_minimal()

plot_mutation_accumulation_all_per_day <- plot_mutation_accumulation_all + scale_fill_manual(values=c("mediumseagreen", "dodgerblue2", "coral2")) +
                                          ggtitle("Mutation accumulation per day all mutations") +
                                          ylab("Mutations per day") +
                                          theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), legend.position = "none")+
                                          labs(fill = "Sample") +
                                          scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

### Analysis deletions with microhomology

## Coördinates of microhomology-mediated deletions for manual validation in IGV
indels <- get_mut_type(vcfs, type = "indel")
indel <- get_indel_context(indels, ref_genome)
indel_counts <- count_indel_contexts(indel)

indels_try <- as.data.frame(indel)
Microhomology_variants <- indels_try[grep('*microhomology', indels_try$muttype),]

APH_treated_microhom <- Microhomology_variants[Microhomology_variants$group_name== "2B3-20F3-24A7 + APH",] 
Coordinates_APH_treated_microhom <- APH_treated_microhom[, c("seqnames", "start", "end", "muttype")]

clone_2B3_3E5_microhom <- Microhomology_variants[Microhomology_variants$group_name== "2B3-3E5",] 
Coordinates_2B3_3E5_microhom <- clone_2B3_3E5_microhom[, c("seqnames", "start", "end", "muttype")]

AHH1_G2_microhom <- Microhomology_variants[Microhomology_variants$group_name== "AHH1-G2",] 
Coordinates_AHH1_G2_microhom <- AHH1_G2_microhom[, c("seqnames", "start", "end", "muttype")]

## Analysis validated microhomology-mediated indels

# Loading validated microhomology-mediated deletions
file_indels_microhom <- read.csv('X/20212904_Indels_with_microhomology_FANCC.csv', sep = ';', header = TRUE, fill = TRUE)

# Number of microhomology-mediated deletions per day
Number_indels_microhom_APH_per_day <- sum(file_indels_microhom$clone == '2B3-20F3-24A7 + APH')/(58+43)
Number_indels_microhom_2B3_3E5_per_day <- sum(file_indels_microhom$clone == '2B3-3E5')/(58)
Number_indels_microhom_AHH1_G2_per_day <- sum(file_indels_microhom$clone == 'AHH1-G2')/(45)

Indels_per_day <- data.frame(Sample=c('2B3-20F3-24A7 + APH', '2B3-3E5', 'AHH1-G2'), Indels_per_day = c(0.4356436,  0.1551724, 0.08888889))
Indels_per_day$Sample <- factor(Indels_per_day$Sample, levels = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

plot_mhm_deletions_per_day <- ggplot(data=Indels_per_day, aes(x=Sample, y=Indels_per_day, fill=Sample)) +
                              geom_bar(stat="identity") + 
                              theme_minimal() +
                              theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), legend.position = "none")+
                              ylab("Number of mhm-deletions per day") +
                              ggtitle("Number of deletions with microhomology per day") +
                              scale_x_discrete(limits = c("AHH1-G2","2B3-3E5", "2B3-20F3-24A7 + APH")) +
                              scale_fill_manual(values=c("mediumseagreen","dodgerblue2", "coral2")) 

# Percentage of microhomology-mediated deletions located in genes
Indels_microhom_in_genes <- is.na(file_indels_microhom$gene) != TRUE
df_indels_microhom_in_genes <- mutate(file_indels_microhom, In_gene = Indels_microhom_in_genes)

plot_mhm_deletions_in_genes <- ggplot(data=df_indels_microhom_in_genes, aes(x=clone)) +
                              geom_bar(aes(fill=In_gene), position = 'fill', stat = "count") +
                              theme_minimal() +
                              theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
                              ylab("Percentage") +
                              labs(fill = "In gene") +
                              ggtitle("Percentage of indels with microhomology located in gene") +
                              scale_x_discrete(limits = c("AHH1-G2","2B3-3E5", "2B3-20F3-24A7 + APH")) +
                              scale_y_continuous(labels = c("0", "25", "50", "75", "100"))

plot_number_mhm_deletions_in_genes <- ggplot(data=df_indels_microhom_in_genes, aes(x=clone)) +
                                      geom_bar(aes(fill=In_gene), position = 'stack', stat = "count") +
                                      theme_minimal() +
                                      theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank())+
                                      ylab("Number of indels with microhomology") +
                                      labs(fill = "In gene") +
                                      ggtitle("Number of indels with microhomology located in gene") +
                                      scale_x_discrete(limits = c("AHH1-G2","2B3-3E5", "2B3-20F3-24A7 + APH")) +
                                      scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), labels = c("0", "5", "10", "15", "20", "25", "30", "35", "40", "45"), limits = c(0, 45))

# Making dataframe classifying indels like in Zou et al. (2018) paper

df_indels <- indels_try[, c("group_name", "seqnames", "start", "end", "width", "muttype")]
df_indels_2 <- df_indels %>% mutate(indel_type = c("none")) %>% mutate(class = c("none"))

df_indels_2$group_name[df_indels_2$group_name == "2B3_20F3_24A7_APH"] <- '2B3-20F3-24A7 + APH' 
df_indels_2$group_name[df_indels_2$group_name == "2B3_3E5"] <- '2B3-3E5' 
df_indels_2$group_name[df_indels_2$group_name == "AHH1_G2"] <- 'AHH1-G2' 

df_indels_2$indel_type[df_indels_2$muttype == 'C_deletion' | df_indels_2$muttype == 'T_deletion' ] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == '2bp_deletion'] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == '3bp_deletion'] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == '4bp_deletion'] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == '5bp_deletion'] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == '8bp_deletion'] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == '14bp_deletion'] <- "Other Del."
df_indels_2$indel_type[df_indels_2$muttype == 'C_insertion' | df_indels_2$muttype == 'T_insertion'] <- "Insertion=1bp"
df_indels_2$indel_type[grepl('bp_insertion', df_indels_2$muttype)] <- "Insertion>=2bp" 
df_indels_2$indel_type[grepl('with_microhomology', df_indels_2$muttype) & df_indels_2$width %in% 0:6] <- "Mh-mediated Del.0-5bp"
df_indels_2$indel_type[grepl('with_microhomology', df_indels_2$muttype) & df_indels_2$width %in% 7:11] <- "Mh-mediated Del.6-10bp"
df_indels_2$indel_type[grepl('with_microhomology', df_indels_2$muttype) & df_indels_2$width %in% 12:21] <- "Mh-mediated Del.11-20bp"
df_indels_2$indel_type[grepl('with_microhomology', df_indels_2$muttype) & df_indels_2$width %in% 22:31] <- "Mh-mediated Del.21-30bp"
df_indels_2$indel_type[grepl('with_microhomology', df_indels_2$muttype) & df_indels_2$width %in% 32:41] <- "Mh-mediated Del.31-40bp"

df_indels_2$class[grepl('Mh-mediated', df_indels_2$indel_type)] <- "Mh-mediated deletions"
df_indels_2$class[grepl('Other', df_indels_2$indel_type)] <- "Other deletions"
df_indels_2$class[grepl('Insertion', df_indels_2$indel_type)] <- "Insertions"

df_indels_2_APH <- df_indels_2[df_indels_2$group_name == '2B3-20F3-24A7 + APH',]
df_indels_2_3E5 <- df_indels_2[df_indels_2$group_name == '2B3-3E5',]
df_indels_2_AHH1 <- df_indels_2[df_indels_2$group_name == 'AHH1-G2',]

# Plot types of indels in own samples
plot_type_indels_samples <- ggplot(data=df_indels_2, aes(x=group_name)) + geom_bar(aes(fill=class), position = 'fill', stat = "count") +
                            theme_minimal() +
                            ggtitle("Types of indels in samples ") +
                            ylab("Percentage") +
                            labs(fill = "Sample") +
                            theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
                            scale_fill_manual(values=c("steelblue1", "violetred2", "gold1")) +
                            scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1.0), labels = c("0", "25","50", "75", "100"), limits = c(0, 1)) +
                            scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

## Analysis microhomology-mediated indels comparison with Zou et al. (2018)

# Loading Zou et al. (2018) data and making dataframe
input_file_zou <- "X/denovo_indels.txt"

raw_ind_data_zou <- read.delim(file = input_file_zou, stringsAsFactors = F)
raw_ind_data_zou_FANCC <- raw_ind_data_zou[raw_ind_data_zou$knockout == "FANCC",]
ind_data_zou_wanted_cols <- raw_ind_data_zou_FANCC[, c(1, 4:9, 12:26)]

ind_data_zou_wanted_cols$indeltype[ind_data_zou_wanted_cols$indeltype != "Microhomology-mediated" & ind_data_zou_wanted_cols$indeltype != 'Ins'] <- "Other deletions"
ind_data_zou_wanted_cols$indeltype[ind_data_zou_wanted_cols$indeltype == "Microhomology-mediated"] <- 'Mh-mediated deletions' 
ind_data_zou_wanted_cols$indeltype[ind_data_zou_wanted_cols$indeltype == 'Ins'] <- "Insertions"

ind_data_zou_wanted_cols_child <- ind_data_zou_wanted_cols[ind_data_zou_wanted_cols$clonelevel == "child",]

# Zou et al. (2018) indel spectrum (parental clone not included for analysis)
plot_indel_type_Zou <- ggplot(data=ind_data_zou_wanted_cols_child, aes(x=Sample)) + geom_bar(aes(fill=indeltype), position = "fill", stat = "count") +
                      theme_minimal() +
                      ggtitle("Types of indels in samples Zou ") +
                      ylab("Percentage") +
                      labs(fill = "Indel type") +
                      theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
                      scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1.0), labels = c("0", "25","50", "75", "100"), limits = c(0, 1)) +
                      scale_fill_manual(values=c("steelblue1", "violetred2", "gold1")) 

# Comparing own data and Zou et al. (2018) data

# Making combined dataframe own data and Zou et al. (2018) data
df_indels_3 <- df_indels_2
names(df_indels_3)[names(df_indels_3) == "group_name"] <- "sample"
df_indels_3 <- mutate(df_indels_3, size = df_indels_3$width - 1)

ind_data_zou_wanted_cols_names <- ind_data_zou_wanted_cols
names(ind_data_zou_wanted_cols_names)[names(ind_data_zou_wanted_cols_names) == "indeltype"] <- "class"
names(ind_data_zou_wanted_cols_names)[names(ind_data_zou_wanted_cols_names) == "Sample"] <- "sample"

df_bind_Zou_own_indels <- rbind(
  data.frame(c(df_indels_3, sapply(setdiff(names(ind_data_zou_wanted_cols_names), names(df_indels_3)), function(x) NA))),
  data.frame(c(ind_data_zou_wanted_cols_names, sapply(setdiff(names(df_indels_3), names(ind_data_zou_wanted_cols_names)), function(x) NA)))
)

ind_data_zou_wanted_cols_child_names <- ind_data_zou_wanted_cols_child
names(ind_data_zou_wanted_cols_child_names)[names(ind_data_zou_wanted_cols_child_names) == "indeltype"] <- "class"
names(ind_data_zou_wanted_cols_child_names)[names(ind_data_zou_wanted_cols_child_names) == "Sample"] <- "sample"

df_bind_Zou_child_own_indels <- rbind(
  data.frame(c(df_indels_3, sapply(setdiff(names(ind_data_zou_wanted_cols_child_names), names(df_indels_3)), function(x) NA))),
  data.frame(c(ind_data_zou_wanted_cols_child_names, sapply(setdiff(names(df_indels_3), names(ind_data_zou_wanted_cols_child_names)), function(x) NA)))
)

# Plots own and Zou indel data
plot_indel_type_number_both <- ggplot(data= df_bind_Zou_child_own_indels, aes(x=sample)) + geom_bar(aes(fill=class), position = "stack", stat = "count") +
                              theme_minimal() +
                              ggtitle("Types of indels in samples own and Zou") +
                              ylab("Number of events") +
                              labs(fill = "Indel type") +
                              theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
                              scale_fill_manual(values=c("steelblue1", "violetred2", "gold1")) +
                              scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", "", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3",
                                                        "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

plot_indel_type_percentage_both <- ggplot(data= df_bind_Zou_child_own_indels, aes(x=sample)) + geom_bar(aes(fill=class), position = "fill", stat = "count") +
                                  theme_minimal() +
                                  ggtitle("Types of indels in samples own and Zou") +
                                  ylab("Percentage") +
                                  labs(fill = "Indel type") +
                                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
                                  scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1.0), labels = c("0", "25","50", "75", "100"), limits = c(0, 1)) + 
                                  scale_fill_manual(values=c("steelblue1", "violetred2", "gold1")) +
                                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", "", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3",
                                                            "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

# Dataframe percentage in gene for own data and Zou et al. data
Indels_microhom_in_genes_own <- is.na(file_indels_microhom$gene) != TRUE
df_indels_microhom_in_genes_own <- mutate(file_indels_microhom, In_gene = Indels_microhom_in_genes_own)

Indels_microhom_in_genes_Zou <- ind_data_zou_wanted_cols_names$Gene != "-"
df_indels_microhom_in_genes_Zou <- mutate(ind_data_zou_wanted_cols_names, In_gene = Indels_microhom_in_genes_Zou) 
df_indels_microhom_in_genes_Zou <- df_indels_microhom_in_genes_Zou[grepl("Mh-mediated", df_indels_microhom_in_genes_Zou$class), c(2, 4, 5, 8, 13, 16, 23)]

names(df_indels_microhom_in_genes_Zou)[names(df_indels_microhom_in_genes_Zou) == "Gene"] <- "gene"
names(df_indels_microhom_in_genes_Zou)[names(df_indels_microhom_in_genes_Zou) == "sample"] <- "clone"
names(df_indels_microhom_in_genes_Zou)[names(df_indels_microhom_in_genes_Zou) == "Chrom"] <- "chrom"

df_bind_Zou_own_mhm <- rbind(
  data.frame(c(df_indels_microhom_in_genes_own, sapply(setdiff(names(df_indels_microhom_in_genes_Zou), names(df_indels_microhom_in_genes_own)), function(x) NA))),
  data.frame(c(df_indels_microhom_in_genes_Zou, sapply(setdiff(names(df_indels_microhom_in_genes_own), names(df_indels_microhom_in_genes_Zou)), function(x) NA)))
)

df_bind_Zou_own_mhm$In_gene[df_bind_Zou_own_mhm$In_gene == "TRUE"] <- "Coding"
df_bind_Zou_own_mhm$In_gene[df_bind_Zou_own_mhm$In_gene == "FALSE"] <- "Non-coding"

df_bind_Zou_own_mhm$In_gene <- factor(df_bind_Zou_own_mhm$In_gene, levels = c("Non-coding", "Coding"))

# Plot percentage of microhomology-mediated deletions in genes
percentage_mhm_deletions_genes <- ggplot(data= df_bind_Zou_own_mhm, aes(x=clone)) +
                                  geom_bar(aes(fill=In_gene), position = 'fill', stat = "count") +
                                  theme_minimal() +
                                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
                                  ylab("Percentage") +
                                  labs(fill = "Coding status") +
                                  ggtitle("") +
                                  scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1.0), labels = c("0", "25","50", "75", "100"), limits = c(0, 1)) +
                                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", " ","HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                            "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

# Dataframe with information on microhomology-mediated located in which gene for percentage in transcribed region
df_bind_Zou_own_mhm$gene[df_bind_Zou_own_mhm$gene == "-"] <- NA
list_genes <- df_bind_Zou_own_mhm$gene[is.na(df_bind_Zou_own_mhm$gene) != TRUE]

file_genes_indels <- "X/mart_export_indels.txt"
df_genes_indels <- read.delim(file_genes_indels, header = TRUE, sep = "", stringsAsFactors = F)
genes_indels <- df_genes_indels$Gene.ID

izo <- full_join(df_bind_Zou_own_mhm, df_genes_indels)
names(izo)[names(izo) == "Gene.ID"] <- "gene_id"
izo <- izo[, c(1, 2, 3, 4, 5, 6, 8, 10)] 

# Dataframe with gene and expression for microhomology-mediated deletions 
#(Expression data comes from the B-lymphocyte cell line GM12878 total RNA-seq dataset ENCSR820PHH generated by the Barbara Wold lab at Caltech and was downloaded from the Encode portal (https://www.encodeproject.org/)) 

input_data_gene_expression_GM12878 <- "X/ENCFF873VWU.tsv" 
data_gene_expression <- read.delim(input_data_gene_expression_GM12878, 
                                   header = TRUE, sep = "", dec = ".", stringsAsFactors = F) 
data_expression <- data_gene_expression[650:59429,]

new_data_expression <- mutate(data_expression, found = c("FALSE"))
data_expression_2 <- separate(new_data_expression, gene_id, into=c('gene_id','gene_id_version'), sep = '\\.')
data_expression_3 <- data_expression_2[, c(1, 7, 8, 19)]
data_expression_3$found[data_expression_3$gene_id %in% genes_indels] <- "TRUE"

epdi <- left_join(izo, data_expression_3)

# Calculating percentages transcribed
percentages_transcribed <- vector()
samples <- c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", "HAP1_FANCC_27", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
             "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11")

for(i in samples) {               
  percentage <- sum(epdi$TPM != 0.00 & epdi$clone == i, na.rm = TRUE)/ (sum(epdi$TPM == 0.00 & epdi$clone == i , na.rm = TRUE) + sum(is.na(epdi$TPM & epdi$clone == i)) + sum(epdi$TPM != 0.00 & epdi$clone == i, na.rm = TRUE)) * 100 
  percentages_transcribed <- c(percentages_transcribed, percentage)    
}

# Dataframe and plot percentage of microhomology-mediated deletions in transcribed regions
transcribed_mhm_deletions <- data.frame(sample = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH",
                                         "HAP1_FANCC_27" ,"HAP1_FANCC_27-1", "HAP1_FANCC_27-3", 
                                         "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", "HAP1_FANCC_27-6", 
                                         "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"),
                                          Transcribed = percentages_transcribed,
                                          Not_transcribed = 100 - percentages_transcribed) 

df_transcribed_mhm <- rbind(
  data.frame(sample = transcribed_mhm_deletions$sample, percentage = transcribed_mhm_deletions$Transcribed , type = "Transcribed"),
  data.frame(sample = transcribed_mhm_deletions$sample, percentage = transcribed_mhm_deletions$Not_transcribed, type ="Non-transcribed")
)

df_transcribed_mhm$type <- factor(df_transcribed_mhm$type, levels = c("Non-transcribed", "Transcribed"))

plot_mhm_deletions_transcribed <- ggplot(data=df_transcribed_mhm, aes(x=sample, y = percentage)) + geom_bar(aes(fill=type), position = 'fill', stat = "identity") +
                                  theme_minimal() +
                                  ggtitle("Percentage of mhm-indels in transcribed region") +
                                  ylab("Percentage") +
                                  labs(fill = "Transcription status") +
                                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
                                  scale_fill_manual(values=c("firebrick2", "cyan3")) +
                                  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), labels = c("0", "25", "50", "75", "100")) +
                                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", " ", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                            "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

