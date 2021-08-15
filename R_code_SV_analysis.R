#### Teuntje Poortvliet
#### Structural Variant analysis

### Preparation for analysis

## Loading required packages and reference genome 
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

## Loading own data
input_data <- "X/SVs_24A7_APH_3E5.csv" 
data_frame_own <- read.csv(input_data, header = TRUE, sep = ';', fill = TRUE)

### SV analysis own data (AHH1-G2 = FANCC WT, 2B3-3E5 = FANCC KO & 2B3-20F3-24A7 + APH = FANCC KO + APH)

## Number of events per sample in total
Number_events_2B3_20F3_24A7_APH <- sum(data_frame_own$clone == "2B3-20F3-24A7 + APH")
Number_events_2B3_3E5 <- sum(data_frame_own$clone == "2B3-3E5")
Number_events_AHH1 <- 0

Number_events_clones <- c(28, 5, 0)
df_number_events <- data_frame(Sample=c("2B3-20F3-24A7 + APH", "2B3-3E5", "AHH1-G2"), 
                               number_of_events=Number_events_clones)
df_number_events$Sample <- factor(df_number_events$Sample, levels = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

plot_number_events <- ggplot(data=df_number_events, aes(x=Sample, y=number_of_events, fill=Sample)) +
                      geom_bar(stat="identity") + 
                      theme_minimal()

plot_total_events <- plot_number_events + scale_fill_manual(values=c("mediumpurple1", "dodgerblue2", "coral2")) +
                    ggtitle("Total SV events per sample") +
                    ylab("Number of SV events") +
                    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") +
                    scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30), labels = c("0", "5", "10", "15", "20", "25", "30"), limits=c(0, 30)) +
                    scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

## Plot number of events per sample per day
Number_events_2B3_20F3_24A7_APH_day <- sum(data_frame_own$clone == "2B3-20F3-24A7 + APH")/ (58+43)
Number_events_2B3_3E5_day <- sum(data_frame_own$clone == "2B3-3E5")/ 58
Number_events_AHH1_day <- 0

Number_events_clones_day <- c(0.28, 0.09, 0)
df_number_events_day <- data_frame(Sample=c("2B3-20F3-24A7 + APH", "2B3-3E5", "AHH1-G2"), 
                                   number_of_events_day=Number_events_clones_day)
df_number_events_day$Sample <- factor(df_number_events_day$Sample, levels = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

plot_number_events_day <- ggplot(data=df_number_events_day, aes(x=Sample, y=number_of_events_day, fill=Sample)) +
                          geom_bar(stat="identity") + 
                          theme_minimal()

plot_sv_per_day <- plot_number_events_day_comp <- plot_number_events_day + scale_fill_manual(values=c("mediumpurple1", "dodgerblue" ,"coral2")) +
                  ggtitle("SV events per sample per day") +
                  xlab("Sample") + ylab("Number of SV events per day") +
                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
                  scale_y_continuous(breaks=c(0, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30), labels = c("0", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"), limits=c(0, 0.30)) +
                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))

## Plot type events per sample
plot_number_event_types <- ggplot(data=data_frame_own, aes(x=clone)) + geom_bar(aes(fill=svclass), position = 'fill', stat = "count") +
                          theme_minimal() +
                          ggtitle("SV types per sample") +
                          ylab("Contribution per event type") +
                          labs(fill = "SV type") +
                          theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +  
                          scale_y_continuous(breaks = c(0, 0.25, .50, .75, 1.0), labels = c("0", "25","50", "75", "100"), limits = c(0, 1)) +
                          scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH")) +
                          scale_fill_manual(values=c("coral2", "gold3", "limegreen", "darkturquoise", "cornflowerblue", "magenta"))

plot_number_event_types_numbers <- ggplot(data=data_frame_own, aes(x=clone)) + geom_bar(aes(fill=svclass), position = 'stack', stat = "count") +
                                  theme_minimal() + 
                                  ggtitle("SV types per sample") +
                                  ylab("Count different event types") +
                                  labs(fill = "SV type") +
                                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +  
                                  scale_y_continuous(breaks = c(0, 10, 20, 30), labels = c("0", "10","20", "30"), limits = c(0, 30)) +
                                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH")) +
                                  scale_fill_manual(values=c("coral2", "gold3", "limegreen", "darkturquoise", "cornflowerblue", "magenta"))

## Plot size SVs
plot_SV_size <- ggplot(data_frame_own, aes(x=clone, y=size, fill=clone)) +
                geom_boxplot() +
                theme_minimal() +
                ggtitle("SV sizes in the different samples") +
                theme(plot.title = element_text(hjust = 0.5)) + 
                ylab("SV sizes") +
                labs(fill = "Sample") +
                scale_fill_manual(values=c("mediumpurple1", "coral2", "dodgerblue2")) + 
                theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), legend.position = "none") +  
                scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH"))  +
                scale_y_continuous(breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 500000), labels = c("0" ,"1","10", "100", "1K", "10K", "100K", "500K"), limits=c(0.1, 500000), trans = "log10") +
                geom_jitter(shape=16, position=position_jitter(0.2))

## Percentage SVs in genes
vector_in_gene <- is.na(data_frame_own$gene) != TRUE
gene_df <- mutate(data_frame_own, In_gene = vector_in_gene)

plot_percentage_in_gene <- ggplot(data=gene_df, aes(x=clone)) + geom_bar(aes(fill=In_gene), position = 'fill', stat = "count") +
                          theme_minimal() +
                          ggtitle("Percentage of events in gene") +
                          ylab("Percentage of events") +
                          labs(fill = "In gene") +
                          scale_y_continuous(breaks=c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75", "100"), limits=c(0, 1.0)) +
                          theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +  
                          scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH")) 

plot_number_in_gene <- ggplot(data=gene_df, aes(x=clone)) + geom_bar(aes(fill=In_gene), position = 'stack', stat = "count") +
                      theme_minimal() +
                      ggtitle("Number of events in gene") +
                      ylab("Number of events") +
                      labs(fill = "In gene") +
                      theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +  
                      scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30), labels = c("0", "5", "10", "15", "20", "25", "30"), limits=c(0, 30)) +
                      scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH")) 

### SV analysis comparison with Zou et al. (2018) data

## Loading data Zou et al. and preparing dataframe with Zou et al. and own data
input_file <- "X/denovo_rgs.txt"
raw_data_zou <- read.delim(file = input_file, stringsAsFactors = F)

df_wanted_rows_zou <- raw_data_zou[grep("*FANCC*", raw_data_zou$knockout),]
wanted_cols <- c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 16, 17, 18)
df_wanted_rows_cols_zou <- df_wanted_rows_zou[, wanted_cols]
names(df_wanted_rows_cols_zou)[names(df_wanted_rows_cols_zou) == "bkdist"] <- "size"

data_frame_own_names <- data_frame_own
names(data_frame_own_names)[names(data_frame_own_names) == "clone"] <- "sample"
names(data_frame_own_names)[names(data_frame_own_names) == "chrom"] <- "chr1"

df_bind <- rbind(
  data.frame(c(data_frame_own_names, sapply(setdiff(names(df_wanted_rows_cols_zou), names(data_frame_own_names)), function(x) NA))),
  data.frame(c(df_wanted_rows_cols_zou, sapply(setdiff(names(data_frame_own_names), names(df_wanted_rows_cols_zou)), function(x) NA)))
)

## Plot SV types own data and Zou et al data combined (excluding parental clone)
all_samples_zou <- unique(df_wanted_rows_cols_zou$sample)
order_samples_zou <- c("HAP1_FANCC_27", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                       "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11")

plot_types_in_samples_both <- ggplot(data=df_bind, aes(x=sample)) + geom_bar(aes(fill=svclass), position = 'stack', stat = "count") +
                              theme_minimal() +
                              ggtitle("Types of events in samples ") +
                              ylab("Number of events") +
                              labs(fill = "SV type") +
                              theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
                              scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", " ", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                        "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11")) +
                              scale_y_continuous(breaks=c(0, 10, 20, 30, 40), labels = c("0", "10", "20", "30", "40"), limits=c(0, 45)) +
                              scale_fill_manual(values=c("coral2", "gold3", "limegreen", "darkturquoise", "cornflowerblue", "magenta")) 

## Plot number of SVs per day
SVs_Zou_per_sample_per_day <- c(sum(df_bind$sample == "HAP1_FANCC_27"), sum(df_bind$sample == "HAP1_FANCC_27-1"), sum(df_bind$sample == "HAP1_FANCC_27-3"),
                                sum(df_bind$sample == "HAP1_FANCC_27-4"), sum(df_bind$sample == "HAP1_FANCC_27-5"), sum(df_bind$sample == "HAP1_FANCC_27-6"),
                                sum(df_bind$sample == "HAP1_FANCC_27-8"), sum(df_bind$sample == "HAP1_FANCC_27-11"))/30
SVs_own_per_day <- c(0, 0.0862069, 0.2772277)

df_number_events_day_both <- data_frame(Sample=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", "HAP1_FANCC_27", "HAP1_FANCC_27-1", 
                                                 "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                 "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"), 
                                        number_of_events_day=c(SVs_own_per_day, SVs_Zou_per_sample_per_day))
df_number_events_day_both$Sample <- factor(df_number_events_day_both$Sample, levels = c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", " ", "HAP1_FANCC_27", "HAP1_FANCC_27-1", 
                                                                                        "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                                                        "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

plot_number_SVs_in_samples_both <- ggplot(data=df_number_events_day_both, aes(x=Sample, y=number_of_events_day, fill=Sample)) +
                                  geom_bar(stat="identity") + 
                                  theme_minimal() +
                                  ggtitle("Number of SVs per day ") +
                                  ylab("Number of events") +
                                  labs(fill = "Sample") +
                                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
                                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", " ", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                            "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

### Analysis SVs of Zou et al. (2018) located in genes and SVs located in transcribed regions both

## Loading data gene expression (data comes from Encode: https://www.encodeproject.org/experiments/ENCSR820PHH/)
input_data_gene_expression_GM12878 <- "X/ENCFF873VWU.tsv"
data_gene_expression <- read.delim(input_data_gene_expression_GM12878, 
                                   header = TRUE, sep = "", dec = ".", stringsAsFactors = F) 
data_expression <- data_gene_expression[650:59429,]

## Coupling gene names to Ensembl IDs
# Getting names of genes affected and using Biomart online to get the Ensembl ID's
List_genes_affected_APH <- data_frame_own$gene[is.na(data_frame_own$gene) != TRUE & data_frame_own$clone == "2B3-20F3-24A7 + APH"] 
List_genes_affected_2B3_3E5 <- data_frame_own$gene[is.na(data_frame_own$gene) != TRUE & data_frame_own$clone == "2B3-3E5"]

vector_ensembl_IDs_APH <- c('ENSG00000151468', 	
                            'ENSG00000120594',
                            'ENSG00000130021',	
                            'ENSG00000152611',	
                            'ENSG00000173442',	
                            'ENSG00000067798',	
                            'ENSG00000101336',	
                            'ENSG00000189283',	
                            'ENSG00000169744',	
                            'ENSG00000137831',	
                            'ENSG00000153310',	
                            'ENSG00000183230',	
                            'ENSG00000186153',	
                            'ENSG00000078328',	
                            'ENSG00000070366',	
                            'ENSG00000170160',	
                            'ENSG00000114279',	
                            'ENSG00000181234')

vector_ensembl_IDs_3E5 <- c('ENSG00000157766',
                            'ENSG00000170325',	
                            'ENSG00000177084')

file_APH <- "X/mart_export_APH.txt"
df_names_id_genes_APH <- read.delim(file_APH, header = TRUE, sep = "", stringsAsFactors = F)

file_3E5 <- "X/mart_export_3E5.txt"
df_names_id_genes_3E5 <- read.delim(file_3E5, header = TRUE, sep = "", stringsAsFactors = F)

file_Zou <- "X/mart_export_Zou.txt"
df_names_id_genes_Zou <- read.delim(file_Zou, header = TRUE, sep = "", stringsAsFactors = F)
vector_ensembl_IDs_Zou <- df_names_id_genes_Zou$gene_id

## Labelling genes in expression dataframe when affected by SV in APH-treated clone and adding gene name
new_data_expression <- mutate(data_expression, found = c("FALSE"))
data_expression_3 <- separate(new_data_expression, gene_id, into=c('gene_id','gene_id_version'), sep = '\\.')
data_expression_3$found[data_expression_3$gene_id %in% vector_ensembl_IDs_APH] <- "TRUE"

ranked_data_expression_3 <- mutate(data_expression_3, ranked = dense_rank(data_expression_3$TPM))
expression_data_2B3_23F3_24A7_APH <- left_join(ranked_data_expression_3, df_names_id_genes_APH)

## Labelling genes in expression dataframe when affected by SV in 2B3-3E5 FANCC KO clone and adding gene name
new_data_expression <- mutate(data_expression, found = c("FALSE"))
data_expression_4 <- separate(new_data_expression, gene_id, into=c('gene_id','gene_id_version'), sep = '\\.')
data_expression_4$found[data_expression_4$gene_id %in% vector_ensembl_IDs_3E5] <- "TRUE"

ranked_data_expression_4 <- mutate(data_expression_4, ranked = dense_rank(data_expression_4$TPM))
expression_data_3E5 <- left_join(ranked_data_expression_4, df_names_id_genes_3E5)

## Labelling genes in expression dataframe when affected by SV in Zou et al. (2018) data and adding gene name
new_data_expression <- mutate(data_expression, found = c("FALSE"))
data_expression_6 <- separate(new_data_expression, gene_id, into=c('gene_id','gene_id_version'), sep = '\\.')
data_expression_6$found[data_expression_6$gene_id %in% vector_ensembl_IDs_Zou] <- "TRUE"

ranked_data_expression_6 <- mutate(data_expression_6, ranked = dense_rank(data_expression_6$TPM))
expression_data_Zou <- left_join(ranked_data_expression_6, df_names_id_genes_Zou)

## Loading SV data Zou et al.(2018) 
input_file <- "X/denovo_rgs.txt"
raw_data_zou <- read.delim(file = input_file, stringsAsFactors = F)
FANCC_Zou <- raw_data_zou[raw_data_zou$knockout == "FANCC",]
fdz <- FANCC_Zou

## Put Zou et al. (2018) data into GRanges
zou_granges1 <- GRanges(seqnames = fdz$chr1, ranges = IRanges(start = fdz$start1, end = fdz$end1))

## Load gene locations
genes <- read.table("X/X/gencode.v19.annotation.prot.gene.bed", header = F, sep = "\t")
colnames(genes) <- c("chr", "start", "end", "strand", "type", "gene", "annotation")
genes_granges <- GRanges(seqnames = gsub("chr", "", genes$chr), ranges = IRanges(start = genes$start, end = genes$end))
genome(zou_granges1) = genome(genes_granges) = 'hg19'

## Overlap the gene file with the Zou data coordinates of SVs
hits1 <- findOverlaps(query = genes_granges, subject = zou_granges1)

## Extract subject (Zou) hits and determine percentage in gene
overlap1 <- sort(subjectHits(hits1))
length(overlap1)/length(zou_granges1)

## Extract object hits and get list with affected genes
hits1_1 <- findOverlaps(query = zou_granges1 , subject = genes_granges)

overlap1_1 <- sort(subjectHits(hits1_1))
list_genes_hits_1 <- genes[overlap1_1, c(6)]

## Make dataframe with Zou et al. and own data with information whether SV is located in gene
Zou_data_SV_genes <- mutate(fdz, In_gene = "FALSE")
Zou_data_SV_genes$In_gene[overlap1] <- "TRUE"  

data_frame_own <- read.csv(input_data, header = TRUE, sep = ';', fill = TRUE)
vector_in_gene <- is.na(data_frame_own$gene) != TRUE
gene_df <- mutate(data_frame_own, In_gene = vector_in_gene)

names(gene_df)[names(gene_df) == "clone"] <- "sample"
names(gene_df)[names(gene_df) == "chrom"] <- "chr1"

df_bind_Zou_own_SV <- rbind(
  data.frame(c(gene_df, sapply(setdiff(names(Zou_data_SV_genes), names(gene_df)), function(x) NA))),
  data.frame(c(Zou_data_SV_genes, sapply(setdiff(names(gene_df), names(Zou_data_SV_genes)), function(x) NA)))
)

df_bind_Zou_own_SV$In_gene <- as.character(df_bind_Zou_own_SV$In_gene)
df_bind_Zou_own_SV$In_gene[df_bind_Zou_own_SV$In_gene == "TRUE"] <- "Coding"
df_bind_Zou_own_SV$In_gene[df_bind_Zou_own_SV$In_gene == "FALSE"] <- "Non-coding"
df_bind_Zou_own_SV$In_gene <- factor(df_bind_Zou_own_SV$In_gene, levels = c("Non-coding", "Coding"))

## Plot SVs from own data and data Zou et al. (2018) located in gene
plot_SV_percentage_in_gene <- ggplot(data=df_bind_Zou_own_SV, aes(x=sample)) + geom_bar(aes(fill=In_gene), position = 'fill', stat = "count") +
  theme_minimal() +
  ggtitle("Percentage of SVs in gene") +
  ylab("Percentage") +
  labs(fill = "Coding status") +
  scale_y_continuous(breaks=c(0, .25, .50, .75, 1.00), labels = c("0", "25", "50", "75", "100"), limits=c(0, 1.0)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(limits=c("AHH1", "2B3-3E5", "2B3-20F3-24A7 + APH", " ", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                            "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

## Making dataframe for Zou et al. (2018) data in which SV is coupled to the gene it is in
df_hits_fancc <- fdz[sort(subjectHits(hits1)),]
hits_genes <- genes[subjectHits(hits1_1), 6]

df_hits_fancc <- mutate(df_hits_fancc, genes = hits_genes) 
df_hits_fancc <- df_hits_fancc[c(1:10, 12:36, 38:114),]

df_Zou_genes_coupled <- full_join(fdz, df_hits_fancc)

df_Zou_genes_coupled$genes <- as.character(df_Zou_genes_coupled$genes)
df_Zou_genes_coupled$genes[df_Zou_genes_coupled$genes == "C7orf10"] <- "SUGCT"
df_Zou_genes_coupled$genes[df_Zou_genes_coupled$genes == "C6orf106"] <- "ILRUN"

dzg <- df_Zou_genes_coupled[, c(1, 2, 3, 11, 12, 16, 17, 19)]

## Determining the percentage of SVs located in transcribed genes for own and Zou et al. (2018) data

# Dataframe including info for gene ID, expression, SV present and affected gene with annotated if present in Zou et al. (2018) data
new_data_expression <- mutate(data_expression, found = c("FALSE"))

data_expression_1 <- separate(new_data_expression, gene_id, into=c('gene_id','gene_id_version'), sep = '\\.')
data_expression_1$found[data_expression_1$gene_id %in% vector_ensembl_IDs_Zou] <- "TRUE"

ranked_data_expression_1 <- mutate(data_expression_1, ranked = dense_rank(data_expression_1$TPM))
expression_data_Zou <- left_join(ranked_data_expression_1, df_names_id_genes_Zou)

edz <- expression_data_Zou[, c(1, 7, 8, 19, 21)]
names(edz)[names(edz) == "gene_name"] <- "genes"

## Dataframe Zou et al. (2018) including names genes affected by SV and expression data
wanted_rows_expression_data <- edz[edz$gene_id %in% vector_ensembl_IDs_Zou,]
df_all_info <- full_join(dzg, wanted_rows_expression_data)
df_all_info <- df_all_info[1:211,]

## Calculating percentages transcribed Zou
percentages_transcribed_Z <- vector()
samples_Z <- c("HAP1_FANCC_27", "HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
               "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11")

for(i in samples_Z) {                  # Head of for-loop
  percentage <- sum(df_all_info$TPM != 0.00 & df_all_info$sample == i, na.rm = TRUE)/ (sum(df_all_info$TPM == 0.00 & df_all_info$sample == i , na.rm = TRUE) + sum(is.na(df_all_info$TPM & df_all_info$sample == i)) + sum(df_all_info$TPM != 0.00 & df_all_info$sample == i, na.rm = TRUE)) * 100 
  percentages_transcribed_Z <- c(percentages_transcribed_Z, percentage)    
}

## Dataframe own SV data including names genes affected and expression data
epd_APH <- expression_data_2B3_23F3_24A7_APH[,c(1, 7, 8, 19, 21)]
epd_E35 <- expression_data_3E5[,c(1, 7, 8, 19, 21)]
epd_APH_3E5 <- full_join(epd_APH, epd_E35)

wanted_rows_expression_data_2 <- epd_APH_3E5[is.na(epd_APH_3E5$Gene_name) != TRUE,]
names(wanted_rows_expression_data_2)[names(wanted_rows_expression_data_2) == "Gene_name"] <- "gene"

data_frame_own$gene <- as.character(data_frame_own$gene)
data_frame_own$gene[data_frame_own$gene== "FAM49B"] <- "CYRIB"

df_all_info_2 <- full_join(data_frame_own, wanted_rows_expression_data_2)

## Calculating percentages transcribed own data
percentages_transcribed_own <- vector()
samples_own <- c("2B3-3E5", "2B3-20F3-24A7 + APH")

for(i in samples_own) {                  
  percentage <- sum(df_all_info_2$TPM != 0.00 & df_all_info_2$clone == i, na.rm = TRUE)/ (sum(df_all_info_2$TPM == 0.00 & df_all_info_2$clone == i , na.rm = TRUE) + sum(is.na(df_all_info_2$TPM & df_all_info_2$clone == i)) + sum(df_all_info_2$TPM != 0.00 & df_all_info_2$clone == i, na.rm = TRUE)) * 100 
  percentages_transcribed_own <- c(percentages_transcribed_own, percentage)    
}

## Dataframe SV located in transcribed genes both own and Zou
percentages_transcribed_both <- c(percentages_transcribed_own, percentages_transcribed_Z)
length(percentages_transcribed_both)

transcribed <- data.frame(sample = c("2B3-3E5", "2B3-20F3-24A7 + APH",
                                     "HAP1_FANCC_27" ,"HAP1_FANCC_27-1", "HAP1_FANCC_27-3", 
                                     "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", "HAP1_FANCC_27-6", 
                                     "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"),
                          transcribed = percentages_transcribed_both,
                          not_transcribed = 100 - percentages_transcribed_both) 

df_transcribed <- rbind(
  data.frame(sample = transcribed$sample, percentage = transcribed$transcribed , type = "Transcribed"),
  data.frame(sample = transcribed$sample, percentage = transcribed$not_transcribed, type ="Non-transcribed")
)

df_transcribed$type <- factor(df_transcribed$type, levels = c("Non-transcribed", "Transcribed"))

# Plot SVs located in transcribed genes
plot_percentage_SV_transcribed <- ggplot(data=df_transcribed, aes(x=sample, y = percentage)) + geom_bar(aes(fill=type), position = 'fill', stat = "identity") +
                                  theme_minimal() +
                                  ggtitle("Percentage of SVs in transcribed region") +
                                  ylab("Percentage") +
                                  labs(fill = "Transcription status") +
                                  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
                                  scale_fill_manual(values=c("firebrick2", "cyan3")) +
                                  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), labels = c("0", "25", "50", "75", "100")) +
                                  scale_x_discrete(limits=c("AHH1-G2", "2B3-3E5", "2B3-20F3-24A7 + APH", " ","HAP1_FANCC_27-1", "HAP1_FANCC_27-3", "HAP1_FANCC_27-4", "HAP1_FANCC_27-5", 
                                                            "HAP1_FANCC_27-6", "HAP1_FANCC_27-8", "HAP1_FANCC_27-11"))

