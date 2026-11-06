library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(tidyverse)
library(optparse)

option_list <- list(
  make_option("--root", type = "character", help = "Root input directory containing species subdirectories."),
  make_option("--outdir", type = "character", help = "Directory for saved plots."),
  make_option("--species", type = "character", help = "Comma-separated species subdirectory names."),
  make_option("--species-names", type = "character", help = "Comma-separated species display names.")
)
opt <- parse_args(OptionParser(option_list = option_list))

root <- opt$root
outdir <- opt$outdir
species_list <- strsplit(opt$species, ",", fixed = TRUE)[[1]]
species_name_list <- strsplit(opt$`species-names`, ",", fixed = TRUE)[[1]]
if (length(species_list) != length(species_name_list)) {
  stop("--species and --species-names must contain the same number of values.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

j <- 1
for (j in 1:length(species_list))
{
  species <- species_list[[j]]
  species_name <- species_name_list[[j]]
  indir <- file.path(root, species)
  
  # read in gene frequencies and generate dictionary
  absolute_gene_IDs <- TRUE
  
  gene_freq_file <- file.path(indir, "train_gene_counts_abs.txt")
  gene.freq.df <- read.csv(gene_freq_file, header = 1, sep = "\t")
  pseudolikelihood_file <- file.path(indir, "all_genomes_train_with_unk_per_gene.txt")
  pseudolikelihood_random_file <- file.path(indir, "all_genomes_train_with_unk_random_per_gene.txt")
  pseudolikelihood.df <- read.csv(pseudolikelihood_file, header = 1, sep = "\t")
  pseudolikelihood.random.df <- read.csv(pseudolikelihood_random_file, header = 1, sep = "\t")
  gene.freq.df$Gene_ID <- as.numeric(gene.freq.df$Gene_ID)
  pseudolikelihood.df$Gene_ID <- as.numeric(pseudolikelihood.df$Gene_ID)
  pseudolikelihood.random.df$Gene_ID <- as.numeric(pseudolikelihood.random.df$Gene_ID)
  
  if (absolute_gene_IDs) {
    Total_genomes_tmp <- gene.freq.df$Total_genomes[1]
    Total_genes_tmp <- gene.freq.df$Total_genes[1]
    gene.freq.df$Gene_ID <- abs(gene.freq.df$Gene_ID)
    gene.freq.df_tmp <- gene.freq.df
    gene.freq.df <- gene.freq.df %>% 
      group_by(Gene_ID) %>% 
      summarize(Genome_count = sum(Genome_count), 
                Total_count = sum(Total_count))
    
    gene.freq.df$Genome_freq <- gene.freq.df$Genome_count / Total_genomes_tmp
    gene.freq.df$Total_freq <- gene.freq.df$Total_count / Total_genes_tmp
    
    
    pseudolikelihood.df$Gene_ID <- abs(pseudolikelihood.df$Gene_ID)
    pseudolikelihood.random.df$Gene_ID <- abs(pseudolikelihood.random.df$Gene_ID)
  }
  summary.df.pseudolikelihood <- pseudolikelihood.df %>% 
    group_by(Gene_ID) %>% 
    summarize(Avg_log_pseudolikelihood = mean(log_pseudolikelihood), Std_log_pseudolikelihood = sd(log_pseudolikelihood))
  
  summary.df.pseudolikelihood <- summary.df.pseudolikelihood[!is.na(summary.df.pseudolikelihood$Gene_ID),]
  summary.df.pseudolikelihood$Type <- "Non-random"
  
  summary.df.pseudolikelihood.random <- pseudolikelihood.random.df %>% 
    group_by(Gene_ID) %>% 
    summarize(Avg_log_pseudolikelihood = mean(log_pseudolikelihood), Std_log_pseudolikelihood = sd(log_pseudolikelihood))
  summary.df.pseudolikelihood.random$Type <- "Random"
  
  summary.df.pseudolikelihood.random <- summary.df.pseudolikelihood.random[!is.na(summary.df.pseudolikelihood.random$Gene_ID),]
  
  combined.df.nonrandom <- merge(summary.df.pseudolikelihood, gene.freq.df, by = c("Gene_ID"))
  combined.df.random <- merge(summary.df.pseudolikelihood.random, gene.freq.df, by = c("Gene_ID"))
  
  all.merged <- rbind(combined.df.nonrandom, combined.df.random)
  all.merged$Species <- species_name
  
  if (j == 1){
    all.merged.total <- all.merged
  } else {
    all.merged.total <- rbind(all.merged.total, all.merged)
  }
}

# plot average
p <- ggplot(all.merged.total, aes(x = Genome_freq, y=Avg_log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  facet_grid(Species ~ Type) +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file=file.path(outdir, "log_pseudolikelihood_vs_gene_freq_training_scatter.png"), plot=p, height = 6, width = 8)

# plot average
p <- ggplot(all.merged.total, aes(x = Genome_freq, y=Avg_log_pseudolikelihood,)) +
  geom_point(alpha=0.2) +
  facet_grid(Species ~ Type) +
  geom_smooth(method="loess", colour="blue", fill = "red", linewidth = 0.5) +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood")
p

ggsave(file=file.path(outdir, "log_pseudolikelihood_vs_gene_freq_training_loess.png"), plot=p, height = 6, width = 8)


breaks <- seq(0, 1, by = 0.1)
p <- ggplot(all.merged.total, aes(x = cut(Genome_freq, breaks = breaks), y=Avg_log_pseudolikelihood,)) +
  geom_boxplot() +
  facet_grid(Species ~ Type) +
  theme_light() +
  labs(x = "Gene Frequency", y = "Gene log pseudolikelihood") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(file=file.path(outdir, "log_pseudolikelihood_vs_gene_freq_training_boxplot.png"), plot=p, height = 6, width = 11)
