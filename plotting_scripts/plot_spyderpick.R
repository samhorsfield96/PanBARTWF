library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(ggrepel)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Directory containing Spyderpick outlier files."),
  make_option("--outdir", type = "character", help = "Directory for generated plots.")
)
opt <- parse_args(OptionParser(option_list = option_list))

indir <- opt$indir
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
direct_extreme <- list.files(indir, pattern = "_direct_extreme_outliers\\.tsv$", full.names = TRUE)
direct_normal <- list.files(indir, pattern = "_direct_outliers\\.tsv$", full.names = TRUE)
indirect_extreme <- list.files(indir, pattern = "_indirect_extreme_outliers\\.tsv$", full.names = TRUE)
indirect_normal <- list.files(indir, pattern = "_indirect_outliers\\.tsv$", full.names = TRUE)

direct_extreme_df <- read.csv(direct_extreme, header = TRUE, sep = "\t")
direct_extreme_df$Direct <- TRUE
direct_extreme_df$Extreme <- TRUE
direct_extreme_df$Label <- "Direct extreme outlier"

direct_normal_df <- read.csv(direct_normal, header = TRUE, sep = "\t")
direct_normal_df$Direct <- TRUE
direct_normal_df$Extreme <- FALSE
direct_normal_df$Label <- "Direct outlier"

indirect_extreme_df <- read.csv(indirect_extreme, header = TRUE, sep = "\t")
indirect_extreme_df$Direct <- FALSE
indirect_extreme_df$Extreme <- TRUE
indirect_extreme_df$Label <- "Indirect extreme outlier"

indirect_normal_df <- read.csv(indirect_normal, header = TRUE, sep = "\t")
indirect_normal_df$Direct <- FALSE
indirect_normal_df$Extreme <- FALSE
indirect_normal_df$Label <- "Indirect outlier"

total_df <- rbind(direct_extreme_df, direct_normal_df, indirect_extreme_df, indirect_normal_df)

total_df$Direct[total_df$Direct == TRUE] <- "Direct"
total_df$Direct[total_df$Direct == FALSE] <- "Indirect"
total_df$Extreme[total_df$Extreme == TRUE] <- "Extreme outlier"
total_df$Extreme[total_df$Extreme == FALSE] <- "Outlier"

hypothetical.label <- "?"
unnamed.label <- "?"
NA.label <- "?"
total_df$Gene[total_df$Gene == "" & total_df$Product == "hypothetical protein"] <- hypothetical.label
total_df$Gene.1[total_df$Gene.1 == "" & total_df$Product.1 == "hypothetical protein"] <- hypothetical.label
total_df$Gene[total_df$Gene == "" & total_df$Product != "hypothetical protein"] <- unnamed.label
total_df$Gene.1[total_df$Gene.1 == "" & total_df$Product.1 != "hypothetical protein"] <- unnamed.label
total_df$Gene[total_df$Gene == "NA" | is.na(total_df$Gene)] <- NA.label
total_df$Gene.1[total_df$Gene.1 == "NA" | is.na(total_df$Gene.1)] <- NA.label
total_df$gene.pairs <- paste0(total_df$Gene, "~", total_df$Gene.1)

to.remove <- c(hypothetical.label, unnamed.label, "NA")
total_df$gene.pairs[total_df$Gene %in% to.remove | total_df$Gene.1 %in% to.remove] <- ""

colours <- c("Outlier" = "#696969", 
             "Extreme outlier" = "#E64B35FF")

extreme.outlier.threshold <- 0.6255980000000001
outlier.threshold <- 0.39962337500000006

pos <- position_jitter(width = 0.2)

p.all <- ggplot(total_df, aes(x = Direct, y = MI, colour = Extreme)) +
  geom_hline(yintercept = extreme.outlier.threshold, colour = "black", linetype = "dashed") +
  geom_jitter(size=2, alpha=0.9, position = position_jitter(seed = 1)) +
  guides(colour=guide_legend(title="Outlier type")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=14,face="bold"), strip.text.y = element_text(size = 14), 
        legend.title=element_text(size=12,face="bold"), legend.text=element_text(size=12)) +
  scale_color_manual(values = colours) + 
  geom_text_repel(
    data = total_df,
    aes(label = ifelse(Extreme %in% c("Extreme outlier"), gene.pairs, "")),
    color = "black",
    segment.color = "gray20",
    max.overlaps = 100,
    position = position_jitter(seed = 1),
    min.segment.length = unit(0, 'lines'), 
  ) +
  labs(x = "Connection Type", y = "Mutual Information (MI)") 
p.all

ggsave(file=file.path(outdir, "E_coli_spyderpick_direct_only_MI_None.png"), plot=p.all, height = 6, width = 8)
ggsave(file=file.path(outdir, "E_coli_spyderpick_direct_only_MI_None.svg"), plot=p.all, height = 6, width = 8)
