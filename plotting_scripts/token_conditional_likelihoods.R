library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyverse)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Directory containing token pseudolikelihood files."),
  make_option("--outdir", type = "character", help = "Directory for generated plots.")
)
opt <- parse_args(OptionParser(option_list = option_list))

colours <- c("ST131 B/C (+ve)" = "#696969", 
             "ST131 B/C (-ve)" = "#E64B35FF",
             "ST131 A (+ve)" = "#4DBBD5FF",
             "ST131 A (-ve)" = "#4DBBD5FF",
             "ST95 (+ve)" = "#00A087FF",
             "ST95 (-ve)" = "#00A087FF",
             "ST73 (+ve)" = "#3C5488FF",
             "ST73 (-ve)" = "#3C5488FF")

parse_filename <- function(filename) {
  list(
    bla_pos = if (grepl("\\+ve", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST131_BC = if (grepl("ST131B_C", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST131_A = if (grepl("ST131A", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST73 = if (grepl("CC73", filename)) {
      TRUE
    } else {
      FALSE
    },
    ST95 = if (grepl("CC95", filename)) {
      TRUE
    } else {
      FALSE
    },
    train = if (grepl("train", filename)) {
      TRUE
    } else {
      FALSE
    },
    test = if (grepl("test", filename)) {
      TRUE
    } else {
      FALSE
    },
    val = if (grepl("val", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

indir <- opt$indir
outdir <- paste0(opt$outdir, .Platform$file.sep)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(indir, pattern = "locus_pseudolikelihood\\.txt$", full.names = TRUE)

# non pres
labels <- c("ST131 B/C (-ve)", "ST131 B/C (+ve)", "ST131 A (-ve)", "ST131 A (+ve)", "ST73 (-ve)", "ST73 (+ve)", "ST95 (-ve)", "ST95 (+ve)")
#pres



i <- 15
for (i in 1:length(files))
{
  filename <- files[[i]]
  df.raw <- read.csv(filename, sep = "\t")
  df.rows <- df.raw[, c(1, 2)]
  df <- df.raw[, -c(1, 2)]
  df[df == 0] <- NA
  
  parsed_values <- parse_filename(filename)
  
  lineage <- NA
  lineage <- if (parsed_values$ST131_BC & is.na(lineage)) "ST131 B/C" else lineage
  lineage <- if (parsed_values$ST131_A & is.na(lineage)) "ST131 A" else lineage
  lineage <- if (parsed_values$ST73 & is.na(lineage)) "ST73" else lineage
  lineage <- if (parsed_values$ST95 & is.na(lineage)) "ST95" else lineage
  
  bla_pos <- parsed_values$bla_pos
  bla_pos_str <- if (parsed_values$bla_pos) "+ve" else "-ve"
  
  data_type <- NA
  data_type <- if (parsed_values$train & is.na(data_type)) "Training" else data_type
  data_type <- if (parsed_values$val & is.na(data_type)) "Validation" else data_type
  data_type <- if (parsed_values$test & is.na(data_type)) "Testing" else data_type
  
  quantiles.val <- t(apply(df,1,quantile,c(0.25,0.50,0.75), na.rm= TRUE))
  avg.val <- apply(df,1,mean, na.rm= TRUE)
  max.val <- apply(df,1,max, na.rm= TRUE)
  min.val <- apply(df,1,min, na.rm= TRUE)
  
  label <- paste0(lineage, " (", bla_pos_str, ")")
  
  to.append <- data.frame(Lineage = lineage, bla_pos = bla_pos, label = label, data_type = data_type, Dataset = filename, GeneID = df.rows$Gene_ID, GenomeID = df.rows$Genome_Index, Quantile25 = quantiles.val[,1], Quantile50 = quantiles.val[,2], Quantile75 = quantiles.val[,3], Min = min.val, Max = max.val, Mean = avg.val)
  
  if (i == 1){
    final.df <- to.append
  } else
  {
    final.df <- rbind(final.df, to.append)
  }
}

# ordering
ordering <- c("ST131 B/C (+ve)", "ST131 B/C (-ve)", "ST131 A (+ve)", "ST131 A (-ve)", "ST95 (+ve)", "ST95 (-ve)", "ST73 (+ve)", "ST73 (-ve)")
final.df$label <- factor(final.df$label, levels = ordering)
final.df$data_type <- factor(final.df$data_type, levels = c("Training", "Testing"))

# negative genomes
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST131 A (+ve)"),
                             c("ST131 B/C (+ve)", "ST131 A (-ve)"), c("ST131 B/C (+ve)", "ST95 (+ve)"),
                             c("ST131 B/C (+ve)", "ST95 (-ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"))
p.all <- ggplot(final.df, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(data_type~ .) +
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = "log pseudolikelihood", x = "Lineage")
p.all

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_all.png"), plot=p.all, height = 11, width = 10)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_all.svg"), plot=p.all, height = 11, width = 10)

# positive genomes
subset.df.pos <- subset(final.df, bla_pos == TRUE)
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 A (+ve)"), c("ST131 B/C (+ve)", "ST95 (+ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"))
p.pos <- ggplot(subset.df.pos, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(data_type~ .) +
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = "log pseudolikelihood", x = "Lineage")
p.pos

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_pos_only.png"), plot=p.pos, height = 11, width = 10)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_pos_only.svg"), plot=p.pos, height = 11, width = 10)

# negative genomes
subset.df.neg <- subset(final.df, bla_pos == FALSE | label == "ST131 B/C (+ve)")
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST131 A (-ve)"), c("ST131 B/C (+ve)", "ST95 (-ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"))
p.neg <- ggplot(subset.df.neg, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(data_type~ .) +
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = "log pseudolikelihood", x = "Lineage")
p.neg

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_neg_only.png"), plot=p.neg, height = 11, width = 10)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_neg_only.svg"), plot=p.neg, height = 11, width = 10)

# for presentation
colours_pres <- c("ST131 B/C (+ve)" = "#800000", 
             "ST131 B/C (-ve)" = "#FA8072",
             "ST73 (+ve)" = "#00008B",
             "ST73 (-ve)" = "#4169E1")

subset.df.pres <- subset(final.df, Lineage == "ST131 B/C" | Lineage == "ST73")
subset.df.pres <- subset(subset.df.pres, data_type == "Testing")
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"))
p.pres <- ggplot(subset.df.pres, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours_pres) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = "log pseudolikelihood", x = "Lineage")
p.pres

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_pres.png"), plot=p.pres, height = 7, width = 11)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_pres.svg"), plot=p.pres, height = 7, width = 11)

# for presentation, all pairwise comps
colours_pres <- c("ST131 B/C (+ve)" = "#800000", 
                  "ST131 B/C (-ve)" = "#FA8072",
                  "ST73 (+ve)" = "#00008B",
                  "ST73 (-ve)" = "#4169E1")

subset.df.pres <- subset(final.df, Lineage == "ST131 B/C" | Lineage == "ST73")
subset.df.pres <- subset(subset.df.pres, data_type == "Testing")
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"), c("ST73 (+ve)", "ST131 B/C (-ve)"), c("ST73 (+ve)", "ST73 (-ve)"))
p.pres <- ggplot(subset.df.pres, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours_pres) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = "log pseudolikelihood", x = "Lineage")
p.pres

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_pres_all_comps.png"), plot=p.pres, height = 7, width = 11)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_pres_all_comps.png"), plot=p.pres, height = 7, width = 11)


subset.df <- subset(final.df, label == "ST73 (+ve)")
ST73.df.ordered <- subset.df %>%
  arrange(desc(data_type), desc(Max))

subset.df <- subset(final.df, label == "ST95 (+ve)")
ST95.df.ordered <- subset.df %>%
  arrange(desc(data_type), desc(Max))


# testing only
subset.df.pos <- subset(final.df, bla_pos == TRUE & data_type == "Testing")
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 A (+ve)"), c("ST131 B/C (+ve)", "ST95 (+ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"))
p.pos.test <- ggplot(subset.df.pos, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = NULL, x = NULL) +
  scale_y_continuous(limits=c(2,-15))
p.pos.test

subset.df.neg <- subset(final.df, data_type == "Testing" & (bla_pos == FALSE | label == "ST131 B/C (+ve)"))
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST131 A (-ve)"), c("ST131 B/C (+ve)", "ST95 (-ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"))
p.neg.test <- ggplot(subset.df.neg, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif") +
  labs(y = NULL, x = NULL) +
  scale_y_continuous(limits=c(2,-15))
p.neg.test

fig <- ggarrange(p.pos.test, p.neg.test,
                 labels = "AUTO",
                 ncol = 2, nrow = 1)

fig <- annotate_figure(fig,
                bottom = text_grob("Lineage", size = 20, hjust = 0, face = "bold"),
                left   = text_grob("log pseudolikelihood", size = 20, rot = 90, hjust = 0, face = "bold"))

fig

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_testing_only_all.svg"), plot=fig, height = 7, width = 11)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_testing_only_all.png"), plot=fig, height = 7, width = 11)

# training only
subset.df.pos <- subset(final.df, bla_pos == TRUE & data_type == "Training")
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 A (+ve)"), c("ST131 B/C (+ve)", "ST95 (+ve)"), c("ST131 B/C (+ve)", "ST73 (+ve)"))
p.pos.test <- ggplot(subset.df.pos, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = NULL, x = NULL) +
  scale_y_continuous(limits=c(2,-15))
p.pos.test

subset.df.neg <- subset(final.df, data_type == "Training" & (bla_pos == FALSE | label == "ST131 B/C (+ve)"))
pairwise_comparisons <- list(c("ST131 B/C (+ve)", "ST131 B/C (-ve)"), c("ST131 B/C (+ve)", "ST131 A (-ve)"), c("ST131 B/C (+ve)", "ST95 (-ve)"), c("ST131 B/C (+ve)", "ST73 (-ve)"))
p.neg.test <- ggplot(subset.df.neg, aes(x = label, y = Max, colour = label)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 16, angle = 45, hjust = 1), axis.text.y = element_text(size = 16), axis.title=element_text(size=20,face="bold"), strip.text.y = element_text(size = 16), 
        legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_manual(values = colours) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, label = "p.signif", p.adjust.method = "bonferroni") +
  labs(y = NULL, x = NULL) +
  scale_y_continuous(limits=c(2,-15))
p.neg.test

fig <- ggarrange(p.pos.test, p.neg.test,
                 labels = "AUTO",
                 ncol = 2, nrow = 1)

fig <- annotate_figure(fig,
                       bottom = text_grob("Lineage", size = 20, hjust = 0, face = "bold"),
                       left   = text_grob("log pseudolikelihood", size = 20, rot = 90, hjust = 0, face = "bold"))

fig

ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_training_only_all.svg"), plot=fig, height = 7, width = 11)
ggsave(file=paste0(outdir, "Ecoli_blaCTX-M_max_value_training_only_all.png"), plot=fig, height = 7, width = 11)
