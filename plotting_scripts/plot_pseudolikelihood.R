library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(optparse)

option_list <- list(
  make_option("--infile", type="character", help="Path to pseudolikelihood results .xlsx file."),
  make_option("--sheet", type="character", default="All_Pneumo", help="Sheet name to read from the .xlsx file. Default = 'All_Pneumo'"),
  make_option("--outpref", type="character", default="output", help="Output prefix for saved plots. Default = 'output'")
)
opt <- parse_args(OptionParser(option_list=option_list))

# plot trained genomes as random
data <- read.xlsx(opt$infile, opt$sheet)
data$log_pseudolikelihood <- as.numeric(data$log_pseudolikelihood)
data$Type <- factor(data$Type, levels=c("Random", "Diff. Species", "Trained", "Similar", "Novel"))
data$Inference <- factor(data$Inference, levels=c("Encoder", "Encoder-decoder"))
data$Random <- factor(data$Random, levels=c("TRUE", "FALSE"))
data$pseudolikelihood <- exp(data$log_pseudolikelihood)

pairwise_comparisons <- list(c("Trained", "Similar"), c("Trained", "Novel"), c("Similar", "Novel"), 
                             c("Similar", "Random"), c("Novel", "Random"), c("Trained", "Random"))
data_subset <- subset(data, Inference == "Encoder-decoder")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=2, alpha=0.9) +
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_color_npg() + 
  stat_summary(fun=mean, colour="black", geom="text", show_guide = FALSE, vjust=-1, size = 5, aes( label=round(after_stat(y), digits=2))) + 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=2, show_guide = FALSE) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p

ggsave(file=paste0(opt$outpref, "_comparison_encoder_decoder_all_plot.svg"), plot=p, height = 6, width = 10)
ggsave(file=paste0(opt$outpref, "_comparison_encoder_decoder_all_plot.png"), plot=p, height = 6, width = 10)

# for presentation
colour_vec <- c("Trained" = "#CD4309", "Similar" = "#2f47ed",  "Novel" = "#fcb43a", "Random" = "#0a895e", "Diff. Species" = "#8e6f48")
pairwise_comparisons <- list(c("Trained", "Similar"), c("Trained", "Novel"), c("Similar", "Novel"))
data_subset <- subset(data, Inference == "Encoder-decoder" & Type != "Random" & Type != "Diff. Species")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_colour_manual(values = colour_vec) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, label = "p.signif", comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p
ggsave(file=paste0(opt$outpref, "_comparison_presentation_trained.png"), plot=p, height = 6, width = 8)

pairwise_comparisons <- list(c("Random", "Similar"), 
                             c("Novel", "Random"), c("Novel", "Similar"))
data_subset <- subset(data, Inference == "Encoder-decoder" & Type != "Trained" & Type != "Diff. Species")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_colour_manual(values = colour_vec) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, label = "p.signif", comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p
ggsave(file=paste0(opt$outpref, "_comparison_presentation_random.png"), plot=p, height = 6, width = 8)

# all comparisons
pairwise_comparisons <- list(c("Diff. Species", "Random"), c("Random", "Similar"), 
                             c("Novel", "Random"))
data_subset <- subset(data, Inference == "Encoder-decoder" & Type != "Trained")
p <- ggplot(data_subset, aes(x = Type, y = log_pseudolikelihood, colour = Type)) +
  geom_boxplot(outlier.shape = NA) + 
  guides(colour=guide_legend(title="Alignment")) +
  theme_light() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), strip.text.x = element_text(size = 11), legend.title=element_text(size=14,face="bold"), legend.text=element_text(size=12), legend.position = "none") +
  scale_colour_manual(values = colour_vec) + 
  stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, label = "p.signif", comparisons = pairwise_comparisons) +
  labs(y = "log pseudolikelihood", x = "Prompted Genome")
p
ggsave(file=paste0(opt$outpref, "_comparison_diff_species_random.png"), plot=p, height = 6, width = 8)
