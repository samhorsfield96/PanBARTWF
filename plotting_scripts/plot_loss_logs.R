library(ggplot2)
library(ggsci)
library(ggpubr)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Directory containing Pneumo/ and Ecoli/ log folders."),
  make_option("--outdir", type = "character", help = "Directory for generated plots.")
)
opt <- parse_args(OptionParser(option_list = option_list))

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("_sampled_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("_held_out_", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

plot_results <- function(df_paths, outpref) {
  df_all <- data.frame(Type = c(), Epoch = c(), Loss = c(), Perplexity = c(), Learning_rate = c(), Accuracy = c(), Precision = c(), Recall = c(), F1 = c(), Sampled = c(), Held_out = c())
  
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "", header = TRUE)
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    df_all <- rbind(df_all, df)
  }
  
  df_all['Data'] <- NA
  df_all['Lineages'] <- NA
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10%"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100%"
  
  df_all['Lineages'][df_all['Held_out'] == TRUE] <- "Subset lineages"
  df_all['Lineages'][df_all['Held_out'] == FALSE] <- "All lineages"
  
  #df_all['Data'] <- factor(df_all['Data'], levels = c("10% data", "100% data"))
  #df_all['Lineages'] <- factor(df_all['Lineages'], levels = c("Subset lineages", "All lineages"))
  
  data_subset = subset(df_all, Type != "Test")
  p_loss <- ggplot(data_subset, aes(x = Epoch, y=Loss, colour = Data)) +
    geom_line(linewidth=1.5) +
    scale_colour_npg() +
    theme_light() +
    facet_grid(rows = vars(Type), cols = vars(Lineages)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title="Training\ndataset")) +
    scale_y_log10() +
    labs(x = "Epoch", y = "Cross Entropy Loss")
  p_loss
  ggsave(file=paste(outpref, "_loss.svg", sep = ""), plot=p_loss, height = 6, width = 8)
  ggsave(file=paste(outpref, "_loss.png", sep = ""), plot=p_loss, height = 6, width = 8)
  
  data_subset = subset(df_all, Type == "Validation")
  p_precision <- ggplot(data_subset, aes(x = Epoch, y=Precision, colour = Data)) +
    geom_line(linewidth=1.5) +
    scale_colour_npg() +
    theme_light() +
    facet_grid(rows = vars(Type), cols = vars(Lineages)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title="Training\ndataset")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Epoch", y = "Precision")
  p_precision
  ggsave(file=paste(outpref, "_precision.svg", sep = ""), plot=p_precision, height = 6, width = 8)
  ggsave(file=paste(outpref, "_precision.png", sep = ""), plot=p_precision, height = 6, width = 8)
  
  p_recall <- ggplot(data_subset, aes(x = Epoch, y=Recall, colour = Data)) +
    geom_line(linewidth=1.5) +
    scale_colour_npg() +
    theme_light() +
    facet_grid(rows = vars(Type), cols = vars(Lineages)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title="Training\ndataset")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Epoch", y = "Recall")
  p_recall
  ggsave(file=paste(outpref, "_recall.svg", sep = ""), plot=p_recall, height = 6, width = 8)
  ggsave(file=paste(outpref, "_recall.png", sep = ""), plot=p_recall, height = 6, width = 8)
  
  p_f1 <- ggplot(data_subset, aes(x = Epoch, y=F1, colour = Data)) +
    geom_line(linewidth=1.5) +
    scale_colour_npg() +
    theme_light() +
    facet_grid(rows = vars(Type), cols = vars(Lineages)) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title="Training\ndataset")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Epoch", y = "F1 score")
  p_f1
  ggsave(file=paste(outpref, "_F1.svg", sep = ""), plot=p_f1, height = 6, width = 8)
  ggsave(file=paste(outpref, "_F1.png", sep = ""), plot=p_f1, height = 6, width = 8)
  
  full_plot <- ggarrange(
    p_loss + xlab(""), 
    p_precision + xlab(""),
    p_recall + xlab(""),
    p_f1 + xlab(""), 
    nrow = 4,
    ncol = 1,
    align = "hv",
    common.legend = TRUE,
    legend = "right",
    labels = "AUTO",
    heights = c(2, 1, 1, 1))
  full_plot
  
  full_plot <- annotate_figure(full_plot,
                                bottom = text_grob("Epoch", color = "black", face = "bold", size = 18))
  full_plot
  ggsave(file=paste(outpref, "_full.svg", sep = ""), plot=full_plot, height = 14, width = 8)
  ggsave(file=paste(outpref, "_full.png", sep = ""), plot=full_plot, height = 14, width = 8)
}

outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species_dirs <- c(
  S_pneumoniae = file.path(opt$indir, "Pneumo"),
  E_coli = file.path(opt$indir, "Ecoli")
)
for (species in names(species_dirs)) {
  indir <- species_dirs[[species]]
  df_paths <- list.files(indir, pattern = "_logs\\.txt$", full.names = TRUE)
  if (length(df_paths) == 0) {
    stop("No log files found in ", indir)
  }
  plot_results(df_paths, file.path(outdir, species))
}

