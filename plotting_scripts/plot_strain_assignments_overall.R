library(ggplot2)
library(ggsci)
library(ggpubr)
library(rlang)
library(dplyr)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Root directory containing species assignment subdirectories."),
  make_option("--outdir", type = "character", help = "Directory for generated plots."),
  make_option("--species", type = "character", default = "S_pneumoniae,E_coli",
              help = "Comma-separated species subdirectory names.")
)
opt <- parse_args(OptionParser(option_list = option_list))

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("_sampled_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("_held_out_removed_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out_only = if (grepl("_held_out_only_", filename)) {
      TRUE
    } else {
      FALSE
    },
    With_unk = if (grepl("with_unk", filename)) {
      TRUE
    } else {
      FALSE
    },
    Sketchlib = if (grepl("sketchlib", filename)) {
      TRUE
    } else {
      FALSE
    },
    Bacformer = if (grepl("bacformer", filename)) {
      TRUE
    } else {
      FALSE
    },
    Bacformer_complete = if (grepl("bacformer_complete", filename)) {
      TRUE
    } else {
      FALSE
    },
    Species = if (grepl("E_coli", filename)) {
      "E_coli"
    } else {
      "S_pneumoniae"
    }
  )
}

assign_database_name <- function(df_all){
  
  df_all_sample <- df_all[, c("Species", "Data", "Lineages")]
  
  df_all_sample$Species_parsed <- case_when(
    df_all_sample$Species == "E_coli"       ~ "Ecoli",
    df_all_sample$Species == "S_pneumoniae" ~ "Pneumo",
    .default = NA_character_
  )
  df_all_sample$Data_parsed <- gsub("%", "", df_all_sample$Data)
  
  df_all_sample$Subset_parsed <- case_when(
    df_all_sample$Lineages == "Subset lineages"       ~ "Sub",
    df_all_sample$Lineages == "All lineages" ~ "All",
    df_all_sample$Lineages == "Held-out lineages" ~ "HO",
    .default = NA_character_
  )
  
  DatasetName <- paste0(df_all_sample$Species_parsed, df_all_sample$Data_parsed, df_all_sample$Subset_parsed)
  DatasetName <- factor(DatasetName, levels <- c("Pneumo10Sub", "Pneumo100Sub", "Pneumo10All", "Pneumo100All", "Pneumo10HO", "Pneumo100HO", "Ecoli10Sub", "Ecoli100Sub", "Ecoli10All", "Ecoli100All", "Ecoli10HO", "Ecoli100HO"))
  DatasetName
}

generate_bal_plots <- function(data_subset, outpref, legend_title, facets)
{
  rows <- "Species_long"
  
  if (facets == "Lineages") {
    p_bal_accuracy <- ggplot(data_subset, aes(x = K, y=Accuracy, colour = Classifier)) +
      geom_smooth(method = "loess", 
                  method.args = list(control = loess.control(
                    surface = "direct"
                  ))) +
      scale_colour_npg() +
      theme_light() +
      facet_wrap(DatasetName ~ ., ncol=4) +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        strip.text.y = element_text(size = 12)
      ) +
      guides(colour=guide_legend(title="Classifier")) +
      labs(x = "k-value", y = "Cluster Assignment Accuracy") +
      scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0)) +
      geom_hline(yintercept = 1, linetype="dashed")
    p_bal_accuracy
    ggsave(file=paste(outpref, "bal_acc_lineages.svg", sep = ""), plot=p_bal_accuracy, height = 5, width = 10)
    ggsave(file=paste(outpref, "bal_acc_lineages.png", sep = ""), plot=p_bal_accuracy, height = 5, width = 10)
  }
  
  if (facets == "Unknown_token") {
    p_bal_accuracy <- ggplot(data_subset, aes(x = K, y=Accuracy, colour = Classifier, linetype = Data)) +
      geom_smooth(method = "loess") +
      scale_colour_npg() +
      theme_light() +
      facet_grid(Species_long ~ Unknown_token) +
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 12),
        strip.text.y = element_text(size = 12, face = "italic")
      ) +
      guides(colour=guide_legend(title="Classifier")) +
      labs(x = "k-value", y = "Cluster Assignment Accuracy") +
      scale_y_continuous(breaks = c(0.0,0.25,0.5,0.75,1.0)) +
      geom_hline(yintercept = 1, linetype="dashed")
    p_bal_accuracy
    ggsave(file=paste(outpref, "bal_acc_unknown_token.svg", sep = ""), plot=p_bal_accuracy, height = 5, width = 8)
    ggsave(file=paste(outpref, "bal_acc_unknown_token.png", sep = ""), plot=p_bal_accuracy, height = 5, width = 8)
  }
 
  
  # p_bal_accuracy
  # p_bal_precision
  # p_bal_recall
  #   p_bal_precision + xlab(""), 
  #   p_bal_recall + xlab(""),
  #   p_bal_accuracy + xlab(""),
  # full_plot
  # full_plot
}

generate_unbal_plots <- function(data_subset, outpref, legend_title, facets)
{
  p_accuracy <- ggplot(data_subset, aes(x = K, y=Accuracy, colour = Data, linetype = Classifier)) +
    geom_smooth(method = "loess") +
    scale_colour_npg() +
    theme_light() +
    facet_grid(cols = vars(!!sym(facets))) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title=legend_title)) +
    labs(x = "k-value", y = "Accuracy") +
    geom_hline(yintercept = 1, linetype="dashed")
  p_accuracy
  ggsave(file=paste(outpref, "_acc.svg", sep = ""), plot=p_accuracy, height = 6, width = 8)
  ggsave(file=paste(outpref, "_acc.png", sep = ""), plot=p_accuracy, height = 6, width = 8)
  
  p_precision <- ggplot(data_subset, aes(x = K, y=Precision, colour = Data, linetype = Classifier)) +
    geom_smooth(method = "loess") +
    scale_colour_npg() +
    theme_light() +
    facet_grid(cols = vars(!!sym(facets))) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title=legend_title)) +
    labs(x = "k-value", y = "Precision") +
    geom_hline(yintercept = 1, linetype="dashed")
  p_precision
  ggsave(file=paste(outpref, "_precision.svg", sep = ""), plot=p_precision, height = 6, width = 8)
  ggsave(file=paste(outpref, "_precision.png", sep = ""), plot=p_precision, height = 6, width = 8)
  
  p_recall <- ggplot(data_subset, aes(x = K, y=Recall, colour = Data, linetype = Classifier)) +
    geom_smooth(method = "loess") +
    scale_colour_npg() +
    theme_light() +
    facet_grid(cols = vars(!!sym(facets))) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    guides(colour=guide_legend(title=legend_title)) +
    labs(x = "k-value", y = "Recall") +
    geom_hline(yintercept = 1, linetype="dashed")
  p_recall
  ggsave(file=paste(outpref, "_recall.svg", sep = ""), plot=p_recall, height = 6, width = 8)
  ggsave(file=paste(outpref, "_recall.png", sep = ""), plot=p_recall, height = 6, width = 8)
}

plot_results <- function(df_paths, outpref) {
  
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
    df$filename <- filename
    
    if (j == 1) {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all['Data'] <- NA
  df_all['Lineages'] <- NA
  df_all['Unknown_token'] <- NA
  df_all['Classifier'] <- NA
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10%"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100%"
  
  df_all['Lineages'][df_all['Held_out'] == TRUE] <- "Subset lineages"
  df_all['Lineages'][df_all['Held_out'] == FALSE] <- "All lineages"
  df_all['Lineages'][df_all['Held_out_only'] == TRUE] <- "Held-out lineages"
  
  df_all['Unknown_token'][df_all['With_unk'] == TRUE] <- "With unknown"
  df_all['Unknown_token'][df_all['With_unk'] == FALSE] <- "Without unknown"
  
  df_all['Classifier'][df_all['Sketchlib'] == TRUE] <- "Sketchlib"
  df_all['Classifier'][df_all['Sketchlib'] == FALSE & df_all['Bacformer'] == TRUE & df_all['Bacformer_complete'] == FALSE] <- "Bacformer-L (MAG)"
  df_all['Classifier'][df_all['Sketchlib'] == FALSE & df_all['Bacformer'] == TRUE & df_all['Bacformer_complete'] == TRUE] <- "Bacformer-L (Complete)"
  df_all['Classifier'][df_all['Sketchlib'] == FALSE & df_all['Bacformer'] == FALSE & df_all['Bacformer_complete'] == FALSE] <- "PanBART"
  
  df_all$Data <- factor(df_all$Data, levels = c("10%", "100%"))
  df_all$Lineages <- factor(df_all$Lineages, levels = c("All lineages", "Subset lineages", "Held-out lineages"))
  df_all$Classifier <- factor(df_all$Classifier, levels = c("PanBART", "Sketchlib", "Bacformer-L (MAG)", "Bacformer-L (Complete)"))
  df_all$Species_long <- gsub("_", ". ", df_all$Species)
  
  # add final row which assigns database name
  df_all$DatasetName <- assign_database_name(df_all)
  
  legend_title <- "Training\ndataset"
  facets <- "Lineages"
  data_subset = subset(df_all, (Label != "overall") & (Lineages != "Held-out lineages") & ((Held_out_only != TRUE & With_unk == TRUE) | (With_unk != TRUE & (Sketchlib == TRUE | Bacformer == TRUE))))
  generate_bal_plots(data_subset, outpref, legend_title, facets)
  
  facets <- "Unknown_token"
  data_subset = subset(df_all, Label != "overall" & (Held_out_only == TRUE | Sketchlib == TRUE))
  generate_bal_plots(data_subset, paste0(outpref, "held_out_only"), legend_title, facets)
}

indir <- opt$indir
outdir <- opt$outdir
species_list <- strsplit(opt$species, ",", fixed = TRUE)[[1]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

{
  df_paths <- unlist(lapply(species_list, function(sp) {
    list.files(
      path = file.path(indir, sp),
      pattern = "per_label_accuracy\\.tsv$",
      recursive = TRUE,
      full.names = TRUE
    )
  }))
  
  outpref <- paste0(outdir, .Platform$file.sep)
  plot_results(df_paths, outpref)
}
