library(ggplot2)
library(ggsci)
library(ggpubr)
library(rlang)
library(dplyr)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Directory containing unknown-token count CSV files."),
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

plot_results <- function(df_paths, outpref) {
  
  j <- 1
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header = TRUE)
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
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
  df_all['Classifier'][df_all['Sketchlib'] == FALSE] <- "PanBART"
  
  
  df_all$Data <- factor(df_all$Data, levels = c("10%", "100%"))
  df_all$Lineages <- factor(df_all$Lineages, levels = c("All lineages", "Subset lineages", "Held-out lineages"))
  df_all$Classifier <- factor(df_all$Classifier, levels = c("PanBART", "Sketchlib"))
  df_all$Species_long <- gsub("_", ". ", df_all$Species)
  
  # add final row which assigns database name
  df_all$DatasetName <- assign_database_name(df_all)
  df_all_unknown <- subset(df_all, Token == "<unk>")
  
  p <- ggplot(df_all_unknown, aes(x = Count, fill = Species_long)) +
    geom_histogram() +
    facet_grid(Data ~ Lineages, scales = "free_y") +
    theme_light() +
    theme(legend.text  = element_text(face = "italic")) +
    labs(x = "Unknown Count", y = "No. Genomes", fill = "Species") +
    scale_fill_npg()
  p
  ggsave(file=file.path(outpref, "unknown_token_count.png"), plot=p, height = 4, width = 7)
}


indir <- opt$indir
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

df_paths <- list.files(
  path = indir,
  pattern = "unknown_tokens_counts\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
plot_results(df_paths, outdir)
