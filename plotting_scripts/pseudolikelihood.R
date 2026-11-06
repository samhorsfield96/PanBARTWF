library(ggplot2)
library(ggsci)
library(ggpubr)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Root directory containing species pseudolikelihood subdirectories."),
  make_option("--outdir", type = "character", help = "Directory for generated plots."),
  make_option("--species", type = "character", default = "S_pneumoniae", help = "Species subdirectory name.")
)
opt <- parse_args(OptionParser(option_list = option_list))

parse_filename <- function(filename) {
  list(
    Sampled = if (grepl("_sampled_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Train = if (grepl("_train_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Test = if (grepl("_test_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("_held_out_w", filename)) {
      TRUE
    } else {
      FALSE
    },
    Random = if (grepl("_random_", filename)) {
      TRUE
    } else {
      FALSE
    },
    With_unk = if (grepl("_with_unk_", filename)) {
      TRUE
    } else {
      FALSE
    },
    recomb_2 = if (grepl("_2_transfers_", filename)) {
      TRUE
    } else {
      FALSE
    },
    recomb_10 = if (grepl("_10_transfers_", filename)) {
      TRUE
    } else {
      FALSE
    },
    recomb_20 = if (grepl("_20_transfers_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Smitis = if (grepl("_mitis_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Species = "S_pneumoniae",
    Subset = "Subset Lineages"
  )
}

assign_database_name <- function(df_all){
  
  df_all_sample <- df_all[, c("Species", "Data", "Subset")]
  
  df_all_sample$Species_parsed <- case_when(
    df_all_sample$Species == "E_coli"       ~ "Ecoli",
    df_all_sample$Species == "S_pneumoniae" ~ "Pneumo",
    .default = NA_character_
  )
  df_all_sample$Data_parsed <- gsub("% Data", "", df_all_sample$Data)
  df_all_sample$Subset_parsed <- substr(df_all_sample$Subset, 1, 3)
  
  DatasetName <- paste0(df_all_sample$Species_parsed, df_all_sample$Data_parsed, df_all_sample$Subset_parsed)
  DatasetName <- factor(DatasetName, levels <- c("Pneumo10Sub", "Pneumo100Sub", "Pneumo10All", "Pneumo100All", "Ecoli10Sub", "Ecoli100Sub", "Ecoli10All", "Ecoli100All"))
  DatasetName
}

plot_results <- function(df_paths, outpref)
{
 
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
    
    if (j == 1)
    {
      df_all <- df
    } else {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all['Unknown_token'] <- NA
  df_all['Query'] <- NA
  df_all['Data'] <- NA
  
  df_all['Unknown_token'][df_all['With_unk'] == TRUE] <- "With unknown"
  df_all['Unknown_token'][df_all['With_unk'] == FALSE] <- "Without unknown"
  
  df_all['Query'][df_all['Smitis'] == TRUE] <- "Diff. species"
  df_all['Query'][df_all['Train'] == TRUE] <- "Training"
  df_all['Query'][df_all['Test'] == TRUE] <- "Testing"
  df_all['Query'][df_all['Held_out'] == TRUE] <- "Novel"
  df_all['Query'][df_all['Random'] == TRUE] <- "Random"
  df_all['Query'][df_all['recomb_2'] == TRUE] <- "2 transfers"
  df_all['Query'][df_all['recomb_10'] == TRUE] <- "10 transfers"
  df_all['Query'][df_all['recomb_20'] == TRUE] <- "20 transfers"
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10% Data"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100% Data"
  df_all$DatasetName <- assign_database_name(df_all)
  
  df_all$Query <- factor(df_all$Query, levels = c("Training", "Testing", "2 transfers", "10 transfers", "20 transfers", "Novel", "Diff. species", "Random"))
  
  pairwise_comparisons <- list(c("Training", "Testing"), c("Testing", "Novel"), c("Novel", "Diff. species"))
  data_subset <- subset(df_all, Unknown_token == "With unknown" & (Query != "Random" & Query != "2 transfers" & Query != "10 transfers" & Query != "20 transfers"))
  p <- ggplot(data_subset, aes(x = Query, y=sum, colour = Query)) +
    geom_boxplot() +
    scale_colour_npg() +
    facet_grid(DatasetName ~ .) +
    theme_light() +
    stat_compare_means(
      comparisons = pairwise_comparisons,
      method = "wilcox.test",
      paired = FALSE,
      label = "p.signif",
      p.adjust.method = "bonferroni",
      size = 5.5,
      label.x = 1.0,
      symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
        symbols = c("****", "***", "**", "*", "ns")
      )
    ) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
    guides(colour=guide_legend(title="Dataset size")) +
    labs(x = "Query type", y = "Log Pseudolikelihood") +
    geom_hline(yintercept = 0, linetype="dashed")
  p
  ggsave(file=paste(outpref, "_pseudolikelihood_comparison_w_unknown_no_random.svg", sep = ""), plot=p, height = 10, width = 8)
  ggsave(file=paste(outpref, "_pseudolikelihood_comparison_w_unknown_no_random.png", sep = ""), plot=p, height = 10, width = 8)
  
  pairwise_comparisons <- list(c("Training", "Testing"), c("Testing", "2 transfers"), c("2 transfers", "10 transfers"), c("10 transfers", "20 transfers"), c("20 transfers", "Novel"), c("Novel", "Diff. species"), c("Diff. species", "Random"))
  data_subset <- subset(df_all, Unknown_token == "With unknown")
  p <- ggplot(data_subset, aes(x = Query, y=sum, colour = Query)) +
    geom_boxplot() +
    scale_colour_npg() +
    facet_grid(DatasetName ~ .) +
    theme_light() +
    stat_compare_means(
      comparisons = pairwise_comparisons,
      method = "wilcox.test",
      paired = FALSE,
      label = "p.signif",
      p.adjust.method = "bonferroni",
      size = 5.5,
      label.x = 1.0,
      symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
        symbols = c("****", "***", "**", "*", "ns")
      )
    ) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust= 1), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
    guides(colour=guide_legend(title="Dataset size")) +
    labs(x = "Query type", y = "Log Pseudolikelihood") +
    geom_hline(yintercept = 0, linetype="dashed")
  p
  ggsave(file=paste(outpref, "_pseudolikelihood_comparison_w_unknown_random.svg", sep = ""), plot=p, height = 10, width = 8)
  ggsave(file=paste(outpref, "_pseudolikelihood_comparison_w_unknown_random.png", sep = ""), plot=p, height = 10, width = 8)
  
  data_subset <- subset(df_all, Unknown_token != "With unknown")
  p <- ggplot(data_subset, aes(x = Query, y=sum, colour = Query)) +
    geom_boxplot() +
    scale_colour_npg() +
    facet_grid(Data ~ .) +
    theme_light() +
    stat_compare_means(paired = FALSE, size = 5.5, label.x = 1.0, comparisons = pairwise_comparisons, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
    guides(colour=guide_legend(title="Dataset size")) +
    labs(x = "Query type", y = "Log Pseudolikelihood") +
    geom_hline(yintercept = 0, linetype="dashed")
  p
  ggsave(file=paste(outpref, "_pseudolikelihood_comparison_wo_unknown.svg", sep = ""), plot=p, height = 10, width = 8)
  ggsave(file=paste(outpref, "_pseudolikelihood_comparison_wo_unknown.png", sep = ""), plot=p, height = 10, width = 8)
}

species <- opt$species
indir <- file.path(opt$indir, species)
outdir <- file.path(opt$outdir, species)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


{
  df_paths <- list.files(path = indir, pattern = "*_pseudolikelihood.txt", full.names = TRUE)
  outpref <- file.path(outdir, species)
  plot_results(df_paths, outpref)
}