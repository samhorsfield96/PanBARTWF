library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(ggpubr)
library(reporter)
library(tools)
library(dplyr)
library(forcats)
library(reshape2)
library(tidyr)
library(dplyr)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Input directory containing species subdirectories."),
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
    Held_out = if (grepl("_held_out_removed", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out_only = if (grepl("_held_out_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Novel = if (grepl("_query_", filename)) {
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
    Smitis = if (grepl("_mitis_", filename)) {
      TRUE
    } else {
      FALSE
    },
    With_unk = if (!grepl("_without_unk_", filename)) {
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

plot_individual <- function(df_paths, df_paths_silhouette, outpref, method)
{
  j <- 22
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header = TRUE)
    # hold original cluster
    df$Cluster_ori <- df$Cluster
    
    # get silhouette values
    filename_silhouette <- df_paths_silhouette[j]
    df_silhouette <- read.table(filename_silhouette, sep = "\t", comment.char = "", header = TRUE)
    # remove rows with count 1
    df_silhouette <- df_silhouette[!df_silhouette$Count <= 1, ]
    df_silhouette <- df_silhouette[!df_silhouette$Label == "overall", ]
    median_silhouette <- median(df_silhouette$Silhouette_width)
    IQR_silhouette <- IQR(df_silhouette$Silhouette_width)
    mean_silhouette <- median(df_silhouette$Silhouette_width)
    stddev_silhouette <- sd(df_silhouette$Silhouette_width)
    
    df$Silhouette_median <- median_silhouette
    df$Silhouette_IQR <- IQR_silhouette
    df$Silhouette_mean <- mean_silhouette
    df$Silhouette_stddev <- stddev_silhouette
    
    parsed_filename <- file_path_sans_ext(basename(filename))
    df$parsed_filename <- parsed_filename
    
    # top 10 cluster sizes per dataset  
    n_top <- 10
    
    # sort columns
    top_classes <- df %>%
      count(Cluster, sort = TRUE) %>%
      top_n(n_top, n) %>%
      pull(Cluster)
    
    df <- df %>%
      mutate(Cluster = ifelse(Cluster %in% top_classes, Cluster, "0"),
             Cluster = factor(Cluster))
    
    df_silhouette <- df_silhouette %>%
      mutate(Label = ifelse(Label %in% top_classes, Label, "0"),
             Label = factor(Label))
    # Put "0" last instead of first
    df$Cluster <- fct_relevel(df$Cluster, "0", after = Inf)
    
    #Sort rows so "0" comes last (plotted on top)
    df <- df[order(df$Cluster != "0"), ]
    
    # get stats for top clusters only
    top_df_silhouette <- subset(df_silhouette, Label != "0")
    median_silhouette <- median(top_df_silhouette$Silhouette_width)
    IQR_silhouette <- IQR(top_df_silhouette$Silhouette_width)
    mean_silhouette <- median(top_df_silhouette$Silhouette_width)
    stddev_silhouette <- sd(top_df_silhouette$Silhouette_width)
    
    df$Silhouette_median_top <- median_silhouette
    df$Silhouette_IQR_top <- IQR_silhouette
    df$Silhouette_mean_top <- mean_silhouette
    df$Silhouette_stddev_top <- stddev_silhouette
    
    # Build a named colour vector:
    n_top <- nlevels(df$Cluster) - 1
    # brewer colours
    my_cols <- c(
      setNames(brewer.pal(min(12, n_top), "Set3")[1:n_top], levels(df$Cluster)[-nlevels(df$Cluster)]),
      "0" = "gray30"
    )
    
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) + geom_point(size=0.01, alpha =0.5) +
      theme_light() +
      scale_color_manual(values=my_cols) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
      annotate(
        "label",
        x = Inf, y = Inf,   # top right corner
        label = paste0("Median silhouette = ", round(median_silhouette, 2), "\n",
                       "IQR silhouette = ", round(IQR_silhouette, 2)),
        hjust = 1.1, vjust = 1.1,  # nudges inward from the edge
        label.size = 0.5,          # border thickness
        label.r = unit(0.15, "lines"), # corner rounding
        fill = "white",            # background color
        color = "black",            # text + border color
        size = 3
      ) +
      theme(
        panel.background = element_rect(fill = "black"), # plot background
        panel.grid.major = element_line(color = "gray40"), # major grid
        panel.grid.minor = element_line(color = "gray20")  # minor grid
      )
    p
    
    # add metadata
    parsed_values <- parse_filename(filename)
    for (name in names(parsed_values)) {
      df[[name]] <- parsed_values[[name]]
    }
    
    species <- df$Species[1]
    
    dir.create(file.path(outpref, species), recursive = TRUE, showWarnings = FALSE)
    ggsave(file=paste0(outpref, "/", species, "/", species, "_", method, "_", parsed_filename, ".png"), plot=p, height = 6, width = 6)
    
    if (j == 1) {
      df_all <- df
    } else 
    {
      df_all <- rbind(df_all, df)
    }
  }
  
  df_all['Unknown_token'] <- NA
  df_all['Data'] <- NA
  df_all['Subset'] <- NA
  df_all['Query'] <- NA
  
  df_all['Subset'][df_all['Held_out'] != TRUE ] <- "All lineages"
  df_all['Subset'][df_all['Held_out'] == TRUE | df_all['Held_out_only'] == TRUE] <- "Subset lineages"
  
  df_all['Unknown_token'][df_all['With_unk'] == TRUE] <- "With unknown"
  df_all['Unknown_token'][df_all['With_unk'] == FALSE] <- "Without unknown"
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10% Data"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100% Data"
  
  df_all['Query'][df_all['Smitis'] == TRUE] <- "S. mitis"
  df_all['Query'][df_all['Train'] == TRUE] <- "Training"
  df_all['Query'][df_all['Test'] == TRUE] <- "Testing"
  df_all['Query'][df_all['Held_out'] != TRUE & df_all['Held_out_only'] == TRUE] <- "Novel"
  df_all['Query'][df_all['Novel'] == TRUE & df_all['Train'] == TRUE] <- "Novel vs. Training"
  df_all['Query'][df_all['Novel'] == TRUE & df_all['Train'] != TRUE] <- "Novel"
  df_all['Query'][df_all['Test'] == TRUE & df_all['Train'] == TRUE] <- "Testing vs. Training"
  df_all['Query'][df_all['Smitis'] == TRUE & df_all['Train'] == TRUE ] <- "S. mitis vs. Training"
  df_all['Query'][df_all['Smitis'] == TRUE & df_all['Train'] != TRUE] <- "S. mitis"
  
  # add final row which assigns database name
  df_all$DatasetName <- assign_database_name(df_all)
  
  # Get one row per unique 'group'
  test_df <- df_all %>%
    distinct(parsed_filename, .keep_all = TRUE)
  
  df_all
}

plot_all <- function(df, outpref, plot_height, plot_width, point_size, point_alpha, svg=FALSE) {
    # top 10 cluster sizes per dataset  
    n_top <- 10
  
    # sort columns
    top_classes <- df %>%
      group_by(DatasetName) %>%
      count(Cluster, sort = TRUE) %>%
      top_n(n_top, n) %>%
      pull(Cluster)
    
    df <- df %>%
      group_by(DatasetName) %>%
      mutate(Cluster = ifelse(Cluster_ori %in% top_classes, Cluster_ori, "0"))
    
    # Put "0" last instead of first
    df$Cluster <- fct_relevel(df$Cluster, "0", after = Inf)
    
    # brewer colours
    df_ranked <- df %>%
      filter(Cluster != "0") %>%   # <-- add this
      group_by(Data, Classifier, Subset, Species, Cluster) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(Data, Classifier, Subset, Species) %>%
      slice_max(count, n = n_top) %>%
      mutate(rank = row_number()) %>%
      ungroup() %>%
      select(Data, Classifier, Subset, Species, Cluster, rank)
    
    # Join rank back onto main df
    df <- df %>%
      left_join(df_ranked, by = c("Data", "Classifier", "Subset", "Species", "Cluster")) %>%
      mutate(rank = ifelse(is.na(rank), 0, rank))  # non-top clusters get rank 0
    
    df <- df %>% arrange(desc(rank == 0), rank)
    
    set3_cols <- brewer.pal(12, "Set3")
    
    rank_cols <- c(
      setNames(
        set3_cols[(seq_along(1:n_top) - 1) %% 12 + 1],  # cycles through 12 colours
        as.character(1:n_top)
      ),
      "0" = "gray30"
    )
    
    silhouette <- df %>% 
      group_by(DatasetName, Classifier) %>%
      summarise(Silhouette_median_top = if (is.na(unique(Silhouette_median_top))) unique(Silhouette_median) else unique(Silhouette_median_top),
                Silhouette_IQR_top =if (is.na(unique(Silhouette_IQR_top))) unique(Silhouette_IQR) else unique(Silhouette_IQR_top),
                Silhouette_median_all = unique(Silhouette_median), Silhouette_IQR_all = unique(Silhouette_IQR))
    
    
    p <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = as.character(rank))) + geom_point(size=point_size, alpha =point_alpha) +
      theme_light() +
      facet_grid(DatasetName ~ Classifier) +
      scale_color_manual(values=rank_cols) +
      theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12), legend.position = "none") +
      geom_label(
        data = silhouette,
        aes(
          x = Inf, y = Inf,
          label = paste0(          "Median silhouette = ", round(Silhouette_median_top, 2), "\n",
          "IQR silhouette = ", round(Silhouette_IQR_top, 2)),
        ),
        inherit.aes = FALSE,
        hjust = 1.1, vjust = 1.1,
        fill = "white",
        color = "black",
        size = 3
      ) +
      theme(
        panel.background = element_rect(fill = "black"), # plot background
        panel.grid.major = element_line(color = "gray40"), # major grid
        panel.grid.minor = element_line(color = "gray20")  # minor grid
      )
    p
    
    ggsave(file=paste0(outpref, ".png"), plot=p, height = plot_height, width = plot_width)
    if (svg)
    {
      ggsave(file=paste0(outpref, ".svg"), plot=p, height = plot_height, width = plot_width)
    }
}

plot_results <- function(df_paths_panbart, df_paths_sketchlib, df_paths_panbart_silhouette, df_paths_sketchlib_silhouette, 
                         df_paths_bacformer, df_paths_bacformer_silhouette, df_paths_bacformer_complete, df_paths_bacformer_complete_silhouette, outpref)
{
  df_all_panbart <- plot_individual(df_paths_panbart, df_paths_panbart_silhouette, outpref, "panbart")
  df_all_panbart$Classifier <- "PanBART"
  df_all_sketchlib <- plot_individual(df_paths_sketchlib, df_paths_sketchlib_silhouette, outpref, "sketchlib")
  df_all_sketchlib$Classifier <- "Sketchlib"
  df_all_bacformer <- plot_individual(df_paths_bacformer, df_paths_bacformer_silhouette, outpref, "bacformer")
  df_all_bacformer$Classifier <- "Bacformer-L (MAG)"
  df_all_bacformer_complete <- plot_individual(df_paths_bacformer_complete, df_paths_bacformer_complete_silhouette, outpref, "bacformer")
  df_all_bacformer_complete$Classifier <- "Bacformer-L (Complete)"
  
  
  merged_df <- rbind(df_all_sketchlib, df_all_panbart, df_all_bacformer, df_all_bacformer_complete)
  merged_df$Classifier <- factor(merged_df$Classifier, levels = c("PanBART", "Sketchlib", "Bacformer-L (MAG)", "Bacformer-L (Complete)"))
  
  # training only, subset - fine
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset == "Subset lineages" & Smitis == FALSE)
  outpref_plot <- paste0(outpref, "/training_subset_lineages")
  plot_all(subset_df, outpref_plot, 14, 10, 0.5, 0.2)
  
  # plot 3 training only, all, 100% data - fine
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset != "Subset lineages" & Data == "100% Data" & Smitis == FALSE)
  outpref_plot <- paste0(outpref, "/training_all_lineages_100_data")
  plot_all(subset_df, outpref_plot, 10, 18, 0.001, 0.2, svg = TRUE)
  
  # training only, all, 10% data - fine
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset != "Subset lineages" & Data == "10% Data" & Smitis == FALSE)
  outpref_plot <- paste0(outpref, "/training_all_lineages_10_data")
  plot_all(subset_df, outpref_plot, 10, 10, 0.5, 0.2)
  
  # training only, all - fine
  subset_df <- subset(merged_df, Query == "Training" & Unknown_token == "With unknown" & Subset != "Subset lineages" & Smitis == FALSE)
  outpref_plot <- paste0(outpref, "/training_all_lineages")
  plot_all(subset_df, outpref_plot, 14, 10, 0.001, 0.2)
  
  # testing only, subset - fine
  subset_df <- subset(merged_df, Query == "Testing" & Unknown_token == "With unknown" & Subset == "Subset lineages" & Smitis == FALSE)
  outpref_plot <- paste0(outpref, "/testing_subset_lineages")
  plot_all(subset_df, outpref_plot, 14, 10, 0.5, 0.2)
  
  # testing only, all - fine
  subset_df <- subset(merged_df, Query == "Testing" & Unknown_token == "With unknown" & Subset != "Subset lineages")
  outpref_plot <- paste0(outpref, "/testing_all_lineages")
  plot_all(subset_df, outpref_plot, 14, 10, 0.5, 0.2)
  
  # held-out only, with unknown - fine
  subset_df <- subset(merged_df, Query == "Novel" & Unknown_token == "With unknown" & Subset == "Subset lineages" & Smitis == FALSE)
  outpref_plot <- paste0(outpref, "/novel_subset_lineages_with_unknown")
  plot_all(subset_df, outpref_plot, 14, 10, 0.5, 0.5)
  
  # held only, without unknown - fine
  subset_df <- subset(merged_df, (Query == "Novel" & Smitis == FALSE & Unknown_token != "With unknown" & Subset == "Subset lineages") | (Query == "Novel" & Classifier == "Sketchlib" & Subset == "Subset lineages"))
  outpref_plot <- paste0(outpref, "/novel_subset_lineages_without_unknown")
  plot_all(subset_df, outpref_plot, 14, 10, 0.75, 0.5)
  
  # held out vs training - fine
  subset_df <- subset(merged_df, (Query == "Novel vs. Training" & Unknown_token != "With unknown" & Smitis == FALSE))
  outpref_plot <- paste0(outpref, "novel_v_training_without_unknown")
  plot_all(subset_df, outpref_plot, 14, 5, 0.75, 0.5)
  
  # held out vs training - fine
  subset_df <- subset(merged_df, (Query == "Novel vs. Training" & Unknown_token == "With unknown" & Smitis == FALSE))
  outpref_plot <- paste0(outpref, "novel_v_training_with_unknown")
  plot_all(subset_df, outpref_plot, 14, 5, 0.75, 0.5)
  
  # S mitis vs training
  subset_df <- subset(merged_df, (Query == "S. mitis vs. Training" & Unknown_token != "With unknown"))
  outpref_plot <- paste0(outpref, "/mitis_v_training_without_unknown")
  plot_all(subset_df, outpref_plot, 7, 5, 0.001, 0.2)
  
}

indir <- opt$indir
outdir <- opt$outdir
species_list <- strsplit(opt$species, ",", fixed = TRUE)[[1]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

df_paths_panbart <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "panbart"),
    pattern = "*_UMAP.csv",
    full.names = TRUE
  )
}))
df_paths_panbart_silhouette <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "panbart"),
    pattern = "*_silhouette.tsv",
    full.names = TRUE
  )
}))
df_paths_sketchlib <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "sketchlib"),
    pattern = "*_UMAP.csv",
    full.names = TRUE
  )
}))
df_paths_sketchlib_silhouette <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "sketchlib"),
    pattern = "*_silhouette.tsv",
    full.names = TRUE
  )
}))
df_paths_bacformer <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "bacformer"),
    pattern = "*_UMAP.csv",
    full.names = TRUE
  )
}))
df_paths_bacformer_silhouette <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "bacformer"),
    pattern = "*_silhouette.tsv",
    full.names = TRUE
  )
}))
df_paths_bacformer_complete <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "bacformer_complete"),
    pattern = "*_UMAP.csv",
    full.names = TRUE
  )
}))
df_paths_bacformer_complete_silhouette <- unlist(lapply(species_list, function(sp) {
  list.files(
    path = file.path(indir, sp, "bacformer_complete"),
    pattern = "*_silhouette.tsv",
    full.names = TRUE
  )
}))

outpref <- outdir

plot_results(df_paths_panbart, df_paths_sketchlib, df_paths_panbart_silhouette, df_paths_sketchlib_silhouette, 
             df_paths_bacformer, df_paths_bacformer_silhouette, df_paths_bacformer_complete, df_paths_bacformer_complete_silhouette, outpref)
