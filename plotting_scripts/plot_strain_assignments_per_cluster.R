library(ggplot2)
library(ggsci)
library(ggpubr)
library(rlang)
library(tools)
library(mgcv)
library(optparse)

option_list <- list(
  make_option("--indir", type = "character", help = "Root directory containing species assignment subdirectories."),
  make_option("--species", type = "character", default = "E_coli", help = "Species subdirectory name."),
  make_option("--outpref", type = "character", default = "output", help = "Output prefix for saved plots.")
)
opt <- parse_args(OptionParser(option_list = option_list))

search_string <- function(filename, string){
  search_rg <- paste0(".*", string, "_(-?[0-9.]+(?:[eE][-+]?[0-9]+)?)(?:_|\\.|$).*")
  label = if (grepl(string, filename)) {
    as.numeric(sub(search_rg, "\\1", filename))
  } else {
    NA
  }
  return(label)
}

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
    Query = if (grepl("_query_", filename)) {
      TRUE
    } else {
      FALSE
    },
    S_mitis = if (grepl("_S_mitis_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Sketchlib = if (grepl("sketchlib", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

parse_df <- function(df_paths) {
  
  j <- 27
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = "\t", comment.char = "", header = TRUE)
    
    colnames(df)[which(names(df) == "Label")] <- "Cluster"
    df <- df[!df$Cluster == "overall", ]
    df <- df[!df$Cluster == "overall_balanced", ]
    
    df$filename <- parsed_filename <- file_path_sans_ext(basename(filename))
    
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
  
  df_all['Data'] <- NA
  df_all['Lineages'] <- NA
  df_all['Unknown_token'] <- NA
  df_all['Classifier'] <- NA
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10%"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100%"
  
  df_all['Lineages'][df_all['Held_out'] == TRUE] <- "Subset lineages"
  df_all['Lineages'][df_all['Held_out'] == FALSE] <- "All lineages"
  df_all['Lineages'][df_all['Held_out_only'] == TRUE] <- "Held-out lineages"
  df_all['Lineages'][df_all['Held_out'] == FALSE & df_all['Query'] == TRUE ] <- "All lineages (Query)"
  df_all['Lineages'][df_all['Held_out'] == FALSE & df_all['S_mitis'] == TRUE ] <- "All lineages (S. mitis)"
  
  df_all['Unknown_token'][df_all['With_unk'] == TRUE] <- "With unknown"
  df_all['Unknown_token'][df_all['With_unk'] == FALSE] <- "Without unknown"
  
  df_all['Classifier'][df_all['Sketchlib'] == TRUE] <- "Sketchlib"
  df_all['Classifier'][df_all['Sketchlib'] == FALSE] <- "PanBART"
  
  
  df_all$Data <- factor(df_all$Data, levels = c("10%", "100%"))
  df_all$Lineages <- factor(df_all$Lineages, levels = c("All lineages", "Subset lineages", "Held-out lineages", "All lineages (Query)", "All lineages (S. mitis)"))
  df_all$Classifier <- factor(df_all$Classifier, levels = c("PanBART", "Sketchlib"))
  
  df_all
}

generate_plots <- function(data_subset, outpref)
{
  gam_results <- data_subset %>%
    filter(
      is.finite(Test_count),
      Test_count > 0,
      is.finite(Accuracy)
    ) %>%
    group_by(K, Data, Classifier) %>%
    group_modify(~ {
      
      # No meaningful smooth-term test if Accuracy is constant
      if (n_distinct(.x$Accuracy) < 2) {
        return(tibble(
          R2 = NA_real_,
          p_value = NA_real_
        ))
      }
      
      # Fit GAM and safely handle errors
      tryCatch({
        model <- gam(
          Accuracy ~ s(log10(Test_count), k = 4),
          data = .x,
          method = "REML"
        )
        
        sm <- summary(model)
        
        tibble(
          R2 = sm$r.sq,
          p_value = sm$s.table[1, "p-value"]
        )
        
      }, error = function(e) {
        warning(
          "GAM failed for K = ", .y$K,
          ", Data = ", .y$Data,
          ", Classifier = ", .y$Classifier,
          ": ", conditionMessage(e)
        )
        
        tibble(
          R2 = NA_real_,
          p_value = NA_real_
        )
      })
    }) %>%
    ungroup()
  
  # Create annotation labels
  gam_labels <- gam_results %>%
    filter(is.finite(R2), is.finite(p_value)) %>%
    group_by(K, Data) %>%
    arrange(Classifier, .by_group = TRUE) %>%
    mutate(
      label = sprintf("R² = %.3f, p = %.3g", R2, p_value),
      y_position = 0.02 + (row_number() - 1) * 0.5
    ) %>%
    ungroup()
  
  p <- ggplot(data_subset, aes(x = Test_count, y=Accuracy, colour = Classifier)) +
    geom_point(alpha=0.1) +
    geom_smooth(method = "gam")+
    scale_colour_npg() +
    theme_light() +
    facet_grid(K ~ Data) +
    geom_text(
      data = gam_labels,
      aes(
        x = Inf,
        y = y_position,
        label = label,
        colour = Classifier
      ),
      inherit.aes = FALSE,
      hjust = 1.1,
      vjust = 0,
      size = 3
    ) +
    scale_x_log10() +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=16,face="bold"), legend.title=element_text(size=16,face="bold"), legend.text=element_text(size=14), strip.text = element_text(size=12)) +
    labs(x = "Cluster size", y = "Accuracy") +
    geom_hline(yintercept = 1, linetype="dashed")
  p
  ggsave(file=paste(outpref, "_frequency.png", sep = ""), plot=p, height = 10, width = 8)
  
  p_bal_accuracy <- ggplot(data_subset, aes(x = K, y=Accuracy, colour = Classifier)) +
    geom_smooth(method = "loess") +
    scale_colour_npg() +
    theme_light() +
    facet_wrap(Data ~ ., ncol=4) +
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
    geom_hline(yintercept = 1, linetype="dashed")
  p_bal_accuracy
  ggsave(file=paste(outpref, "_LOESS_accuracy.png", sep = ""), plot=p_bal_accuracy, height = 3, width = 6)
}

plot_results <- function(df_paths_panbart, df_paths_sketchlib, outpref) {
  df_panbart <- parse_df(df_paths_panbart)
  df_sketchlib <- parse_df(df_paths_sketchlib)
  
  df_merged <- rbind(df_panbart, df_sketchlib)
  
  df_test <- subset(df_merged, Classifier == "Sketchlib")
  
  # all data, without unknown
  
  # all data, with unknown
  data_subset = subset(df_merged, Lineages == "All lineages" & ((Classifier == "PanBART" & Unknown_token != "Without unknown") | (Classifier == "Sketchlib" & Unknown_token == "Without unknown")))
  outpref_plot <- paste0(outpref, "_all_lineages_w_unknown")
  generate_plots(data_subset, outpref_plot)
  
  # Subset lineages without unknown
  
  # Subset lineages with unknown
  data_subset = subset(df_merged, Lineages == "Subset lineages" & ((Classifier == "PanBART" & Unknown_token != "Without unknown") | (Classifier == "Sketchlib" & Unknown_token == "Without unknown")))
  outpref_plot <- paste0(outpref, "_held_out_removed_w_unknown")
  generate_plots(data_subset, outpref_plot)
  
  # Held-out lineages without unknown
  
  # Held-out lineages with unknown
}

species <- opt$species
indir <- file.path(opt$indir, species)
outpref <- opt$outpref
dir.create(dirname(outpref), recursive = TRUE, showWarnings = FALSE)

{
  df_paths_panbart <- list.files(path = file.path(indir, "panbart"), pattern = "per_label_accuracy\\.tsv$", full.names = TRUE)
  df_paths_sketchlib <- list.files(path = file.path(indir, "sketchlib"), pattern = "per_label_accuracy\\.tsv$", full.names = TRUE)
  df_paths_bacformer <- list.files(path = file.path(indir, "bacformer"), pattern = "per_label_accuracy\\.tsv$", full.names = TRUE)
  df_paths_bacformer_complete <- list.files(path = file.path(indir, "bacformer_complete"), pattern = "per_label_accuracy\\.tsv$", full.names = TRUE)
  plot_results(df_paths_panbart, df_paths_sketchlib, outpref)
}
