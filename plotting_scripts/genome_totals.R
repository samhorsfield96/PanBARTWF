library(dplyr)
library(ggplot2)
library(optparse)

option_list <- list(
  make_option("--root", type = "character", help = "Root input directory containing species subdirectories."),
  make_option("--outdir", type = "character", help = "Directory for generated tables and plots."),
  make_option("--species", type = "character", default = "E_coli,S_pneumoniae",
              help = "Comma-separated species subdirectory names."),
  make_option("--species-names", type = "character", default = "E. coli,S. pneumoniae",
              help = "Comma-separated species display names.")
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
    Validation = if (grepl("_val_", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out = if (grepl("_held_out", filename)) {
      TRUE
    } else {
      FALSE
    },
    Held_out_removed = if (grepl("_held_out_removed", filename)) {
      TRUE
    } else {
      FALSE
    },
    Smitis = if (grepl("mitis", filename)) {
      TRUE
    } else {
      FALSE
    }
  )
}

plot_results <- function(df_paths, outdir, species)
{
  j <- 5
  for (j in 1:length(df_paths))
  {
    filename <- df_paths[j]
    df <- read.table(filename, sep = ",", comment.char = "", header = TRUE)
    
    if (nrow(df) == 0) {
      next
    }
    
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
  
  df_all["Species_model"] <- species
  df_all['Type'] <- NA
  df_all['Data'] <- NA
  df_all['Lineages'] <- NA
  
  df_all['Type'][df_all['Train'] == TRUE] <- "Training"
  df_all['Type'][df_all['Test'] == TRUE] <- "Testing"
  df_all['Type'][df_all['Validation'] == TRUE] <- "Validation"
  df_all['Type'][df_all['Held_out'] == TRUE & df_all['Held_out_removed'] == FALSE] <- "Novel"
  df_all['Type'][df_all['Smitis'] == TRUE] <- "S_mitis"
  
  df_all['Data'][df_all['Sampled'] == TRUE] <- "10% Data"
  df_all['Data'][df_all['Sampled'] == FALSE] <- "100% Data"
  
  df_all['Lineages'][df_all['Held_out'] == TRUE] <- "Subset"
  df_all['Lineages'][df_all['Held_out'] == FALSE] <- "All"
  
  df_all = subset(df_all, select = -c(Sampled,Train,Test,Validation,Held_out,Smitis,Held_out_removed))
  
  df_all <- df_all[, c(3, 1, 2, 4, 5, 6)]
  df_all
}

root <- opt$root
outdir <- opt$outdir
species_list <- strsplit(opt$species, ",", fixed = TRUE)[[1]]
species_name_list <- strsplit(opt$`species-names`, ",", fixed = TRUE)[[1]]
if (length(species_list) != length(species_name_list)) {
  stop("--species and --species-names must contain the same number of values.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species_results <- lapply(seq_along(species_list), function(i) {
  indir <- file.path(root, species_list[[i]])
  df_paths <- list.files(path = indir, pattern = "_clusters\\.txt$", full.names = TRUE, recursive = TRUE)
  plot_results(df_paths, outdir, species_name_list[[i]])
})
df_total <- do.call(rbind, species_results)

df_sum <- df_total %>%
  group_by(Species_model, Type, Data, Lineages) %>%
  summarise(Count = n())

write.csv(df_total, file = file.path(outdir, "sample_ids.csv"), row.names = FALSE)

write.csv(df_sum, file = file.path(outdir, "count_summary.csv"), row.names = FALSE)

# generate figures showing cluster size distribution
df_cluster_sum <- df_total %>%
  group_by(Species_model, Cluster, Type, Data, Lineages) %>%
  summarise(Count = n())

df_cluster_sum_subset <- subset(df_cluster_sum, Lineages == "All" & Type != "S_mitis")

df_cluster_sum_subset <- df_cluster_sum_subset %>%
  group_by(Species_model, Cluster, Data) %>%
  summarise(Count = sum(Count))

df_cluster_sum_subset_max <- df_cluster_sum_subset %>%
  group_by(Species_model, Data) %>%
  summarise(Max = max(Count),
            Min = min(Count),
            Median = median(Count))

p <- ggplot(data = df_cluster_sum_subset, aes(x = Count, after_stat(count))) +
  geom_density() +
  facet_wrap(Species_model ~ Data) +
  scale_x_log10() +
  theme_light() +
  labs(x = "Cluster size", y = "Count") 
p

ggsave(file=file.path(outdir, "cluster_size_density.png"), plot=p, height = 6, width = 8)
