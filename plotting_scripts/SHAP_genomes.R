library(optparse)

option_list <- list(
  make_option("--infile-genome", type = "character", help = "Path to SHAP genomes file."),
  make_option("--infile-clusters", type = "character", help = "Path to cluster assignment CSV file."),
  make_option("--outfile", type = "character", default = "SHAP_genomes.csv",
              help = "Output CSV filename. Default = 'SHAP_genomes.csv'.")
)
opt <- parse_args(OptionParser(option_list = option_list))

genome.df <- read.csv(opt$`infile-genome`, header = FALSE)
genome.df <- data.frame(Taxon = unique(genome.df$V1))
cluster.assignment <- read.csv(opt$`infile-clusters`)

merged.df <- merge(genome.df, cluster.assignment, by = "Taxon")
merged.df <- subset(merged.df, Lineages == "All")

write.csv(merged.df, file = opt$outfile, row.names = FALSE)
