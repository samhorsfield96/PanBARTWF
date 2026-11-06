library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(optparse)

option_list <- list(
  make_option("--infile", type = "character", help = "Path to attention summary CSV file."),
  make_option("--outdir", type = "character", help = "Directory for generated plots.")
)
opt <- parse_args(OptionParser(option_list = option_list))

top_n_genes <- 20
# how to aggregate attention_from_locus per gene: sum or mean
agg_fun <- sum

df <- read_csv(opt$infile, show_col_types = FALSE,
               col_types = cols(.default = col_guess(), annot_Gene = col_character()))

# keep only entries with a gene annotation
df <- df %>% filter(!is.na(annot_Gene), annot_Gene != "")

# overall most frequent genes, used as a common ordered x axis for all panels
gene_order <- df %>%
  group_by(annot_Gene) %>%
  summarise(n = agg_fun(attention_from_locus, na.rm = TRUE)) %>%
  arrange(desc(n)) %>%
  slice_head(n = top_n_genes) %>%
  pull(annot_Gene)

counts <- df %>%
  filter(annot_Gene %in% gene_order) %>%
  group_by(attention_type, layer, attention_head, annot_Gene) %>%
  summarise(n = agg_fun(attention_from_locus, na.rm = TRUE), .groups = "drop") %>%
  mutate(annot_Gene = factor(annot_Gene, levels = gene_order),
         layer = factor(layer),
         attention_head = factor(attention_head))

plot_type <- function(type) {
  ggplot(filter(counts, attention_type == type),
         aes(x = annot_Gene, y = n)) +
    geom_col() +
    facet_grid(layer ~ attention_head, labeller = label_both) +
    labs(x = "Gene", y = "Attention from locus") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

types <- c("cross", "encoder", "decoder")
plots <- lapply(types, plot_type)
names(plots) <- types

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
for (type in types) {
  ggsave(file.path(opt$outdir, paste0("attention_summary_", type, ".png")), plots[[type]],
         width = 18, height = 10, dpi = 300)
}

combined <- ggarrange(plotlist = plots, ncol = 1, nrow = 3)
ggsave(file.path(opt$outdir, "attention_summary_all.png"), combined,
       width = 14, height = 30, dpi = 300, limitsize = FALSE)
