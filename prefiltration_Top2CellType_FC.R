library(dplyr)
data = readRDS("~/Bureau/Work/PhD/conferences/HADACA3/starting_kit_phase1/data/mixes1_SDE5_pdac.rds")
reference = readRDS("~/Bureau/Work/PhD/conferences/HADACA3/starting_kit_phase1/data/reference_pdac.rds")

mix_rna = data[["mix_rna"]] %>% as.matrix()
mix_met = data[["mix_met"]] %>% as.matrix()
ref_bulkRNA = reference[["ref_bulkRNA"]] %>% as.matrix()
ref_met = reference[["ref_met"]] %>% as.matrix()

# Identify genes with highest FC between first and second most expressed cell types. 
# Calculate the fold change between the top two cell types for each gene
calculate_fold_change <- function(x) {
  sorted_values <- sort(x, decreasing = TRUE) # Sort expression values for the gene
  if (length(sorted_values) > 1) {
    return(sorted_values[1] / sorted_values[2]) # Fold change: top value / second top value
  } else {
    return(NA) # Handle cases where there are less than 2 cell types
  }
}

# Apply the function to each row of the matrix
fold_changes <- apply(ref_bulkRNA, 1, calculate_fold_change)

# Create a dataframe of genes and their fold changes
fold_change_df <- data.frame(
  Gene = rownames(ref_bulkRNA),
  FoldChange = fold_changes
)

# Filter for non-NA fold changes and sort by fold change in descending order
top_genes <- fold_change_df[!is.na(fold_change_df$FoldChange), ]
top_genes <- top_genes[order(-top_genes$FoldChange), ]

# Select the top 500 genes
top_genes_FC <- head(top_genes, 500)

genelist_top500FC = top_genes_FC$Gene

# Identify CpG Sites with highest variation

filter_cpg <- data.frame("site" = rownames(ref_met), 
                         "variance" = met_variances,
                         "varrank" = rank(-met_variances)) %>% 
  arrange(varrank)
metlist_highvariance = (filter_cpg %>% filter(varrank < 100))$site

metlist_mediumvariance = (filter_cpg %>% filter(varrank < 300))$site

metlist_lowvariance = (filter_cpg %>% filter(varrank < 500))$site