library(dplyr)
# data = readRDS("~/Bureau/Work/PhD/conferences/HADACA3/starting_kit_phase1/data/mixes1_SDE5_pdac.rds")
reference = readRDS("~/Bureau/Work/PhD/conferences/HADACA3/starting_kit_phase1/data/reference_pdac.rds")

# rna_df = data1[["mix_rna"]] %>% as.matrix()
# met_df = data1[["mix_met"]] %>% as.matrix()
rna_ref = reference[["ref_bulkRNA"]] %>% as.matrix()
met_ref = reference[["ref_met"]] %>% as.matrix()

# Identify the fraction of cell types where each gene is expressed (positive expression)
unexpressed_fraction <- apply(rna_ref, 1, function(x) {
  sum(x == 0) / length(x)  # Count positive values and divide by total cell types
})

gene_variances <- apply(rna_ref %>% scale(), 1, var)

# Select genes only expressed in a low fraction of cell types (most discriminatory?)
filter_genes <- data.frame("gene" = names(unexpressed_fraction), 
                           "fraction" = unexpressed_fraction, 
                           "variance" = gene_variances,
                           "varrank" = rank(-gene_variances)) %>% 
  arrange(desc(fraction))
# High variance : top 1000 by variance
# Medium variance : top 3000 by variance
# Low variance : top 5000 by variance
# High fraction : expressed in only one cell type
# Medium fraction : expressed in up to 2 cell types
# Low fraction : expressed in up to 4 cell types
genelist_highfraction_highvariance = (filter_genes %>% filter((fraction > 0.6 & fraction < 1) | varrank < 1000) %>% arrange(varrank))$gene

genelist_mediumfraction_highvariance = (filter_genes %>% filter((fraction > 0.2 & fraction < 1) | varrank < 1000) %>% arrange(varrank))$gene

genelist_lowfraction_highvariance = (filter_genes %>% filter((fraction > 0 & fraction < 1) | varrank < 1000) %>% arrange(varrank))$gene

genelist_highfraction_mediumvariance = (filter_genes %>% filter((fraction > 0.6 & fraction < 1) | varrank < 3000) %>% arrange(varrank))$gene

genelist_mediumfraction_mediumvariance = (filter_genes %>% filter((fraction > 0.2 & fraction < 1) | varrank < 3000) %>% arrange(varrank))$gene

genelist_lowfraction_mediumvariance = (filter_genes %>% filter((fraction > 0 & fraction < 1) | varrank < 3000) %>% arrange(varrank))$gene

genelist_highfraction_lowvariance = (filter_genes %>% filter((fraction > 0.6 & fraction < 1) | varrank < 5000) %>% arrange(varrank))$gene

genelist_mediumfraction_lowvariance = (filter_genes %>% filter((fraction > 0.2 & fraction < 1) | varrank < 5000) %>% arrange(varrank))$gene

genelist_lowfraction_lowvariance = (filter_genes %>% filter((fraction > 0 & fraction < 1) | varrank < 5000) %>% arrange(varrank))$gene

met_variances <- apply(met_ref %>% scale(), 1, var)

# Filter CpG sites by variance only
# High variance : top 1000 by variance
# Medium variance : top 3000 by variance
# Low variance : top 5000 by variance
filter_cpg <- data.frame("site" = rownames(met_ref), 
                           "variance" = met_variances,
                           "varrank" = rank(-met_variances)) %>% 
  arrange(varrank)
metlist_highvariance = (filter_cpg %>% filter(varrank < 1000))$site

metlist_mediumvariance = (filter_cpg %>% filter(varrank < 3000))$site

metlist_lowvariance = (filter_cpg %>% filter(varrank < 5000))$site
