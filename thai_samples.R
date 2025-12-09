#clear working directory:
rm(list=ls())

#################Load packages+ set directory###################################
library("tidyverse")
library("ggplot2")
library("BiocManager")
library("phyloseq")
library("dplyr")
library("reshape2")
library("ggpubr")
library("microbiome")
library("rstatix")
library("vegan")
library("pals")
library("microViz")#heatmaps,PCA
library("Polychrome")
library("ape")
library("corrplot")
library("microbial")#contains a LEFSE process among other things
library("ggrepel") # this is to make labels for ggplot that stick out
library("devtools")
library("ggforce")
library("gridExtra")
library("tibble")
library("DESeq2")
library("patchwork")
library("DECIPHER")
library("phangorn")
library("ComplexHeatmap")
library("ANCOMBC")
library("pheatmap")
library("RColorBrewer")
library("cowplot")
library("pairwiseAdonis")
library("microbiomeMarker")



setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/thai_samples")


########################IMPORT LIST##############################################
# Read the CSV files- create different object for each project code to be displayed within the plot:
otu_data_rosenbergii <- read.csv("OTU_rosenbergii.csv", row.names = 1, check.names = FALSE)
tax_data_rosenbergii <- read.csv("hiera_BLAST_rosenbergii.csv", row.names = 1, check.names = FALSE)
meta_data_rosenbergii <- read.csv("SraRunTable_rosenbergii.csv", row.names = 1, check.names = FALSE)

otu_data_vannamei <- read.csv("OTU_vannamei.csv", row.names = 1, check.names = FALSE)
tax_data_vannamei <- read.csv("hiera_BLAST_vannamei.csv", row.names = 1, check.names = FALSE)
meta_data_vannamei <- read.csv("SraRunTable_vannamei.csv", row.names = 1, check.names = FALSE)

otu_data_monodon <- read.csv("OTU_monodon.csv", row.names = 1, check.names = FALSE)
tax_data_monodon <- read.csv("hiera_BLAST_monodon.csv", row.names = 1, check.names = FALSE)
meta_data_monodon <- read.csv("SraRunTable_monodon.csv", row.names = 1, check.names = FALSE)



#Build the comp plot:

keep_bacteria <- tax_data_rosenbergii$Domain == "Bacteria"
otu_data_rosenbergii <- otu_data_rosenbergii[keep_bacteria, ]
tax_data_rosenbergii <- tax_data_rosenbergii[keep_bacteria, ]


# Convert OTU table - transform matrix 
otu_mat_rosenbergii <- as.matrix(otu_data_rosenbergii)
otu_table_rosenbergii <- otu_table(otu_mat_rosenbergii, taxa_are_rows = TRUE)

# Build sample data
sample_data_rosenbergii <- sample_data(meta_data_rosenbergii)

#load taxonomy data:
tax_table_rosenbergii <- tax_table(as.matrix(tax_data_rosenbergii))

colnames(otu_table_rosenbergii) <- rownames(sample_data_rosenbergii)

#create phyloseq object:
physeq_rosenbergii <- phyloseq(otu_table_rosenbergii, tax_table_rosenbergii, sample_data_rosenbergii)
physeq_rosenbergii # shows data 


nowater_phyr <- prune_samples(sample_data(physeq_rosenbergii)$Blastocystis != "Water",
                             physeq_rosenbergii)
print(sample_data(nowater_phyr))

#rarefactaion curve
curvedata<-as(t(otu_table(physeq_rosenbergii)), "matrix")
class(curvedata)
curve<-rarecurve(curvedata, step=100, col="hotpink")
print(curve)

#rarefy
physeq1.rosenbergii <- rarefy_even_depth(nowater_phyr, sample.size = 200000, rngseed = 123)

phyloseq_normalised_rosenbergii <- transform_sample_counts(physeq1.rosenbergii, function(otu) otu/sum(otu))

ntaxa(physeq1.rosenbergii)



#Build the comp plot:

keep_bacteria <- tax_data_vannamei$Domain == "Bacteria"
otu_data_vannamei <- otu_data_vannamei[keep_bacteria, ]
tax_data_vannamei <- tax_data_vannamei[keep_bacteria, ]


# Convert OTU table - transform matrix 
otu_mat_vannamei <- as.matrix(otu_data_vannamei)
otu_table_vannamei <- otu_table(otu_mat_vannamei, taxa_are_rows = TRUE)

# Build sample data
sample_data_vannamei <- sample_data(meta_data_vannamei)

#load taxonomy data:
tax_table_vannamei <- tax_table(as.matrix(tax_data_vannamei))

colnames(otu_table_vannamei) <- rownames(sample_data_vannamei)

#create phyloseq object:
physeq_vannamei <- phyloseq(otu_table_vannamei, tax_table_vannamei, sample_data_vannamei)
physeq_vannamei # shows data 


nowater_phyv <- prune_samples(sample_data(physeq_vannamei)$Blastocystis != "Water",
                              physeq_vannamei)
print(sample_data(nowater_phyv))

#rarefactaion curve
curvedata<-as(t(otu_table(physeq_vannamei)), "matrix")
class(curvedata)
curve<-rarecurve(curvedata, step=100, col="purple")
print(curve)

#rarefy
physeq1.vannamei <- rarefy_even_depth(nowater_phyv, sample.size = 200000, rngseed = 123)


phyloseq_normalised_vannamei <- transform_sample_counts(physeq1.vannamei, function(otu) otu/sum(otu))

ntaxa(physeq1.vannamei)


#Build the comp plot:

keep_bacteria <- tax_data_monodon$Domain == "Bacteria"
otu_data_monodon <- otu_data_monodon[keep_bacteria, ]
tax_data_monodon <- tax_data_monodon[keep_bacteria, ]


# Convert OTU table - transform matrix 
otu_mat_monodon <- as.matrix(otu_data_monodon)
otu_table_monodon <- otu_table(otu_mat_monodon, taxa_are_rows = TRUE)

# Build sample data
sample_data_monodon <- sample_data(meta_data_monodon)

#load taxonomy data:
tax_table_monodon <- tax_table(as.matrix(tax_data_monodon))

colnames(otu_table_monodon) <- rownames(sample_data_monodon)

#create phyloseq object:
physeq_monodon <- phyloseq(otu_table_monodon, tax_table_monodon, sample_data_monodon)
physeq_monodon # shows data 


nowater_phym <- prune_samples(sample_data(physeq_monodon)$Blastocystis != "Water",
                              physeq_monodon)
print(sample_data(nowater_phym))

#rarefactaion curve
curvedata<-as(t(otu_table(physeq_monodon)), "matrix")
class(curvedata)
curve<-rarecurve(curvedata, step=100, col="cyan")
print(curve)

#rarefy
physeq1.monodon <- rarefy_even_depth(nowater_phym, sample.size = 200000, rngseed = 123)

phyloseq_normalised_monodon <- transform_sample_counts(physeq1.monodon, function(otu) otu/sum(otu))

ntaxa(physeq1.monodon)










# helper to process a phyloseq to Phylum-level compositional phyloseq 
process_to_phylum <- function(physeq_obj, mean_threshold = 0.005) {
  relab <- phyloseq::transform_sample_counts(physeq_obj, function(x) x / sum(x))
  subs <- phyloseq::filter_taxa(relab, function(x) mean(x) > mean_threshold, prune = FALSE)
  phy_com <- phyloseq::prune_taxa(taxa = subs, x = physeq_obj)
  phy_com_p <- phy_com %>%
    microbiome::aggregate_taxa(level = "Phylum") %>%
    microbiome::transform(transform = "compositional")
  return(phy_com_p)
}

# process each organism
phy_v <- process_to_phylum(physeq1.vannamei)
phy_r <- process_to_phylum(physeq1.rosenbergii)
phy_m <- process_to_phylum(physeq1.monodon)

# 1) get the union of Phylum taxa names across all three processed objects
get_phyla_names <- function(phy) {
  # aggregated phyloseq should have tax_table with Phylum in rank name "Phylum"
  tt <- tax_table(phy)
  # Try to find the correct column named "Phylum" (or first column if different)
  if ("Phylum" %in% colnames(tt)) {
    as.character(tt[, "Phylum"])
  } else {
    # fallback: use taxonomic ranks column 1
    as.character(tt[, 1])
  }
}
phyla_all <- unique(c(get_phyla_names(phy_v), get_phyla_names(phy_r), get_phyla_names(phy_m)))
phyla_all <- phyla_all[!is.na(phyla_all)]
# optional: sort for nicer legend ordering
phyla_all <- sort(phyla_all)

# 2) build a named color vector for all phyla
n_phyla <- length(phyla_all)
# choose a palette - use RColorBrewer palettes and expand if needed
base_cols <- brewer.pal(min(max(3, n_phyla), 12), "Pastel1")  
if (n_phyla > length(base_cols)) {
  my_cols <- colorRampPalette(brewer.pal(12, "Pastel1"))(n_phyla)
} else {
  my_cols <- base_cols[seq_len(n_phyla)]
}

names(my_cols) <- phyla_all

phylum_plot <- function(physeq_obj, title_text, palette_named) {
  p <- plot_composition(
    physeq_obj,
    otu.sort = "abundance",
    average_by = "Organism"
  ) +
    ggtitle(title_text) +
    ylab("% Relative Abundance") +
    xlab("Samples") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 0),
    ) +
    # apply the common fill scale - this is a ggplot layer, so ok to add here
    scale_fill_manual(values = palette_named, na.value = "grey50", drop = FALSE) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
  
  return(p)
}


# 1) Make a legend-only plot that contains all phyla as factor levels
legend_df <- data.frame(
  Phylum = factor(names(my_cols), levels = names(my_cols)),
  x = 1, y = 1
)


legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    name   = "Phylum",
    values = my_cols,
    breaks = names(my_cols),
    drop   = FALSE
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 8)
  )


# extract the legend grob (this legend contains ALL phyla)
legend_grob <- cowplot::get_legend(legend_plot)

# 2) Make your three plots WITHOUT legends (as you already have)
p_v_noleg <- phylum_plot(phy_v, "vannamei",    my_cols) + theme(legend.position = "none")
p_r_noleg <- phylum_plot(phy_r, "rosenbergii", my_cols) + theme(legend.position = "none")
p_m_noleg <- phylum_plot(phy_m, "monodon",     my_cols) + theme(legend.position = "none")

# 3) Combine the three plots and attach the single legend on the right
left_row <- cowplot::plot_grid(p_v_noleg, p_r_noleg, p_m_noleg, ncol = 3, align = "v", axis = "lr")
final_plot <- cowplot::plot_grid(left_row, legend_grob, rel_widths = c(3, 0.7))

print(final_plot)






# --------------------------
# Generalized processing helper
# --------------------------
process_to_level <- function(physeq_obj, tax_level = "Phylum", mean_threshold = 0.005) {
  relab <- phyloseq::transform_sample_counts(physeq_obj, function(x) x / sum(x))
  subs <- phyloseq::filter_taxa(relab, function(x) mean(x) > mean_threshold, prune = FALSE)
  phy_com <- phyloseq::prune_taxa(taxa = subs, x = physeq_obj)
  phy_com_l <- phy_com %>%
    microbiome::aggregate_taxa(level = tax_level) %>%
    microbiome::transform(transform = "compositional")
  return(phy_com_l)
}

process_to_level <- function(physeq_obj, tax_level = "Phylum", mean_threshold = 0.005, verbose = TRUE) {
  # Step 1: convert to relative abundance
  relab <- phyloseq::transform_sample_counts(physeq_obj, function(x) x / sum(x))
  
  # Step 2: filter by mean relative abundance
  taxa_means <- taxa_sums(relab) / nsamples(relab)
  keep_taxa <- taxa_names(relab)[taxa_means > mean_threshold]
  
  # Compute filtering stats
  n_before <- ntaxa(relab)
  n_after  <- length(keep_taxa)
  perc_kept <- (n_after / n_before) * 100
  
  # Average abundance retained (total % of reads kept)
  total_abund_retained <- sum(taxa_means[keep_taxa]) / sum(taxa_means) * 100
  
  # Step 3: prune and aggregate to chosen taxonomic level
  phy_filt <- prune_taxa(keep_taxa, physeq_obj)
  phy_agg <- phy_filt %>%
    microbiome::aggregate_taxa(level = tax_level) %>%
    microbiome::transform(transform = "compositional")
  
  # Step 4: Optionally report summary to console
  if (verbose) {
    message("Processing to ", tax_level, " level")
    message(" - Initial taxa: ", n_before)
    message(" - Retained taxa: ", n_after, " (", round(perc_kept, 1), "%)")
    message(" - Mean relative abundance threshold: ", mean_threshold, " (", mean_threshold*100, "%)")
    message(" - Total mean abundance retained: ", round(total_abund_retained, 1), "% of total abundance")
    message(" - Aggregated taxa (", tax_level, "): ", ntaxa(phy_agg))
    message("-------------------------------------------------------------")
  }
  
  # Step 5: Attach summary as attribute for easy retrieval
  attr(phy_agg, "filter_summary") <- list(
    tax_level = tax_level,
    mean_threshold = mean_threshold,
    n_before = n_before,
    n_after = n_after,
    perc_kept = perc_kept,
    abundance_retained = total_abund_retained,
    n_after_aggregated = ntaxa(phy_agg)
  )
  
  return(phy_agg)
}






# --------------------------
# Generate processed phyloseq objects
# --------------------------

phy_v_family <- process_to_level(physeq1.vannamei, "Family")
phy_r_family <- process_to_level(physeq1.rosenbergii, "Family")
phy_m_family <- process_to_level(physeq1.monodon, "Family")

phy_v_genus   <- process_to_level(physeq1.vannamei,   "Genus")
phy_r_genus   <- process_to_level(physeq1.rosenbergii, "Genus")
phy_m_genus   <- process_to_level(physeq1.monodon,     "Genus")

phy_v_species <- process_to_level(physeq1.vannamei,   "Species")
phy_r_species <- process_to_level(physeq1.rosenbergii, "Species")
phy_m_species <- process_to_level(physeq1.monodon,     "Species")

phy_v_phylum <- process_to_level(physeq1.vannamei,   "Phylum")
phy_r_phylum <- process_to_level(physeq1.rosenbergii, "Phylum")
phy_m_phylum <- process_to_level(physeq1.monodon,     "Phylum")



# --------------------------
# Helper to extract taxa names
# --------------------------
get_tax_names <- function(phy, level) {
  tt <- tax_table(phy)
  if (level %in% colnames(tt)) {
    as.character(tt[, level])
  } else {
    as.character(tt[, 1])
  }
}

# --------------------------
# Function to build a color palette and plot at any level
# --------------------------
make_palette <- function(phy_list, level) {
  taxa_all <- unique(unlist(lapply(phy_list, get_tax_names, level = level)))
  taxa_all <- taxa_all[!is.na(taxa_all)]
  taxa_all <- sort(taxa_all)
  
  n_taxa <- length(taxa_all)
  base_cols <- brewer.pal(min(max(3, n_taxa), 12), "Paired")
  if (n_taxa > length(base_cols)) {
    my_cols <- colorRampPalette(brewer.pal(12, "Paired"))(n_taxa)
  } else {
    my_cols <- base_cols[seq_len(n_taxa)]
  }
  names(my_cols) <- taxa_all
  return(my_cols)
}

plot_compositional <- function(physeq_obj, title_text, palette_named, level) {
  p <- microbiome::plot_composition(
    physeq_obj,
    otu.sort = "abundance", 
    average_by = "Gregarines"
  ) +
    ggtitle(title_text) +
    ylab("% Relative Abundance") +
    xlab("Samples") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 0)
    ) +
    scale_fill_manual(values = palette_named, na.value = "grey50", drop = FALSE) +
    guides(fill = guide_legend(title = level, ncol = 1))
  
  return(p)
}

# --------------------------
# Function to build composite plot (3 samples + shared legend)
# --------------------------
make_composite_plot <- function(phy_v, phy_r, phy_m, level_label) {
  palette_named <- make_palette(list(phy_v, phy_r, phy_m), level_label)
  
  # Create legend plot
  legend_df <- data.frame(
    Taxon = factor(names(palette_named), levels = names(palette_named)),
    x = 1, y = 1
  )
  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      name   = level_label,
      values = palette_named,
      breaks = names(palette_named),
      drop   = FALSE
    ) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 8)
    )
  legend_grob <- cowplot::get_legend(legend_plot)
  
  # Create 3 plots
  p_v <- plot_compositional(phy_v, "vannamei", palette_named, level_label) + theme(legend.position = "none")
  p_r <- plot_compositional(phy_r, "rosenbergii", palette_named, level_label) + theme(legend.position = "none")
  p_m <- plot_compositional(phy_m, "monodon", palette_named, level_label) + theme(legend.position = "none")
  
  # Combine
  left_row <- cowplot::plot_grid(p_v, p_r, p_m, ncol = 3, align = "v", axis = "lr")
  final_plot <- cowplot::plot_grid(left_row, legend_grob, rel_widths = c(3, 0.7))
  
  print(final_plot)
}

# --------------------------
# Generate Genus- and Species-level compositional plots
# --------------------------
make_composite_plot(phy_v_genus, phy_r_genus, phy_m_genus, "Genus")
make_composite_plot(phy_v_species, phy_r_species, phy_m_species, "Species")
make_composite_plot(phy_v_phylum,phy_r_phylum,phy_m_phylum, "Phylum")
make_composite_plot(phy_v_family, phy_r_family, phy_m_family, "Family")























































physeq_combined <- merge_phyloseq(physeq1.rosenbergii, physeq1.monodon, physeq1.vannamei)



# Extract relative abundance OTU table (already compositional from earlier step)
otu_table_pca <- as(otu_table(physeq_combined), "matrix")

# If taxa are rows, transpose:
if (taxa_are_rows(physeq_combined)) {
  otu_table_pca <- t(otu_table_pca)
}

# Run PCA
pca_result <- prcomp(otu_table_pca, center = TRUE, scale. = TRUE)

# Create a dataframe with PCA scores
pca_df <- as.data.frame(pca_result$x)
pca_df$SampleID <- rownames(pca_df)

# Prepare metadata
metadata <- data.frame(sample_data(physeq_combined)) %>%
  tibble::rownames_to_column("SampleID")

# PCA data
pca_df <- as.data.frame(pca_result$x) %>%
  tibble::rownames_to_column("SampleID")

# Merge while preserving all metadata columns
pca_df <- left_join(pca_df, metadata, by = "SampleID")

head(pca_df$Gregarines) 

ggplot(pca_df, aes(x = PC1, y = PC2, color = Blastocystis)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(
    title = "PCA - presences/ absense of greagrines ",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "% variance)")
  )


# Get PCA scores and add sample IDs
pca_df <- as.data.frame(pca_result$x) %>%
  tibble::rownames_to_column("SampleID")

# Prepare metadata (make sure SampleID is a column)
metadata <- data.frame(sample_data(physeq_combined)) %>%
  tibble::rownames_to_column("SampleID")

# Confirm Host column exists
print(colnames(metadata))

# Merge PCA data with metadata
pca_df <- left_join(pca_df, metadata, by = "SampleID")

# Check that Host is in the result
print(colnames(pca_df))
summary(pca_df$Organism)

# Plot PCA by Host - 1 plot combining all species 
ggplot(pca_df, aes(x = PC1, y = PC2, color = Host)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(
    title = "PCA of Microbial Communities by Host Species",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "% variance)")
  )


# Plot PCA by Host - 1 plot per species on same page
ggplot(pca_df, aes(x = PC1, y = PC2, color = Organism)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95, show.legend = FALSE) +
  facet_wrap(~ Organism) +
  theme_minimal() +
  labs(
    title = "PCA of Microbial Communities by Organism",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "% variance)")
  )





################################################################################
################################################################################
################################################################################
################################################################################


#Get table with all diversity indices in it
alpha_diversity <- estimate_richness(physeq_combined)
print(alpha_diversity)

#Concatenate desired columns of diversity table with subtype labels
alpha_diversity$Sample <- rownames(alpha_diversity)
selected_df <- alpha_diversity[, c("Sample", "Observed", "Chao1", "Simpson", "Shannon")]
rownames(selected_df) <- NULL #Gets rid of unnamed column at beginning
metadata$Sample <- rownames(metadata) #Creats a new column named 'Sample'
samples <- metadata$SampleID
organisms <- metadata$Blastocystis
shannonscore <- alpha_diversity$Shannon
simpsonscore <- alpha_diversity$Simpson
chaoscore <- alpha_diversity$Chao1
observed <- alpha_diversity$Observed

statframe <- data.frame(samples, organisms, shannonscore, simpsonscore, chaoscore, observed)

head(statframe [1:5, 1:5]) 




#Check normalisation of data (p>0.05 = normal, p<0.05 = non-normal)
shapiro.shan <- shapiro_test(alpha_diversity$Shannon)$p.value
shapiro.chao <- shapiro_test(alpha_diversity$Chao1)$p.value
shapiro.simp <- shapiro_test(alpha_diversity$Simpson)$p.value
shapiro.observed <- shapiro_test(alpha_diversity$Observed)$p.value

#Can check visually via a histogram
histogram(alpha_diversity$Shannon)

#Will's script for STATISTICS
#Stats tests for subtypes
if(shapiro.shan<0.05){result.kruskal.shannon<-statframe %>% kruskal_test(shannonscore ~ organisms)}
if(shapiro.shan>0.05){result.anova.shannon<-statframe %>% anova_test(shannonscore ~ organisms)}
if(shapiro.simp<0.05){result.kruskal.simpson<-statframe %>% kruskal_test(simpsonscore ~ organisms)}
if(shapiro.simp>0.05){result.anova.simpson<-statframe %>% anova_test(simpsonscore ~ organisms)}
if(shapiro.chao<0.05){result.kruskal.chao<-statframe %>% kruskal_test(chaoscore ~ organisms)}
if(shapiro.chao>0.05){result.anova.chao<-statframe %>% anova_test(chanoscore ~ organisms)}
if(shapiro.observed<0.05){result.kruskal.observed<-statframe %>% kruskal_test(observed ~ organisms)}
if(shapiro.observed>0.05){result.anova.observed<-statframe %>% anova_test(observed ~ organisms)}

#Pairwise comparisons for Organisms
if(exists("result.kruskal.shannon")){pwc.shannon<- statframe %>% dunn_test(shannonscore ~ organisms, p.adjust.method = "bonferroni")}
if(exists("result.anova.shannon")){pwc.shannon<- statframe %>% tukey_hsd(shannonscore ~ organisms)}
if(exists("result.kruskal.simpson")){pwc.simpson<- statframe %>% dunn_test(simpsonscore ~ organisms, p.adjust.method = "bonferroni")}
if(exists("result.anova.simpson")){pwc.simpson<- statframe %>% tukey_hsd(simpsonscore ~ organisms)}
if(exists("result.kruskal.chao")){pwc.chao<- statframe %>% dunn_test(chaoscore ~ organisms, p.adjust.method = "bonferroni")}
if(exists("result.anova.chao")){pwc.chao<- statframe %>% tukey_hsd(chaoscore ~ organisms)}
if(exists("result.kruskal.observed")){pwc.observed<- statframe %>% dunn_test(observed ~ organisms, p.adjust.method = "bonferroni")}
if(exists("result.anova.observed")){pwc.observed<- statframe %>% tukey_hsd(observed~ organisms)}


#################SHANNON##################################

#Plot violin plots with pairwise comparisons
pwc.shannon <- pwc.shannon %>% add_xy_position(x = "organisms")
statplot.shannon<-ggviolin(
  statframe,
  x = "organisms",
  y = "shannonscore",
  title = "Shannon Diversity of Species",
  fill = "organisms",
  palette = "Paired" #Paired color palette from RColorBrewer
) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = ("right" )
  ) +
  stat_pvalue_manual(
    pwc.shannon,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(
    caption = get_pwc_label(pwc.shannon)
  )
if(exists("result.kruskal.shannon")){
  statplot.shannon +
    labs(subtitle = get_test_label(result.kruskal.shannon, detailed = FALSE))
}
if(exists("result.anova.shannon")){
  statplot.shannon +
    labs(subtitle = get_test_label(result.anova.shannon, detailed = FALSE))
}


###############JITTER

pwc.shannon <- pwc.shannon %>% add_xy_position(x = "organisms")

statplot.shannon <- ggviolin(
  statframe,
  x = "organisms",
  y = "shannonscore",
  title = "Shannon Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  geom_jitter(aes(x = organisms, y = shannonscore),
              size = 1.5, width = 0.2, alpha = 0.6) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  stat_pvalue_manual(
    pwc.shannon,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(
    caption = get_pwc_label(pwc.shannon)
  )

if(exists("result.kruskal.shannon")){
  statplot.shannon + labs(subtitle = get_test_label(result.kruskal.shannon))
}
if(exists("result.anova.shannon")){
  statplot.shannon + labs(subtitle = get_test_label(result.anova.shannon))
}


#################SIMPSON##################################

#Plot violin plots with pairwise comparisons
pwc.simpson <- pwc.simpson %>% add_xy_position(x = "organisms")
statplot.simpson<-ggviolin(
  statframe,
  x = "organisms",
  y = "simpsonscore",
  title = "Simpson Diversity of Species",
  fill = "organisms",
  palette = "Paired" #Paired color palette from RColorBrewer
) +
  geom_jitter(aes(x = organisms, y = simpsonscore),
              size = 1.5, width = 0.2, alpha = 0.6) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = ("right" )
  ) +
  stat_pvalue_manual(
    pwc.simpson,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(
    caption = get_pwc_label(pwc.simpson)
  )
if(exists("result.kruskal.simpson")){
  statplot.simpson +
    labs(subtitle = get_test_label(result.kruskal.simpson, detailed = FALSE))
}
if(exists("result.anova.simpson")){
  statplot.simpson +
    labs(subtitle = get_test_label(result.anova.simpson, detailed = FALSE))
}

#################CHAO##################################

#Plot violin plots with pairwise comparisons
pwc.chao <- pwc.chao %>% add_xy_position(x = "organisms")
statplot.chao<-ggviolin(
  statframe,
  x = "organisms",
  y = "chaoscore",
  title = "Chao Diversity of Species",
  fill = "organisms",
  palette = "Paired" #Paired color palette from RColorBrewer
) +
  geom_jitter(aes(x = organisms, y = chaoscore),
              size = 1.5, width = 0.2, alpha = 0.6) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = ("right" )
  ) +
  stat_pvalue_manual(
    pwc.chao,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(
    caption = get_pwc_label(pwc.chao)
  )
if(exists("result.kruskal.chao")){
  statplot.chao +
    labs(subtitle = get_test_label(result.kruskal.chao, detailed = FALSE))
}
if(exists("result.anova.chao")){
  statplot.chao +
    labs(subtitle = get_test_label(result.anova.chao, detailed = FALSE))
}


#################OBSERVED##################################


#Plot violin plots with pairwise comparisons
pwc.observed <- pwc.observed %>% add_xy_position(x = "organisms")
statplot.observed<-ggviolin(
  statframe,
  x = "organisms",
  y = "observed",
  title = "Observed Diversity of Species",
  fill = "organisms",
  palette = "Paired" #Paired color palette from RColorBrewer
) +
  geom_jitter(aes(x = organisms, y = observed),
              size = 1.5, width = 0.2, alpha = 0.6) + 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = ("right" )
  ) +
  stat_pvalue_manual(
    pwc.observed,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(
    caption = get_pwc_label(pwc.chao)
  )
if(exists("result.kruskal.observed")){
  statplot.observed +
    labs(subtitle = get_test_label(result.kruskal.observed, detailed = FALSE))
}
if(exists("result.anova.observed")){
  statplot.observed +
    labs(subtitle = get_test_label(result.anova.observed, detailed = FALSE))
}


# To combine all four on a 2 by 2 grid
library("patchwork")

(statplot.shannon | statplot.simpson) / (statplot.chao | statplot.observed) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")



################WITHOUT PAIRWISE COMPARISON###############
##########################################################
##########################################################


#################SHANNON##################################

#Plot violin plots without pairwise comparisons
shannon <- pwc.shannon %>% add_xy_position(x = "orgnaisms")
statplot.shannon <- ggviolin(
  statframe,
  x = "organisms",
  y = "shannonscore",
  title = "Shannon Diversity of Species",
  fill = "organisms",
  palette = "Paired" #Paired palette from RColorBrewer
) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text = element_text(size=0),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )
if(exists("result.kruskal.shannon")){
  statplot.shannon <- statplot.shannon +
    labs(subtitle = get_test_label(result.kruskal.shannon, detailed = FALSE))
}
if(exists("result.anova.shannon")){
  statplot.shannon <- statplot.shannon+
    labs(subtitle = get_test_label(result.anova.shannon, detailed = FALSE))
}
print(statplot.shannon)

#################SIMPSON##################################

simpson <- pwc.simpson %>% add_xy_position(x = "organisms")
statplot.simpson <- ggviolin(statframe, 
                             x = "organisms",
                             y = "simpsonscore",
                             title = "Simpson Diversity of organisms", 
                             fill = "organisms", palette = "Paired") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.position = "right")

if(exists("result.kruskal.simpson")){
  statplot.simpson <- statplot.simpson + 
    labs(subtitle = get_test_label(result.kruskal.simpson, detailed = FALSE))}

if(exists("result.anova.simpson")){
  statplot.simpson <- statplot.simpson + 
    labs(subtitle = get_test_label(result.anova.simpson, detailed = FALSE))}
print(statplot.simpson)

#################CHAO##################################

chao <- pwc.chao %>% add_xy_position(x = "organisms")
statplot.chao <- ggviolin(statframe, 
                          x = "organisms",
                          y = "chaoscore",
                          title = "Chao Diversity of organisms", 
                          fill = "organisms", palette = "Paired") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.position = "right")

if(exists("result.kruskal.chao")){
  statplot.chao <- statplot.chao + 
    labs(subtitle = get_test_label(result.kruskal.chao, detailed = FALSE))}

if(exists("result.anova.chao")){
  statplot.chao <- statplot.chao + 
    labs(subtitle = get_test_label(result.anova.chao, detailed = FALSE))}
print(statplot.chao)



#################OBSERVED##################################

observed <- pwc.observed %>% add_xy_position(x = "organisms") 
statplot.observed <- ggviolin(statframe,
                              x = "organisms",
                              y = "observed", 
                              title = "Observed Taxa",
                              fill = "organisms", 
                              palette = "Paired") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 0),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"), legend.position = "right")
if (exists("result.kruskal.observed")) {
  statplot.observed <- statplot.observed + 
    labs(subtitle = get_test_label(result.kruskal.observed, detailed = FALSE))}

if (exists("result.anova.observed")) {
  statplot.observed <- statplot.observed + 
    labs(subtitle = get_test_label(result.anova.observed, detailed = FALSE))}
print(statplot.observed)



#############COMBINE##########################################
library("patchwork")

(statplot.shannon | statplot.simpson) / (statplot.chao | statplot.observed) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")


################################################################################
################################################################################
###############################Beta Diversity ##################################
################################# /PCoA: #######################################
################################################################################



########################Just Bray###############################################
ordination <- ordinate(physeq_combined, method = "PCoA", distance = "bray")
dist_matrix <- phyloseq::distance(physeq_combined, method = "bray")

meta<- data.frame(sample_data(physeq_combined))
adonis_host <- adonis2(dist_matrix ~ Blastocystis, data = meta, permutations = 999)

pval_host <- adonis_host$`Pr(>F)`[1]
r2_host   <- adonis_host$R2[1]
label_host <- paste0("PERMANOVA: R² = ", round(r2_host, 3), ", p = ", signif(pval_host, 3))


p1 <- plot_ordination(physeq_combined, ordination, color = "Blastocystis") +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Beta diversity by Host", 
       subtitle = label_host) +
  annotate("text", x = 0.7, y = 0.7, label = label_host, hjust = 0)

adonis_org <- adonis2(dist_matrix ~ Organism, data = meta, permutations = 999)
pval_org <- adonis_org$`Pr(>F)`[1]
r2_org   <- adonis_org$R2[1]
label_org <- paste0("PERMANOVA: R² = ", round(r2_org, 3), ", p = ", signif(pval_org, 3))

p2 <- plot_ordination(physeq_combined, ordination, color = "Organism") +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Beta diversity by Organism", 
       subtitle = label_host) +
  annotate("text", x = 0.7, y = 0.7, label = label_org, hjust = 0)

p1
p2




########################Combined methods########################################


#Settings
group_var <- "Blastocystis"  # change to "Organism" or another metadata field if you prefer
methods <- c("bray", "jaccard_pa", "jaccard_ru", "canberra")

# Helper: get OTU matrix with samples as rows
get_otu_mat <- function(ps) {
  x <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) x <- t(x)
  x
}

# Compute distances (note: jaccard_ru = quantitative Jaccard = Ružička)
compute_distance <- function(ps, method) {
  otu <- get_otu_mat(ps)
  if (method == "jaccard_pa") {
    vegan::vegdist(otu, method = "jaccard", binary = TRUE)
  } else if (method == "jaccard_ru") {
    vegan::vegdist(otu, method = "jaccard", binary = FALSE)  # Ružička
  } else {
    phyloseq::distance(ps, method = method)                  # bray, canberra, etc.
  }
}

# Make one PCoA + PERMANOVA plot for a method
make_pcoa_plot <- function(ps, method, group_var) {
  dist_mat <- compute_distance(ps, method)
  # Subset metadata to exactly the samples in the distance object
  labs <- attr(dist_mat, "Labels")
  meta <- data.frame(sample_data(ps))[labs, , drop = FALSE]
  
  # Drop samples with NA in the grouping variable (avoid adonis2 NA error)
  keep <- stats::complete.cases(meta[[group_var]])
  if (!all(keep)) {
    dist_mat <- as.dist(as.matrix(dist_mat)[keep, keep])
    meta <- meta[keep, , drop = FALSE]
  }
  
  # Ordination
  ord <- ordinate(ps, method = "PCoA", distance = dist_mat)
  
  # PERMANOVA
  form <- stats::as.formula(paste("dist_mat ~", group_var))
  ad <- vegan::adonis2(form, data = meta, permutations = 999)
  r2 <- round(ad$R2[1], 3)
  pv <- signif(ad$`Pr(>F)`[1], 3)
  
  # Nice name
  title_name <- switch(method,
                       "bray"        = "Bray–Curtis",
                       "jaccard_pa"  = "Jaccard (presence/absence)",
                       "jaccard_ru"  = "Jaccard (Ružička, abundance)",
                       "canberra"    = "Canberra",
                       toupper(method)
  )
  
  # Plot with subtitle = PERMANOVA
  p <- plot_ordination(ps, ord, color = group_var) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste("PCoA –", title_name),
      subtitle = paste0("PERMANOVA: R² = ", r2, ", p = ", pv)
    ) +
    theme(legend.position = "right")
  
  list(plot = p, dist = dist_mat)
}

#Build all plots
plots <- list()
dists  <- list()

for (m in methods) {
  res <- make_pcoa_plot(physeq_combined, m, group_var)
  plots[[m]] <- res$plot
  dists[[m]] <- res$dist
}

# Combine with a single legend on the right
combined_plot <- wrap_plots(plots, ncol = 2, guides = "collect") &
  theme(legend.position = "right")
print(combined_plot)

# Compare distance matrices: Mantel tests (all pairs)
meth_names <- names(dists)
mantel_r <- matrix(NA_real_, length(meth_names), length(meth_names),
                   dimnames = list(meth_names, meth_names))
mantel_p <- mantel_r

for (i in seq_along(meth_names)) {
  for (j in seq_along(meth_names)) {
    if (j > i) {
      mt <- vegan::mantel(dists[[meth_names[i]]], dists[[meth_names[j]]], permutations = 999)
      mantel_r[i, j] <- unname(mt$statistic)
      mantel_p[i, j] <- unname(mt$signif)
    } else if (j == i) {
      mantel_r[i, j] <- 1
      mantel_p[i, j] <- 0
    }
  }
}
cat("\nMantel r (upper triangle) and p (lower triangle):\n")
print(mantel_r)
print(mantel_p)

#Procrustes vs Bray baseline (using cmdscale coordinates)
coords <- lapply(dists, function(d) stats::cmdscale(d, k = 2))
baseline <- "bray"
proc_tab <- data.frame(
  compare_to_bray = setdiff(names(coords), baseline),
  r = NA_real_,
  p = NA_real_
)
k <- 1
for (m in setdiff(names(coords), baseline)) {
  # Align row order
  common <- intersect(rownames(coords[[baseline]]), rownames(coords[[m]]))
  X <- coords[[baseline]][common, , drop = FALSE]
  Y <- coords[[m]][common, , drop = FALSE]
  pr <- vegan::protest(X, Y, permutations = 999)
  proc_tab$r[k] <- unname(pr$t0)       # correlation-like statistic
  proc_tab$p[k] <- unname(pr$signif)   # p-value
  k <- k + 1
}
cat("\nProcrustes (vs Bray) results:\n")
print(proc_tab)















methods <- c("bray", "jaccard", "canberra", "ruzicka")
plots <- list()
dist_matrices <- list()

for (method in methods) {
  # Compute distance
  if (method == "ruzicka") {
    otu <- as(otu_table(physeq_combined), "matrix")
    if (taxa_are_rows(physeq_combined)) {
      otu <- t(otu)
    }
    dist_matrix <-vegan::vegdist(otu, method = "ruzicka")
  } else {
    dist_matrix <- phyloseq::distance(physeq_combined, method = method)
    dist_matrices[[method]] <- dist_matrix }
  
  # Run PCoA
  ordination <- ordinate(physeq_combined, method = "PCoA", distance = dist_matrix)
  
  # Metadata
  meta <- as(sample_data(physeq_combined), "data.frame")
  
  # PERMANOVA
  adonis_res <- adonis2(dist_matrix ~ Host, data = meta, permutations = 999)
  r2 <- round(adonis_res$R2[1], 3)
  pval <- signif(adonis_res$`Pr(>F)`[1], 3)
  
  # Plot
  p <- plot_ordination(physeq_combined, ordination, color = "Host") +
    geom_point(size = 2, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = paste("PCoA -", toupper(method), "distance"),
      subtitle = paste0(method, " PERMANOVA: R² = ", r2, ", p = ", pval)
    ) +
    theme(legend.position = "right")
  
  plots[[method]] <- p
}

# Combine into one figure with shared legend
combined_plot <- wrap_plots(plots, ncol = 2, guides = "collect") & theme(legend.position = "right")
print(combined_plot)


#Compare distance matrices
# ----- Mantel test between each pair
mantel_results <- matrix(NA, nrow = length(methods), ncol = length(methods),
                         dimnames = list(methods, methods))

for (i in 1:(length(methods)-1)) {
  for (j in (i+1):length(methods)) {
    m <- mantel(dist_matrices[[methods[i]]], dist_matrices[[methods[j]]], permutations = 999)
    mantel_results[i, j] <- m$statistic
    mantel_results[j, i] <- m$signif
  }
}

mantel_results


# ---- other plots:
make_pcoa_plot <- function(physeq, method, group_var = "Host") {
  ord <- ordinate(physeq, method = "PCoA", distance = method)
  dist_matrix <- phyloseq::distance(physeq, method = method)
  meta <- data.frame(sample_data(physeq_combined))
  form <- as.formula(paste("dist_matrix ~", group_var))
  adonis_res <- adonis2(form, data = meta, permutations = 999)
  r2 <- adonis_res$R2[1]
  pval <- adonis_res$`Pr(>F)`[1]
  label <- paste0(method, " PERMANOVA: R² = ", round(r2, 3), ", p = ", signif(pval, 3))
  p <- plot_ordination(physeq, ord, color = group_var) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste("PCoA -", toupper(method), "distance"),
      subtitle = label
    )
  return(p)
}


plots <- lapply(dist_methods, function(m) make_pcoa_plot(physeq_combined, m, "Host"))


combined_plot <- wrap_plots(plots, ncol = 2, guides = "collect", widths = c(2, 1)) & 
  theme(legend.position = "right")

combined_plot



################################################################################
################################################################################
################################################################################

# Robust heatmap script for phyloseq -> pheatmap
library(phyloseq)
library(pheatmap)
library(dplyr)

# Parameters
physeq_obj <- physeq_combined   # replace if needed
taxrank_to_use <- "Phylum"
top_n <- 30

# 1) relative abundances and tax_glom
ps_rel <- transform_sample_counts(physeq_obj, function(x) if(sum(x) == 0) x else x / sum(x))
ps_tax <- tax_glom(ps_rel, taxrank = taxrank_to_use)

# 2) extract abundance table and taxonomy as data.frames
abund_df <- as.data.frame(otu_table(ps_tax))
# Ensure rows = taxa, cols = samples. If taxa_are_rows = FALSE, transpose:
if (!taxa_are_rows(ps_tax)) {
  abund_df <- t(abund_df)
  abund_df <- as.data.frame(abund_df)
}
tax_df <- as.data.frame(tax_table(ps_tax))

# 3) build safe unique taxon labels
# Use the taxonomic level requested; fallback to Family, then to "Unclassified"
tax_labels <- as.character(tax_df[[taxrank_to_use]])
missing_lab <- is.na(tax_labels) | tax_labels == ""
if (any(missing_lab)) {
  fallback <- if ("Family" %in% colnames(tax_df)) as.character(tax_df$Family) else rep("Unclassified", nrow(tax_df))
  tax_labels[missing_lab] <- fallback[missing_lab]
}
tax_labels_unique <- make.unique(tax_labels)
rownames(abund_df) <- tax_labels_unique

## 4) choose top N taxa by rowSums (handle case with fewer taxa)
#n_taxa_available <- nrow(abund_df)
#n_to_take <- min(top_n, n_taxa_available)
#if (n_taxa_available == 0) stop("No taxa left after aggregation - check your phyloseq object.")
#top_taxa <- names(sort(rowSums(abund_df), decreasing = TRUE))[1:n_to_take]
#abund_top <- abund_df[top_taxa, , drop = FALSE]

# 4) filter taxa by mean relative abundance >= 0.5% (i.e., 0.005)
mean_abund <- rowMeans(abund_df)
keep_taxa <- names(mean_abund[mean_abund >= 0.005])

if (length(keep_taxa) == 0) stop("No taxa have mean relative abundance >= 0.5% - try lowering the threshold.")
abund_top <- abund_df[keep_taxa, , drop = FALSE]

# 5) convert to numeric matrix (pheatmap likes numeric matrices)
abund_mat <- as.matrix(abund_top)
mode(abund_mat) <- "numeric"

# Optional transformation for visualization (log)
abund_mat_plot <- log10(abund_mat + 1e-5)




# 5) convert to numeric matrix (pheatmap likes numeric matrices)
abund_mat <- as.matrix(abund_top)
mode(abund_mat) <- "numeric"

# Optional transformation for visualization (log)
abund_mat_plot <- log10(abund_mat + 1e-5)

# 6) Build sample annotation and align to matrix columns
# Create annotation from the same phyloseq object used above (ps_tax)
samp_df <- as.data.frame(sample_data(ps_tax))
# Ensure sample names are rownames of samp_df
if (!all(rownames(samp_df) == sample_names(ps_tax))) {
  rownames(samp_df) <- sample_names(ps_tax)
}

# Build annotation frame with the metadata column "Species" (adjust name if your column is different)
meta_colname <- "Organism"
if (!meta_colname %in% colnames(samp_df)) {
  stop(paste0("Metadata column '", meta_colname, "' not found in sample_data(ps_tax)."))
}
annotation_full <- data.frame(Species = as.character(samp_df[[meta_colname]]))
rownames(annotation_full) <- rownames(samp_df)

# Now align annotation rows to columns of the abundance matrix
matrix_samples <- colnames(abund_mat_plot)
# Find intersection and warn if anything missing
matched_idx <- match(matrix_samples, rownames(annotation_full))
if (any(is.na(matched_idx))) {
  missing_samples <- matrix_samples[is.na(matched_idx)]
  warning("The following samples are present in the abundance matrix but missing in sample_data and will be removed: ",
          paste(missing_samples, collapse = ", "))
  # drop missing samples from matrix
  keep_cols <- !colnames(abund_mat_plot) %in% missing_samples
  abund_mat_plot <- abund_mat_plot[, keep_cols, drop = FALSE]
  matrix_samples <- colnames(abund_mat_plot)
  matched_idx <- matched_idx[!is.na(matched_idx)]
  if (ncol(abund_mat_plot) == 0) stop("No samples remain after removing unmatched samples.")
}

# Subset the annotation to the exact ordering of the matrix columns
annotation_col <- annotation_full[matrix_samples, , drop = FALSE]
rownames(annotation_col) <- matrix_samples

# Final sanity checks
if (!identical(colnames(abund_mat_plot), rownames(annotation_col))) {
  stop("Final check failed: column names of abundance matrix and rownames of annotation do not match.")
}

# 7) Plot heatmap
pheatmap(abund_mat_plot,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize = 10,
         main = "Taxa with mean relative abundance across all samples ≥ 0.5%" )




################################################################################
################################################################################
################################################################################


# Input phyloseq object
ps <- physeq_combined   # merged phyloseq

# Convert phyloseq -> DESeq2 (use the correct metadata column 'Host')
dds <- phyloseq_to_deseq2(ps, ~ Host)
dds$Host <- factor(dds$Host)
dds <- DESeq(dds)

# Prepare taxonomy table with a safe label for plotting
tax_df <- as.data.frame(tax_table(ps))
tax_df$taxa <- rownames(tax_df)
# Create a label column: prefer Genus; otherwise Family_Genus or "Unclassified"
if (!"Genus" %in% colnames(tax_df)) tax_df$Genus <- NA
if (!"Family" %in% colnames(tax_df)) tax_df$Family <- NA
tax_df$TaxLabel <- ifelse(!is.na(tax_df$Genus) & tax_df$Genus != "",
                          as.character(tax_df$Genus),
                          paste0(ifelse(is.na(tax_df$Family) | tax_df$Family == "", "Unclassified", as.character(tax_df$Family)),
                                 "_", seq_len(nrow(tax_df))))
# Ensure unique and readable
tax_df$TaxLabel <- make.unique(as.character(tax_df$TaxLabel))

# Create output directory
dir.create("DESeq2_VolcanoPlots", showWarnings = FALSE)

# All pairwise comparisons of Host levels
comparisons <- combn(levels(dds$Host), 2, simplify = FALSE)

for (comp in comparisons) {
  A <- comp[1]; B <- comp[2]
  message("Running comparison: ", A, " vs ", B)
  
  # Get results
  res <- results(dds, contrast = c("Host", A, B))
  res <- as.data.frame(res)
  res$taxa <- rownames(res)
  
  # Merge taxonomy
  res_tax <- left_join(res, tax_df, by = "taxa")
  
  # Replace NA padj with 1 (so -log10(padj) doesn't explode)
  res_tax$padj[is.na(res_tax$padj)] <- 1
  
  # Create Significance label
  res_tax <- res_tax %>%
    mutate(Significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ paste("Up in", A),
      padj < 0.05 & log2FoldChange < -1 ~ paste("Up in", B),
      TRUE ~ "Not significant"
    ))
  
  # Create dynamic color vector safely
  colors_vec <- c("firebrick3", "steelblue3", "grey70")
  names(colors_vec) <- c(paste("Up in", A), paste("Up in", B), "Not significant")
  
  # Choose top taxa for labeling (best padj)
  top_taxa <- res_tax %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    slice_head(n = 10)
  
  # Make the volcano plot
  p <- ggplot(res_tax, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.8, size = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top_taxa, aes(label = TaxLabel), size = 3, max.overlaps = 15) +
    scale_color_manual(values = colors_vec) +
    theme_bw(base_size = 14) +
    labs(title = paste0("Volcano: ", A, " vs ", B),
         x = paste0("Log2 Fold Change (", A, " vs ", B, ")"),
         y = "-log10(adjusted p-value)",
         color = "Enrichment") +
    theme(legend.position = "top")
  
  # Print and save
  print(p)
  ggsave(filename = paste0("DESeq2_VolcanoPlots/Volcano_", gsub(" ", "_", A), "_vs_", gsub(" ", "_", B), ".png"),
         plot = p, width = 8, height = 6, dpi = 300)
  
  # Save result table
  write.csv(res_tax, file = paste0("DESeq2_VolcanoPlots/DESeq2_results_", gsub(" ", "_", A), "_vs_", gsub(" ", "_", B), ".csv"),
            row.names = FALSE)
}

############################################################
#CORE MICROBIOME AT 80% prevalnce in each sample ###
############################################################



#############################Generic code for all data here: ###################
physeq_comp<-physeq_combined %>% microbiome::transform("compositional")
#use aggregate_taxa to make it genus
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
prevalences <- seq(.05, 1, .05)
physeq_comp_ge<-aggregate_taxa(physeq_comp, level = "Phylum")
taxa_names(physeq_comp_ge)
#then create a colour palette e.g.
core_cols <- c("yellow","blue")
#then insert both into plot_core
core1 <- plot_core(physeq_comp_ge,
                   plot.type = "heatmap",
                   prevalences = prevalences, colours = core_cols,
                   detections = detections, 
                   min.prevalence = 0.8) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") 
#the 0.5 min prevalence means any taxa that appears in >50% of

core1

#####Refined by species:########################################################

phy_rosen <- physeq1.rosenbergii %>% microbiome::transform("compositional")
phy_mono  <- physeq1.monodon %>% microbiome::transform("compositional")
phy_van   <- physeq1.vannamei %>% microbiome::transform("compositional")

# Aggregate to phylum (or genus)
phy_rosen_ge <- aggregate_taxa(phy_rosen, "Genus")
phy_mono_ge  <- aggregate_taxa(phy_mono,  "Genus")
phy_van_ge   <- aggregate_taxa(phy_van,   "Genus")

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)
prevalences <- seq(.05, 1, .05)

core_cols <- c("yellow", "blue")

core_rosen <- plot_core(phy_rosen_ge,
                        plot.type="heatmap",
                        prevalences=prevalences,
                        detections=detections,
                        colours=core_cols,
                        min.prevalence=0.5) +
  ggtitle("Rosenbergii")

core_mono <- plot_core(phy_mono_ge,
                       plot.type="heatmap",
                       prevalences=prevalences,
                       detections=detections,
                       colours=core_cols,
                       min.prevalence=0.5) +
  ggtitle("Monodon")

core_van <- plot_core(phy_van_ge,
                      plot.type="heatmap",
                      prevalences=prevalences,
                      detections=detections,
                      colours=core_cols,
                      min.prevalence=0.5) +
  ggtitle("Vannamei")



combined_core_plot <- core_rosen + core_mono + core_van +
  plot_layout(ncol = 1)

combined_core_plot



############################################################
# Co-occurrence Network of Shrimp Microbiomes
############################################################

# =======================================================
# CO-OCCURRENCE NETWORK (GENUS LEVEL, DUPLICATES HANDLED)
# =======================================================

library(phyloseq)
library(igraph)
library(ggraph)
library(tidyverse)
library(Hmisc)
library(ggrepel)

# 1️⃣ Relative abundance transformation
ps_rel <- transform_sample_counts(physeq_combined, function(x) x / sum(x))

# 2️⃣ Aggregate to genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# 3️⃣ Extract taxonomy and abundance data
tax_df <- as.data.frame(tax_table(ps_genus))
otu_df <- as.data.frame(otu_table(ps_genus))

# Make sure taxa are rows
if (!taxa_are_rows(ps_genus)) {
  otu_df <- t(otu_df)
}


# 4️⃣ Clean up Genus names
tax_df$Genus <- as.character(tax_df$Genus)
tax_df$Genus[is.na(tax_df$Genus) | tax_df$Genus == "" | tax_df$Genus == "?"] <- NA

# Make Genus names unique to avoid rowname clash (temporarily)
tax_df$Genus_unique <- make.unique(ifelse(is.na(tax_df$Genus), "Unclassified", tax_df$Genus))

# Apply those as row names (no duplicates now)
rownames(otu_df) <- tax_df$Genus_unique

# 5️⃣ Remove unclassified taxa
otu_df <- otu_df[!grepl("^Unclassified", rownames(otu_df)), , drop = FALSE]

# 6️⃣ Collapse duplicate genera (sum across samples)
otu_df <- otu_df %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Genus")

# 7️⃣ Filter by mean relative abundance (≥0.2%)
mean_abund <- rowMeans(otu_df)
otu_df_filt <- otu_df[mean_abund >= 0.0005, , drop = FALSE]

# 8️⃣ Compute Spearman correlations
abund_mat <- t(otu_df_filt)
storage.mode(abund_mat) <- "double"
cor_res <- Hmisc::rcorr(as.matrix(abund_mat), type = "spearman")
cor_mat <- cor_res$r
p_mat <- cor_res$P

# 9️⃣ Keep strong & significant correlations
r_cutoff <- 0.5
sig_mask <- (abs(cor_mat) >= r_cutoff) & (p_mat < 0.05)
cor_mat[!sig_mask] <- 0

# 🔟 Build network
net <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
V(net)$Genus <- V(net)$name

# 1️⃣1️⃣ Determine shared vs host-specific genera
sample_df <- as.data.frame(sample_data(physeq_combined))
otu_bin <- as.matrix(otu_df_filt) > 0

# Compute presence across hosts
taxa_pres <- apply(otu_bin, 1, function(x) tapply(x, sample_df$Host, any))
if (is.list(taxa_pres)) taxa_pres <- do.call(cbind, taxa_pres)
shared_taxa <- names(which(rowSums(taxa_pres, na.rm = TRUE) > 1))

V(net)$Shared <- ifelse(V(net)$Genus %in% shared_taxa, "Shared", "Host-specific")

# 1️⃣2️⃣ Plot
ggraph(net, layout = "fr") +
  geom_edge_link(aes(alpha = abs(weight)), color = "grey80") +
  geom_node_point(aes(color = Shared, size = degree(net))) +
  geom_node_text(aes(label = Genus),
                 repel = TRUE,
                 size = 3,
                 check_overlap = TRUE) +
  scale_color_manual(values = c("Shared" = "orange", "Host-specific" = "magenta")) +
  scale_size(range = c(2, 8)) +
  theme_void(base_size = 14) +
  ggtitle("Genus-level Co-occurrence Network (≥0.05% abundance)")




# =======================================================
# CO-OCCURRENCE NETWORK (GENUS OR HIGHER TAXA LABELS)
# =======================================================

library(phyloseq)
library(igraph)
library(ggraph)
library(tidyverse)
library(Hmisc)
library(ggrepel)

# 1️⃣ Relative abundance transformation
ps_rel <- transform_sample_counts(physeq_combined, function(x) x / sum(x))

# 2️⃣ Aggregate to genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# 3️⃣ Extract taxonomy and abundance tables
tax_df <- as.data.frame(tax_table(ps_genus))
otu_df <- as.data.frame(otu_table(ps_genus))

# Make sure taxa are rows
if (!taxa_are_rows(ps_genus)) {
  otu_df <- t(otu_df)
}

# 4️⃣ Fill in missing genus names using higher-level taxonomy
tax_df$Genus <- as.character(tax_df$Genus)
tax_df$Family <- as.character(tax_df$Family)
tax_df$Order <- as.character(tax_df$Order)
tax_df$Phylum <- as.character(tax_df$Phylum)

# Replace missing or "?" genera with higher-level taxon + "_?"
tax_df$Genus_fixed <- ifelse(
  is.na(tax_df$Genus) | tax_df$Genus == "" | tax_df$Genus == "?",
  ifelse(!is.na(tax_df$Family) & tax_df$Family != "", paste0(tax_df$Family, "_?"),
         ifelse(!is.na(tax_df$Order) & tax_df$Order != "", paste0(tax_df$Order, "_?"),
                ifelse(!is.na(tax_df$Phylum) & tax_df$Phylum != "", paste0(tax_df$Phylum, "_?"),
                       "Unclassified_?"))),
  tax_df$Genus
)

# Make names unique
tax_df$Genus_unique <- make.unique(tax_df$Genus_fixed)

# Apply genus/higher-level names as row names
rownames(otu_df) <- tax_df$Genus_unique

# 5️⃣ Combine duplicates (sum across ASVs)
otu_df <- otu_df %>%
  as.data.frame() %>%
  rownames_to_column("Taxon") %>%
  group_by(Taxon) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Taxon")

# 6️⃣ Filter taxa by mean relative abundance (≥0.2%)
mean_abund <- rowMeans(otu_df)
otu_df_filt <- otu_df[mean_abund >= 0.002, , drop = FALSE]

# 7️⃣ Compute correlations (Spearman)
abund_mat <- t(otu_df_filt)
storage.mode(abund_mat) <- "double"
cor_res <- Hmisc::rcorr(as.matrix(abund_mat), type = "spearman")

cor_mat <- cor_res$r
p_mat <- cor_res$P

# 8️⃣ Keep strong, significant correlations
r_cutoff <- 0.6
sig_mask <- (abs(cor_mat) >= r_cutoff) & (p_mat < 0.05)
cor_mat[!sig_mask] <- 0

# 9️⃣ Build network
net <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
V(net)$Taxon <- V(net)$name

# 🔟 Determine shared vs host-specific taxa
sample_df <- as.data.frame(sample_data(physeq_combined))
otu_bin <- as.matrix(otu_df_filt) > 0

taxa_pres <- apply(otu_bin, 1, function(x) tapply(x, sample_df$Host, any))
if (is.list(taxa_pres)) taxa_pres <- do.call(cbind, taxa_pres)
shared_taxa <- names(which(rowSums(taxa_pres, na.rm = TRUE) > 1))

V(net)$Shared <- ifelse(V(net)$Taxon %in% shared_taxa, "Shared", "Host-specific")

# 1️⃣1️⃣ Plot
ggraph(net, layout = "fr") +
  geom_edge_link(aes(alpha = abs(weight)), color = "grey80") +
  geom_node_point(aes(color = Shared, size = degree(net))) +
  geom_node_text(aes(label = Taxon),
                 repel = TRUE,
                 size = 3,
                 check_overlap = TRUE) +
  scale_color_manual(values = c("Shared" = "orange", "Host-specific" = "steelblue")) +
  scale_size(range = c(2, 8)) +
  theme_void(base_size = 14) +
  ggtitle("Co-occurrence Network (≥0.2% abundance, genus or higher taxa)")



# =======================================================
# CO-OCCURRENCE NETWORK (ROBUST TAXONOMIC LABELING)
# =======================================================

library(phyloseq)
library(igraph)
library(ggraph)
library(tidyverse)
library(Hmisc)
library(ggrepel)

# 1️⃣ Transform counts to relative abundance
ps_rel <- transform_sample_counts(physeq_combined, function(x) x / sum(x))

# 2️⃣ Aggregate to genus level
ps_genus <- tax_glom(ps_rel, taxrank = "Genus")

# 3️⃣ Extract taxonomy and abundance
tax_df <- as.data.frame(tax_table(ps_genus))
otu_df <- as.data.frame(otu_table(ps_genus))

# Ensure taxa are rows
if (!taxa_are_rows(ps_genus)) {
  otu_df <- t(otu_df)
}

# Make names unique to avoid duplication errors
tax_df$Taxon_unique <- make.unique(tax_df$Taxon_label)

# Apply to abundance matrix
rownames(otu_df) <- tax_df$Taxon_unique

# 5️⃣ Collapse duplicate taxa (sum across samples)
otu_df <- otu_df %>%
  rownames_to_column("Taxon") %>%
  group_by(Taxon) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Taxon")

# 6️⃣ Filter by mean relative abundance (≥0.2%)
mean_abund <- rowMeans(otu_df)
otu_df_filt <- otu_df[mean_abund >= 0.002, , drop = FALSE]

# 7️⃣ Compute Spearman correlations
abund_mat <- t(otu_df_filt)
storage.mode(abund_mat) <- "double"
cor_res <- Hmisc::rcorr(as.matrix(abund_mat), type = "spearman")

cor_mat <- cor_res$r
p_mat <- cor_res$P

# 8️⃣ Keep strong, significant correlations
r_cutoff <- 0.6
sig_mask <- (abs(cor_mat) >= r_cutoff) & (p_mat < 0.05)
cor_mat[!sig_mask] <- 0

# 9️⃣ Build igraph network
net <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
V(net)$Taxon <- V(net)$name

# 🔟 Determine shared vs host-specific taxa
sample_df <- as.data.frame(sample_data(physeq_combined))
otu_bin <- as.matrix(otu_df_filt) > 0

taxa_pres <- apply(otu_bin, 1, function(x) tapply(x, sample_df$Host, any))
if (is.list(taxa_pres)) taxa_pres <- do.call(cbind, taxa_pres)
shared_taxa <- names(which(rowSums(taxa_pres, na.rm = TRUE) > 1))

V(net)$Shared <- ifelse(V(net)$Taxon %in% shared_taxa, "Shared", "Host-specific")

# 1️⃣1️⃣ Plot the network
ggraph(net, layout = "fr") +
  geom_edge_link(aes(alpha = abs(weight)), color = "grey80") +
  geom_node_point(aes(color = Shared, size = degree(net))) +
  geom_node_text(aes(label = Taxon),
                 repel = TRUE,
                 size = 3,
                 check_overlap = TRUE) +
  scale_color_manual(values = c("Shared" = "orange", "Host-specific" = "steelblue")) +
  scale_size(range = c(2, 8)) +
  theme_void(base_size = 14) +
  ggtitle("Co-occurrence Network (≥0.2% abundance, hierarchical taxon labels)")



# Load libraries
library(phyloseq)
library(ggvenn)   # lightweight and beautiful Venns
library(tidyverse)

# ---- 1. Set parameters ----
taxrank_to_use <- "Species"   # You can also use "Family" or "Species"

# ---- 2. Function to extract taxa names from a phyloseq object ----
get_taxa_at_rank <- function(physeq_obj, rank = "Genus") {
  tax_df <- as.data.frame(tax_table(physeq_obj))
  
  # Handle missing taxonomy
  if (!rank %in% colnames(tax_df)) stop(paste("Rank", rank, "not found in tax_table"))
  taxa_names_rank <- as.character(tax_df[[rank]])
  taxa_names_rank <- taxa_names_rank[!is.na(taxa_names_rank) & taxa_names_rank != "" & taxa_names_rank != "?"]
  
  # Optionally remove unclassified taxa
  taxa_names_rank <- gsub("uncultured|unclassified|Unknown|unknown", "", taxa_names_rank)
  taxa_names_rank <- trimws(taxa_names_rank)
  taxa_names_rank[taxa_names_rank != ""]
}

# ---- 3. Extract taxa sets for each shrimp species ----
taxa_vannamei    <- get_taxa_at_rank(physeq1.vannamei, rank = taxrank_to_use)
taxa_rosenbergii <- get_taxa_at_rank(physeq1.rosenbergii, rank = taxrank_to_use)
taxa_monodon     <- get_taxa_at_rank(physeq1.monodon, rank = taxrank_to_use)

# ---- 4. Create a named list for Venn plotting ----
shrimp_sets <- list(
  "L. vannamei"    = unique(taxa_vannamei),
  "M. rosenbergii" = unique(taxa_rosenbergii),
  "P. monodon"     = unique(taxa_monodon)
)

# ---- 5. Plot the Venn diagram ----
ggvenn(
  shrimp_sets,
  fill_color = c("red", "blue", "green"),
  stroke_size = 0.8,
  set_name_size = 5,
  text_size = 5,
  show_percentage = FALSE
) +
  ggtitle(paste("Shared", taxrank_to_use, "across shrimp species")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))


# ============================================================
# VENN DIAGRAM + OVERLAP VALUES ON PLOT (SHRIMP MICROBIOMES)
# ============================================================

library(phyloseq)
library(ggVennDiagram)
library(tidyverse)
library(grid)

# -------------------------------
# 1️⃣ Parameters
# -------------------------------
taxrank_to_use <- "Genus"  # can also be "Family", "Order", etc.

# -------------------------------
# 2️⃣ Helper function to extract valid taxa
# -------------------------------
get_taxa_at_rank <- function(physeq_obj, rank = "Genus") {
  tax_df <- as.data.frame(tax_table(physeq_obj))
  taxa_names_rank <- as.character(tax_df[[rank]])
  taxa_names_rank <- taxa_names_rank[
    !is.na(taxa_names_rank) &
      taxa_names_rank != "" &
      taxa_names_rank != "?" &
      !grepl("uncultured|unclassified|Unknown|unknown", taxa_names_rank)
  ]
  unique(trimws(taxa_names_rank))
}

# -------------------------------
# 3️⃣ Extract taxa per shrimp species
# -------------------------------
taxa_vannamei    <- get_taxa_at_rank(physeq1.vannamei,    taxrank_to_use)
taxa_rosenbergii <- get_taxa_at_rank(physeq1.rosenbergii, taxrank_to_use)
taxa_monodon     <- get_taxa_at_rank(physeq1.monodon,     taxrank_to_use)

shrimp_sets <- list(
  "L. vannamei"    = taxa_vannamei,
  "M. rosenbergii" = taxa_rosenbergii,
  "P. monodon"     = taxa_monodon
)

# -------------------------------
# 4️⃣ Compute pairwise similarity (for annotation)
# -------------------------------
similarity_indices <- function(setA, setB) {
  inter <- length(intersect(setA, setB))
  union <- length(union(setA, setB))
  jaccard  <- inter / union
  sorensen <- (2 * inter) / (length(setA) + length(setB))
  list(Jaccard = round(jaccard, 2), Sorensen = round(sorensen, 2))
}

pairs <- list(
  c("L. vannamei", "M. rosenbergii"),
  c("L. vannamei", "P. monodon"),
  c("M. rosenbergii", "P. monodon")
)
pairwise_vals <- map(pairs, ~ {
  list(
    Pair = paste(.x[1], "vs", .x[2]),
    similarity_indices(shrimp_sets[[.x[1]]], shrimp_sets[[.x[2]]])
  )
})
pairwise_vals

# -------------------------------
# 5️⃣ Create the base Venn plot
# -------------------------------
venn_plot <- ggVennDiagram(
  shrimp_sets,
  label_alpha = 0,  # hide default small labels
  label = "count",
  edge_size = 0.8,
  set_size = 5
) +
  scale_fill_gradient(low = "#e0ecf4", high = "#8856a7") +
  labs(title = paste("Shared", taxrank_to_use, "across shrimp species")) +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -------------------------------
# 6️⃣ Add Jaccard values manually as annotations (optional)
# -------------------------------
# You can tweak coordinates depending on layout
venn_plot <- venn_plot +
  annotate("text", x = 0.0, y = 1.2, label = "Jaccard Similarity:", fontface = "bold", size = 4.5) +
  annotate("text", x = 0.0, y = 1.05,
           label = paste("vannamei vs rosenbergii =",
                         similarity_indices(taxa_vannamei, taxa_rosenbergii)$Jaccard),
           size = 4) +
  annotate("text", x = 0.0, y = 0.9,
           label = paste("vannamei vs monodon =",
                         similarity_indices(taxa_vannamei, taxa_monodon)$Jaccard),
           size = 4) +
  annotate("text", x = 0.0, y = 0.75,
           label = paste("rosenbergii vs monodon =",
                         similarity_indices(taxa_rosenbergii, taxa_monodon)$Jaccard),
           size = 4)

# -------------------------------
# 7️⃣ Print the plot
# -------------------------------
print(venn_plot)



################################################################################
##########################CORE MICORBIOME ######################################
################################################################################









################################################################################
################################################################################

sample_data(physeq_combined)

colnames(tax_table(physeq_combined))
colnames(tax_table(physeq_combined)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")



lefse_res <- run_lefse(
  physeq_combined,
  group = "Organism",       # your grouping variable
  norm = "CPM",         # normalization (Counts Per Million)
  kw_cutoff = 0.05,     # Kruskal–Wallis test cutoff
  lda_cutoff = 2.0,     # minimum LDA score
  multigrp_strat = TRUE # for multiclass data
)

lefse_res

unique(marker_table(lefse_res)$enrich_group)

library(RColorBrewer)
groups <- unique(marker_table(lefse_res)$enrich_group)
colour_map <- setNames(brewer.pal(length(groups), "Serosenbergii")[1:length(groups)], groups)
colour_map


plot_cladogram(lefse_res, color = colour_map)


