## =========================
## 0. SETUP
## =========================

rm(list = ls())

library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(RColorBrewer)

setwd("~/Library/CloudStorage/OneDrive-UniversityofKent/all_data")

required_cols <- c(
  "Organism",
  "Host",
  "Study",
  "Project",
  "Region_amplified"
)

prjna_codes <- c(
  "422950","1087723","495902","291010","417739","1012318","780955",
  "731310","749331","717320","797514","577421","600476","682282",
  "553862","540737","831321","505066","609648","573062","354668",
  "549032","579035"
)

## =========================
## 1. LOAD, ALIGN, RAREFY
## =========================

physeq_list <- list()
physeq_rarefied_list <- list()

for (prjna in prjna_codes) {
  
  cat("\n===============================\n")
  cat("Processing PRJNA", prjna, "\n")
  cat("===============================\n")
  
  otu_file  <- paste0("OTU_", prjna, ".csv")
  tax_file  <- paste0("hiera_BLAST_", prjna, ".csv")
  meta_file <- paste0("SraRunTable_", prjna, ".csv")
  
  otu_df  <- read.csv(otu_file,  row.names = 1, check.names = FALSE)
  tax_df  <- read.csv(tax_file,  row.names = 1, check.names = FALSE)
  meta_df <- read.csv(meta_file, row.names = 1, check.names = FALSE)
  
  ## --- check required metadata columns ---
  missing <- setdiff(required_cols, colnames(meta_df))
  if (length(missing) > 0) {
    stop(
      "❌ ERROR! The following columns are missing in ", meta_file, ":\n",
      paste(missing, collapse = ", "),
      "\nPlease edit the file and re-run.\n"
    )
  }
  
  ## --- make sure OTU and taxonomy rownames align ---
  common_otus <- intersect(rownames(otu_df), rownames(tax_df))
  dropped_otu <- setdiff(rownames(otu_df), common_otus)
  dropped_tax <- setdiff(rownames(tax_df), common_otus)
  
  if (length(dropped_otu) > 0) {
    cat("⚠ PRJNA", prjna, ": Dropping", length(dropped_otu),
        "OTUs not found in taxonomy.\n")
    otu_df <- otu_df[common_otus, , drop = FALSE]
  }
  if (length(dropped_tax) > 0) {
    cat("⚠ PRJNA", prjna, ": Dropping", length(dropped_tax),
        "taxonomy rows not found in OTU table.\n")
    tax_df <- tax_df[common_otus, , drop = FALSE]
  }
  
  ## --- build phyloseq object ---
  otu_tab  <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)
  tax_tab  <- tax_table(as.matrix(tax_df))
  meta_tab <- sample_data(meta_df)
  
  ## ensure OTU samples match metadata
  colnames(otu_tab) <- rownames(meta_df)
  
  physeq <- phyloseq(otu_tab, tax_tab, meta_tab)
  physeq_list[[prjna]] <- physeq
  
  ## --- rarefy ---
  physeq_rarefied <- rarefy_even_depth(
    physeq,
    sample.size = 10000,
    rngseed = 123,
    verbose = FALSE
  )
  
  physeq_rarefied_list[[prjna]] <- physeq_rarefied
  
  cat("Taxa after rarefaction:", ntaxa(physeq_rarefied), "\n")
}

## drop any empty sets (just in case)
physeq_rarefied_filtered <- Filter(
  function(x) ntaxa(x) > 0 && nsamples(x) > 0,
  physeq_rarefied_list
)

## =========================
## 2. FILTER TO BACTERIA, AGGREGATE, COMPOSITIONAL
## =========================

process_tax_level_safe <- function(phy, rank = "Phylum", mean_threshold = 0.002) {
  
  if (ntaxa(phy) == 0 || nsamples(phy) == 0) return(NULL)
  
  tax_df <- as.data.frame(as.matrix(tax_table(phy)))
  
  if (!"Domain" %in% colnames(tax_df)) {
    cat("❌ No Domain column — skipping dataset.\n")
    return(NULL)
  }
  
  ## --- keep only Bacteria ---
  keep_taxa <- rownames(tax_df)[tax_df$Domain == "Bacteria"]
  keep_taxa <- keep_taxa[!is.na(keep_taxa)]
  
  if (length(keep_taxa) < 1) {
    cat("❌ No bacterial taxa left — skipping dataset.\n")
    return(NULL)
  }
  
  phy <- prune_taxa(keep_taxa, phy)
  
  if (ntaxa(phy) < 1) {
    cat("❌ Dataset empty after bacteria filter — skipping.\n")
    return(NULL)
  }
  
  ## --- relative abundance ---
  relab <- transform_sample_counts(phy, function(x) x / sum(x))
  otu_mat <- as.matrix(otu_table(relab))
  otu_mat[is.na(otu_mat)] <- 0
  
  if (nrow(otu_mat) < 1) {
    cat("❌ OTU table empty after RA transform — skipping dataset.\n")
    return(NULL)
  }
  
  ## --- filter by mean abundance ---
  keep <- apply(otu_mat, 1, function(x) mean(x, na.rm = TRUE) > mean_threshold)
  keep[is.na(keep)] <- FALSE
  
  n_keep <- sum(keep)
  cat("   Taxa kept after abundance filter:", n_keep, "/", nrow(otu_mat), "\n")
  
  if (n_keep < 1) {
    cat("❌ All taxa removed by abundance filter — skipping dataset.\n")
    return(NULL)
  }
  
  phy_filt <- prune_taxa(names(keep)[keep], phy)
  
  if (ntaxa(phy_filt) < 1) {
    cat("❌ Dataset empty after prune_taxa — skipping.\n")
    return(NULL)
  }
  
  ## --- aggregate at chosen rank ---
  agg <- tryCatch(
    microbiome::aggregate_taxa(phy_filt, level = rank),
    error = function(e) {
      cat("❌ aggregate_taxa failed:", e$message, "\n")
      NULL
    }
  )
  if (is.null(agg) || ntaxa(agg) < 1) return(NULL)
  
  ## --- compositional (again, but now at rank level) ---
  comp <- tryCatch(
    microbiome::transform(agg, "compositional"),
    error = function(e) {
      cat("❌ transform('compositional') failed:", e$message, "\n")
      NULL
    }
  )
  if (is.null(comp) || ntaxa(comp) < 1) return(NULL)
  
  return(comp)
}

## =========================
## 3. BUILD LONG DATA FOR EACH RANK (NO MERGE OF PHYLOSEQ)
## =========================

build_composition_df <- function(phy_list, rank = "Phylum", mean_threshold = 0.002) {
  
  all_dfs <- list()
  
  for (nm in names(phy_list)) {
    cat("\n=== Processing", nm, "at rank", rank, "===\n")
    
    ps <- process_tax_level_safe(phy_list[[nm]], rank = rank, mean_threshold = mean_threshold)
    if (is.null(ps)) next
    
    ## psmelt: long format with Sample, Abundance, taxonomy, metadata
    df <- phyloseq::psmelt(ps)
    
    if (!all(required_cols %in% colnames(df))) {
      cat("⚠ Skipping", nm, "— missing one of required metadata columns in melted data.\n")
      next
    }
    
    ## rank column (e.g. "Phylum", "Class", ...)
    if (!rank %in% colnames(df)) {
      cat("⚠ Skipping", nm, "— rank", rank, "not found in psmelt().\n")
      next
    }
    
    df_rank <- df %>%
      dplyr::select(
        Sample,
        Organism,
        Host,
        Study,
        Project,
        Region_amplified,
        Taxon = dplyr::all_of(rank),
        Abundance
      )
    
    ## drop NA / empty taxa names
    df_rank <- df_rank %>%
      filter(!is.na(Taxon), Taxon != "")
    
    all_dfs[[nm]] <- df_rank
  }
  
  if (length(all_dfs) == 0) {
    stop("No datasets survived processing at rank ", rank)
  }
  
  bind_rows(all_dfs, .id = "PRJNA")
}

## =========================
## 4. PALETTE HELPER
## =========================

make_tax_palette <- function(taxa_names) {
  uniq <- unique(taxa_names)
  n <- length(uniq)
  if (n <= 12) {
    cols <- brewer.pal(max(n, 3), "Set3")[seq_len(n)]
  } else {
    cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
  }
  names(cols) <- uniq
  cols
}

## =========================
## 5. GENERATE PLOTS FOR EACH RANK
## =========================

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
plots_list <- list()
comp_data_list <- list()   # if you want the tables too

for (rank in tax_levels) {
  
  cat("\n===============================\n")
  cat("Processing rank:", rank, "\n")
  cat("===============================\n")
  
  comp_df <- build_composition_df(physeq_rarefied_filtered, rank = rank, mean_threshold = 0.002)
  comp_data_list[[rank]] <- comp_df
  
  ## average by Organism (normalised so bars sum to 1)
  comp_org <- comp_df %>%
    group_by(Organism, Taxon) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Organism) %>%
    mutate(Abundance = mean_abund / sum(mean_abund)) %>%
    ungroup()
  
  
  ## palette
  pal <- make_tax_palette(comp_org$Taxon)
  
  p <- ggplot(comp_org, aes(x = Organism, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal, drop = FALSE) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    ggtitle(paste("Composition at", rank)) +
    ylab("Mean relative abundance") +
    xlab("Organism")
  
  plots_list[[rank]] <- p
}

## =========================
## 6. SHOW PLOTS
## =========================

plots_list[["Phylum"]]
plots_list[["Class"]]
plots_list[["Order"]]
plots_list[["Family"]]
plots_list[["Genus"]]



####################################################################
######################## PCA #######################################
####################################################################

# physeq_rarefied_filtered = your list of final cleaned phyloseq objects

# Label each phyloseq object with its PRJNA name
phy_list_labeled <- imap(physeq_rarefied_filtered, ~{
  ps <- .x
  proj <- .y
  
  # Ensure sample_data exists
  if (is.null(sample_data(ps, errorIfNULL = FALSE))) {
    sample_data(ps) <- data.frame(row.names = sample_names(ps))
  }
  
  sample_data(ps)$Project <- proj
  ps
})

# Merge EVERYTHING into one phyloseq object
phy_all <- Reduce(merge_phyloseq, phy_list_labeled)


# Extract OTU table as matrix
otu <- as(otu_table(phy_all), "matrix")
if (taxa_are_rows(phy_all)) otu <- t(otu)

# Run PCA
pca <- prcomp(otu, scale. = TRUE, center = TRUE)

# Convert PCA scores to dataframe
pca_df <- as.data.frame(pca$x) %>%
  tibble::rownames_to_column("SampleID")

# Metadata
meta <- sample_data(phy_all) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("SampleID")

# Merge PCA + metadata
pca_df <- left_join(pca_df, meta, by = "SampleID")

# Ensure Project is a clean factor
pca_df$Project <- as.character(pca_df$Project)
pca_df$Project <- factor(pca_df$Project, levels = sort(unique(pca_df$Project)))

library(Polychrome)

# number of PRJNA groups
n_proj <- length(levels(pca_df$Project))
n_proj <- as.integer(n_proj)   # MUST be integer for Polychrome

# generate palette safely
pal_vec <- tryCatch(
  {
    if (n_proj <= 22) {
      Polychrome::kelly(n_proj)
    } else {
      Polychrome::alphabet(n_proj)
    }
  },
  error = function(e) {
    message("⚠ Polychrome failed, switching to fallback palette.")
    grDevices::rainbow(n_proj)  # fallback for safety
  }
)

# assign names for ggplot
names(pal_vec) <- levels(pca_df$Project)



# Variance explained
pca_var <- pca$sdev^2
pca_var_exp <- pca_var / sum(pca_var)
pc1_lab <- paste0("PC1 (", round(100 * pca_var_exp[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * pca_var_exp[2], 1), "%)")


ggplot(pca_df, aes(PC1, PC2, colour = Project)) +
  geom_point(size = 2) +
  facet_wrap(~ Organism) +
  scale_colour_manual(values = pal_vec) +
  theme_minimal() +
  labs(
    title = "PCA of Microbial Communities by Organism",
    color = "PRJNA Project",
    x = pc1_lab,
    y = pc2_lab
  )


####################################################################
######################## Alpha Diversity ###########################
####################################################################


library(phyloseq)
library(dplyr)
library(tibble)
library(ggpubr)
library(rstatix)
library(patchwork)

## ---------- 1. Estimate alpha diversity ----------

# Use your new merged phyloseq object:
#   phy_all = merged, cleaned phyloseq with Organism in sample_data
alpha_diversity <- estimate_richness(phy_all)
print(alpha_diversity)

# Add sample ID column
alpha_diversity$Sample <- rownames(alpha_diversity)

# Select indices of interest
selected_df <- alpha_diversity[, c("Sample", "Observed", "Chao1", "Simpson", "Shannon")]
rownames(selected_df) <- NULL

## ---------- 2. Build metadata + statframe ----------

# metadata from phy_all
metadata <- sample_data(phy_all) %>%
  data.frame(stringsAsFactors = FALSE)
metadata$Sample <- rownames(metadata)

# variables for statframe
samples      <- metadata$Sample
organisms    <- metadata$Organism
shannonscore <- alpha_diversity$Shannon
simpsonscore <- alpha_diversity$Simpson
chaoscore    <- alpha_diversity$Chao1
observed     <- alpha_diversity$Observed

statframe <- data.frame(samples, organisms, shannonscore, simpsonscore, chaoscore, observed)
head(statframe[1:5, 1:5])

## ---------- 3. Normality checks ----------

shapiro.shan     <- shapiro_test(alpha_diversity$Shannon)$p.value
shapiro.chao     <- shapiro_test(alpha_diversity$Chao1)$p.value
shapiro.simp     <- shapiro_test(alpha_diversity$Simpson)$p.value
shapiro.observed <- shapiro_test(alpha_diversity$Observed)$p.value

# quick visual check
histogram(alpha_diversity$Shannon)

## ---------- 4. Global tests (ANOVA or Kruskal) ----------

if (shapiro.shan < 0.05) {
  result.kruskal.shannon <- statframe %>% kruskal_test(shannonscore ~ organisms)
}
if (shapiro.shan > 0.05) {
  result.anova.shannon <- statframe %>% anova_test(shannonscore ~ organisms)
}

if (shapiro.simp < 0.05) {
  result.kruskal.simpson <- statframe %>% kruskal_test(simpsonscore ~ organisms)
}
if (shapiro.simp > 0.05) {
  result.anova.simpson <- statframe %>% anova_test(simpsonscore ~ organisms)
}

if (shapiro.chao < 0.05) {
  result.kruskal.chao <- statframe %>% kruskal_test(chaoscore ~ organisms)
}
if (shapiro.chao > 0.05) {
  result.anova.chao <- statframe %>% anova_test(chaoscore ~ organisms)   # fixed typo
}

if (shapiro.observed < 0.05) {
  result.kruskal.observed <- statframe %>% kruskal_test(observed ~ organisms)
}
if (shapiro.observed > 0.05) {
  result.anova.observed <- statframe %>% anova_test(observed ~ organisms)
}

## ---------- 5. Pairwise tests ----------

if (exists("result.kruskal.shannon")) {
  pwc.shannon <- statframe %>%
    dunn_test(shannonscore ~ organisms, p.adjust.method = "bonferroni")
}
if (exists("result.anova.shannon")) {
  pwc.shannon <- statframe %>%
    tukey_hsd(shannonscore ~ organisms)
}

if (exists("result.kruskal.simpson")) {
  pwc.simpson <- statframe %>%
    dunn_test(simpsonscore ~ organisms, p.adjust.method = "bonferroni")
}
if (exists("result.anova.simpson")) {
  pwc.simpson <- statframe %>%
    tukey_hsd(simpsonscore ~ organisms)
}

if (exists("result.kruskal.chao")) {
  pwc.chao <- statframe %>%
    dunn_test(chaoscore ~ organisms, p.adjust.method = "bonferroni")
}
if (exists("result.anova.chao")) {
  pwc.chao <- statframe %>%
    tukey_hsd(chaoscore ~ organisms)
}

if (exists("result.kruskal.observed")) {
  pwc.observed <- statframe %>%
    dunn_test(observed ~ organisms, p.adjust.method = "bonferroni")
}
if (exists("result.anova.observed")) {
  pwc.observed <- statframe %>%
    tukey_hsd(observed ~ organisms)
}

################################
## 6. PLOTS WITH P AIRWISE P  ##
################################

################# SHANNON #################

pwc.shannon <- pwc.shannon %>% add_xy_position(x = "organisms")

statplot.shannon <- ggviolin(
  statframe,
  x = "organisms",
  y = "shannonscore",
  title = "Shannon Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  stat_pvalue_manual(
    pwc.shannon,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(caption = get_pwc_label(pwc.shannon))

if (exists("result.kruskal.shannon")) {
  statplot.shannon <- statplot.shannon +
    labs(subtitle = get_test_label(result.kruskal.shannon, detailed = FALSE))
}
if (exists("result.anova.shannon")) {
  statplot.shannon <- statplot.shannon +
    labs(subtitle = get_test_label(result.anova.shannon, detailed = FALSE))
}

statplot.shannon

################# SIMPSON #################

pwc.simpson <- pwc.simpson %>% add_xy_position(x = "organisms")

statplot.simpson <- ggviolin(
  statframe,
  x = "organisms",
  y = "simpsonscore",
  title = "Simpson Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  stat_pvalue_manual(
    pwc.simpson,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(caption = get_pwc_label(pwc.simpson))

if (exists("result.kruskal.simpson")) {
  statplot.simpson <- statplot.simpson +
    labs(subtitle = get_test_label(result.kruskal.simpson, detailed = FALSE))
}
if (exists("result.anova.simpson")) {
  statplot.simpson <- statplot.simpson +
    labs(subtitle = get_test_label(result.anova.simpson, detailed = FALSE))
}

################# CHAO #################

pwc.chao <- pwc.chao %>% add_xy_position(x = "organisms")

statplot.chao <- ggviolin(
  statframe,
  x = "organisms",
  y = "chaoscore",
  title = "Chao Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  stat_pvalue_manual(
    pwc.chao,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(caption = get_pwc_label(pwc.chao))

if (exists("result.kruskal.chao")) {
  statplot.chao <- statplot.chao +
    labs(subtitle = get_test_label(result.kruskal.chao, detailed = FALSE))
}
if (exists("result.anova.chao")) {
  statplot.chao <- statplot.chao +
    labs(subtitle = get_test_label(result.anova.chao, detailed = FALSE))
}

################# OBSERVED #################

pwc.observed <- pwc.observed %>% add_xy_position(x = "organisms")

statplot.observed <- ggviolin(
  statframe,
  x = "organisms",
  y = "observed",
  title = "Observed Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  stat_pvalue_manual(
    pwc.observed,
    hide.ns = TRUE,
    step.increase = 0.1
  ) +
  labs(caption = get_pwc_label(pwc.observed))

if (exists("result.kruskal.observed")) {
  statplot.observed <- statplot.observed +
    labs(subtitle = get_test_label(result.kruskal.observed, detailed = FALSE))
}
if (exists("result.anova.observed")) {
  statplot.observed <- statplot.observed +
    labs(subtitle = get_test_label(result.anova.observed, detailed = FALSE))
}

## --- Combine all four in 2x2 grid ---
(statplot.shannon | statplot.simpson) /
  (statplot.chao | statplot.observed) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")














statplot.chao <- ggviolin(
  statframe,
  x = "organisms",
  y = "chaoscore",
  title = "Chao Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +   # <-- JITTER ADDED
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Add test subtitle (ANOVA or KW)
if (exists("result.kruskal.chao")) {
  statplot.chao <- statplot.chao +
    labs(subtitle = get_test_label(result.kruskal.chao, detailed = FALSE))
}
if (exists("result.anova.chao")) {
  statplot.chao <- statplot.chao +
    labs(subtitle = get_test_label(result.anova.chao, detailed = FALSE))
}
statplot.chao

statplot.shannon <- ggviolin(
  statframe,
  x = "organisms",
  y = "shannonscore",
  title = "Shannon Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

if (exists("result.kruskal.shannon")) {
  statplot.shannon <- statplot.shannon +
    labs(subtitle = get_test_label(result.kruskal.shannon, detailed = FALSE))
}
if (exists("result.anova.shannon")) {
  statplot.shannon <- statplot.shannon +
    labs(subtitle = get_test_label(result.anova.shannon, detailed = FALSE))
}
statplot.shannon



statplot.simpson <- ggviolin(
  statframe,
  x = "organisms",
  y = "simpsonscore",
  title = "Simpson Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

if (exists("result.kruskal.simpson")) {
  statplot.simpson <- statplot.simpson +
    labs(subtitle = get_test_label(result.kruskal.simpson, detailed = FALSE))
}
if (exists("result.anova.simpson")) {
  statplot.simpson <- statplot.simpson +
    labs(subtitle = get_test_label(result.anova.simpson, detailed = FALSE))
}

statplot.simpson


statplot.observed <- ggviolin(
  statframe,
  x = "organisms",
  y = "observed",
  title = "Observed Diversity of Species",
  fill = "organisms",
  palette = "Paired"
) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

if (exists("result.kruskal.observed")) {
  statplot.observed <- statplot.observed +
    labs(subtitle = get_test_label(result.kruskal.observed, detailed = FALSE))
}
if (exists("result.anova.observed")) {
  statplot.observed <- statplot.observed +
    labs(subtitle = get_test_label(result.anova.observed, detailed = FALSE))
}
statplot.observed



(statplot.shannon | statplot.simpson) /
  (statplot.chao    | statplot.observed) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


####################################################################
######################## PCoA ######################################
####################################################################



library(phyloseq)
library(vegan)
library(ggplot2)

## use your merged phyloseq object
ps <- phy_all   # or physeq_combined if you kept that name

## Bray–Curtis distance
dist_bray <- phyloseq::distance(ps, method = "bray")

## PCoA ordination
ord_bray <- ordinate(ps, method = "PCoA", distance = dist_bray)

## metadata
meta <- data.frame(sample_data(ps))

## ----- PERMANOVA by Host (if column exists) -----
if ("Host" %in% colnames(meta)) {
  # make sure rows of metadata line up with distance labels
  samples <- attr(dist_bray, "Labels")
  meta_host <- meta[samples, , drop = FALSE]
  
  adonis_host <- adonis2(dist_bray ~ Host, data = meta_host, permutations = 999)
  pval_host   <- adonis_host$`Pr(>F)`[1]
  r2_host     <- adonis_host$R2[1]
  label_host  <- paste0("PERMANOVA: R² = ",
                        round(r2_host, 3),
                        ", p = ",
                        signif(pval_host, 3))
  
  p_bray_host <- plot_ordination(ps, ord_bray, color = "Host") +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title    = "Bray–Curtis PCoA by Host",
      subtitle = label_host
    )
}

## ----- PERMANOVA by Organism (if column exists) -----
if ("Organism" %in% colnames(meta)) {
  samples <- attr(dist_bray, "Labels")
  meta_org <- meta[samples, , drop = FALSE]
  
  adonis_org <- adonis2(dist_bray ~ Organism, data = meta_org, permutations = 999)
  pval_org   <- adonis_org$`Pr(>F)`[1]
  r2_org     <- adonis_org$R2[1]
  label_org  <- paste0("PERMANOVA: R² = ",
                       round(r2_org, 3),
                       ", p = ",
                       signif(pval_org, 3))
  
  p_bray_org <- plot_ordination(ps, ord_bray, color = "Organism") +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title    = "Bray–Curtis PCoA by Organism",
      subtitle = label_org
    )
}

## ---- show them ----
if (exists("p_bray_host")) print(p_bray_host)
if (exists("p_bray_org"))  print(p_bray_org)

library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)

ps <- phy_all       # your merged phyloseq object
group_var <- "Organism"

methods <- c("bray", "jaccard_pa", "canberra", "ruzicka")

make_pcoa_plot <- function(ps, method, group_var = "Organism") {
  
  # ---- Distance matrix ----
  if (method == "jaccard_pa") {
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu <- t(otu)
    dist_mat <- vegan::vegdist(otu, method = "jaccard", binary = TRUE)
    method_title <- "Jaccard (presence/absence)"
    
  } else if (method == "ruzicka") {
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu <- t(otu)
    dist_mat <- vegan::vegdist(otu, method = "jaccard", binary = FALSE)
    method_title <- "Ružička (quantitative Jaccard)"
    
  } else {
    dist_mat <- phyloseq::distance(ps, method = method)
    method_title <- toupper(method)
  }
  
  # ---- Align metadata ----
  samples <- attr(dist_mat, "Labels")
  meta <- data.frame(sample_data(ps))[samples, , drop = FALSE]
  meta$Group <- meta[[group_var]]
  
  # ---- Ordination ----
  ord <- ordinate(ps, method = "PCoA", distance = dist_mat)
  
  # ---- PERMANOVA ----
  ad <- vegan::adonis2(dist_mat ~ Group, data = meta, permutations = 999)
  r2 <- round(ad$R2[1], 3)
  p  <- signif(ad$`Pr(>F)`[1], 3)
  
  # ---- Plot ----
  p_plot <- plot_ordination(ps, ord, color = group_var) +
    geom_point(size = 2.5, alpha = 0.8) +
    theme_minimal() +
    labs(
      title    = paste("PCoA –", method_title),
      subtitle = paste0("PERMANOVA: R² = ", r2, "  p = ", p),
      color    = group_var
    ) +
    theme(
      legend.position = "right",
      plot.title    = element_text(face = "bold")
    )
  
  list(plot = p_plot, dist = dist_mat)
}

# ---- Build plots ----
plots <- list()
dist_mats <- list()

for (m in methods) {
  res <- make_pcoa_plot(ps, m, group_var)
  plots[[m]] <- res$plot
  dist_mats[[m]] <- res$dist
}

# ---- Combine ----
combined_pcoa_plot <- wrap_plots(plots, ncol = 2, guides = "collect") &
  theme(legend.position = "right")

print(combined_pcoa_plot)


###############################################################################
######################### CORE MICROBIOME ANALYSIS ############################
###############################################################################

# Choose rank for core microbiome (you can change this)
core_rank <- "Class"   # or "Phylum", "Order", "Family", "Genus"

# Extract all processed, filtered, aggregated phyloseq objects at this rank
core_ps_list <- processed_phyloseq_by_rank[[core_rank]]

if (length(core_ps_list) == 0) {
  stop("❌ No phyloseq objects available for core microbiome at rank: ", core_rank)
}

cat("\n=== Building merged phyloseq for CORE microbiome at rank:", core_rank, "===\n")

# -------------------------------------------------------------------------
# 1. Union of all taxa
# -------------------------------------------------------------------------
all_taxa <- Reduce(union, lapply(core_ps_list, taxa_names))

# -------------------------------------------------------------------------
# 2. Build merged OTU table
# -------------------------------------------------------------------------
otu_fixed <- lapply(core_ps_list, function(ps) {
  mat <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) mat <- t(mat)
  
  out <- matrix(
    0,
    nrow = length(all_taxa),
    ncol = ncol(mat),
    dimnames = list(all_taxa, colnames(mat))
  )
  out[rownames(mat), ] <- mat
  
  otu_table(out, taxa_are_rows = TRUE)
})

otu_merged_core <- do.call(cbind, otu_fixed)
## =========================
## 3. BUILD LONG DATA FOR EACH RANK (NO MERGE OF PHYLOSEQ)
## =========================

build_composition_df <- function(phy_list, rank = "Phylum", mean_threshold = 0.002) {
  
  all_dfs <- list()
  
  for (nm in names(phy_list)) {
    cat("\n=== Processing", nm, "at rank", rank, "===\n")
    
    ps <- process_tax_level_safe(phy_list[[nm]], rank = rank, mean_threshold = mean_threshold)
    if (is.null(ps)) next
    
    # psmelt: long format with Sample, Abundance, taxonomy, metadata
    df <- phyloseq::psmelt(ps)
    
    # check metadata columns exist
    if (!all(required_cols %in% colnames(df))) {
      cat("⚠ Skipping", nm, "— missing one of required metadata columns in melted data.\n")
      next
    }
    
    # check rank column exists
    if (!rank %in% colnames(df)) {
      cat("⚠ Skipping", nm, "— rank", rank, "not found in psmelt().\n")
      next
    }
    
    df_rank <- df %>%
      dplyr::select(
        Sample,
        Organism,
        Host,
        Study,
        Project,
        Region_amplified,
        Taxon = dplyr::all_of(rank),
        Abundance
      ) %>%
      dplyr::filter(!is.na(Taxon), Taxon != "")
    
    all_dfs[[nm]] <- df_rank
  }
  
  if (length(all_dfs) == 0) {
    stop("No datasets survived processing at rank ", rank)
  }
  
  dplyr::bind_rows(all_dfs, .id = "PRJNA")
}

## =========================
## 4. PALETTE HELPER
## =========================

make_tax_palette <- function(taxa_names) {
  uniq <- unique(taxa_names)
  n <- length(uniq)
  if (n <= 12) {
    cols <- RColorBrewer::brewer.pal(max(n, 3), "Set3")[seq_len(n)]
  } else {
    cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
  names(cols) <- uniq
  cols
}

## =========================
## 5. GENERATE BAR PLOTS + STORE TABLES FOR EACH RANK
## =========================

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")
plots_list      <- list()
comp_data_list  <- list()   # we’ll use this for core microbiome later

for (rank in tax_levels) {
  
  cat("\n===============================\n")
  cat("Processing rank:", rank, "\n")
  cat("===============================\n")
  
  comp_df <- build_composition_df(
    physeq_rarefied_filtered,
    rank = rank,
    mean_threshold = 0.002
  )
  comp_data_list[[rank]] <- comp_df
  
  # average by Organism, normalised so each Organism bar sums to 1
  comp_org <- comp_df %>%
    dplyr::group_by(Organism, Taxon) %>%
    dplyr::summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Organism) %>%
    dplyr::mutate(Abundance = mean_abund / sum(mean_abund)) %>%
    dplyr::ungroup()
  
  pal <- make_tax_palette(comp_org$Taxon)
  
  p <- ggplot2::ggplot(comp_org, ggplot2::aes(x = Organism, y = Abundance, fill = Taxon)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x         = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x  = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle(paste("Composition at", rank)) +
    ggplot2::ylab("Mean relative abundance") +
    ggplot2::xlab("Organism")
  
  plots_list[[rank]] <- p
}

## =========================
## 6. BUILD COMBINED PHYLOSEQ FROM LONG DATA (FOR CORE)
## =========================

build_phyloseq_from_comp <- function(comp_df, rank = "Class") {
  
  # 1. Aggregate abundance per Sample–Taxon (just in case)
  comp_agg <- comp_df %>%
    dplyr::group_by(Sample, Taxon) %>%
    dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  # 2. OTU table: taxa x samples
  otu_wide <- comp_agg %>%
    tidyr::pivot_wider(
      names_from  = Sample,
      values_from = Abundance,
      values_fill = 0
    )
  
  taxa <- otu_wide$Taxon
  otu_mat <- as.matrix(otu_wide[, -1, drop = FALSE])
  rownames(otu_mat) <- taxa
  
  otu_tab <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  
  # 3. Sample metadata from comp_df (one row per Sample)
  meta <- comp_df %>%
    dplyr::select(Sample, dplyr::all_of(required_cols)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ {
          x <- .x
          x <- x[!is.na(x)]
          if (length(x) == 0) NA_character_ else as.character(x[1])
        }
      ),
      .groups = "drop"
    )
  
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df$Sample
  meta_df$Sample <- NULL
  
  sd_tab <- phyloseq::sample_data(meta_df)
  
  # 4. Minimal taxonomy table for the chosen rank
  tax_df <- data.frame(row.names = taxa)
  tax_df[[rank]] <- taxa   # just store the taxon name under this rank
  tax_tab <- phyloseq::tax_table(as.matrix(tax_df))
  
  phyloseq::phyloseq(otu_tab, tax_tab, sd_tab)
}

## =========================
## 7. CORE MICROBIOME AT A GIVEN RANK (e.g. CLASS)
## =========================

# Choose the rank you want for core analysis:
core_rank <- "Class"   # or "Phylum", "Order", "Family", "Genus"

# 1. Get the long compositional data for this rank
comp_df_core <- comp_data_list[[core_rank]]

# 2. Build a combined phyloseq object from it
phy_core <- build_phyloseq_from_comp(comp_df_core, rank = core_rank)

# 3. Ensure compositional (should already be, but safe)
phy_core_comp <- microbiome::transform(phy_core, "compositional")

# 4. Define detection & prevalence thresholds (as in your original script)
detections  <- round(10^seq(log10(0.01), log10(0.2), length = 9), 3)
prevalences <- seq(0.05, 1, 0.05)
core_cols   <- c("yellow", "blue")

# 5. Plot core microbiome heatmap
core1 <- microbiome::plot_core(
  phy_core_comp,
  plot.type    = "heatmap",
  prevalences  = prevalences,
  colours      = core_cols,
  detections   = detections,
  min.prevalence = 0.8
) +
  ggplot2::labs(x = "Detection Threshold\n(Relative Abundance (%))")

core1



## ============================================================
## 8. CORE MICROBIOME SEPARATELY FOR EACH ORGANISM 
## ============================================================

# detection & prevalence settings (same as earlier)
detections  <- round(10^seq(log10(0.01), log10(0.2), length = 9), 3)
prevalences <- seq(0.05, 1, 0.05)
core_cols   <- c("yellow", "blue")

# Extract list of species present
species_list <- unique(sample_data(phy_core_comp)$Organism)

core_plots_by_species <- list()

cat("\n===============================\n")
cat("Generating core microbiome plots per species (Organism)\n")
cat("===============================\n")

for (sp in species_list) {
  
  cat("\n--- Processing species:", sp, "---\n")
  
  # subset phyloseq object to this species
  ps_sp <- subset_samples(phy_core_comp, Organism == sp)
  
  # Drop taxa not present after subsetting
  ps_sp <- prune_taxa(taxa_sums(ps_sp) > 0, ps_sp)
  
  # Skip if fewer than 2 samples remain
  if (nsamples(ps_sp) < 2) {
    cat("⚠ Skipping", sp, ": too few samples for core analysis.\n")
    next
  }
  
  # Core microbiome heatmap for this species
  p <- microbiome::plot_core(
    ps_sp,
    plot.type      = "heatmap",
    prevalences    = prevalences,
    colours        = core_cols,
    detections     = detections,
    min.prevalence = 0.8
  ) +
    ggplot2::ggtitle(paste("Core Microbiome (", core_rank, ") —", sp)) +
    ggplot2::labs(x = "Detection Threshold\n(Relative Abundance (%))")
  
  core_plots_by_species[[sp]] <- p
}

cat("\n✓ DONE! Use core_plots_by_species$`SpeciesName` to view plots.\n")

core_plots_by_species$`L.vannamei gut metagenome`
core_plots_by_species$`C.quadricarinatus gut metagenome`
core_plots_by_species$`C.cainii gut metagenome`
core_plots_by_species$`P.argus gut metagenome`
core_plots_by_species$`P.monodon gut metagenome`
core_plots_by_species$`P.clarkii gut metagenome`
core_plots_by_species$`P.leniusculus gut metagenome`
core_plots_by_species$`M.nipponese gut metagenome`
x


###############################################################################
######################## HEATMAP PER SPECIES (TOP N TAXA) #####################
###############################################################################

library(pheatmap)
library(dplyr)
library(tidyr)

## 1. Choose taxonomic rank + top N taxa
heat_rank <- "Genus"   # or "Phylum", "Class", "Order", "Family"
top_n     <- 30

## 2. Build long composition table at that rank
##    This reuses your existing filtering (bacteria-only, abundance threshold)
comp_df_heat <- build_composition_df(
  physeq_rarefied_filtered,
  rank = heat_rank,
  mean_threshold = 0.002
)

###############################################################################
############### HEATMAP PER SPECIES — TAXA ≥ 1% RELATIVE ABUND ################
###############################################################################

library(pheatmap)
library(dplyr)
library(tidyr)

## Set rank used in your earlier filtering:
heat_rank <- "Genus"     # can be "Phylum", "Class", etc.

## Threshold for mean relative abundance:
abundance_cutoff <- 0.01   # 1%

## comp_df_heat must already exist from:
# comp_df_heat <- build_composition_df(...)

## 1. Function to make heatmap for ONE species
make_species_heatmap_over1 <- function(comp_df, species_name, tax_label = heat_rank,
                                       cutoff = 0.01) {
  
  message("\n=== Building ≥1% heatmap for species: ", species_name, " ===")
  
  # Subset to this species
  sub <- comp_df %>%
    filter(Organism == species_name)
  
  if (nrow(sub) == 0) {
    warning("No rows for species: ", species_name)
    return(NULL)
  }
  
  # Build taxon × sample matrix
  mat_df <- sub %>%
    select(Sample, Taxon, Abundance) %>%
    group_by(Taxon, Sample) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    pivot_wider(
      names_from = Sample,
      values_from = Abundance,
      values_fill = 0
    )
  
  # Extract matrix
  taxa <- mat_df$Taxon
  mat_df$Taxon <- NULL
  abund_mat <- as.matrix(mat_df)
  rownames(abund_mat) <- taxa
  
  # Compute mean abundance per taxon
  mean_abund <- rowMeans(abund_mat)
  
  # Filter ≥ 1% mean relative abundance
  keep_taxa <- names(mean_abund[mean_abund >= cutoff])
  
  if (length(keep_taxa) == 0) {
    warning("No taxa ≥ ", cutoff, " for species ", species_name)
    return(NULL)
  }
  
  abund_mat <- abund_mat[keep_taxa, , drop = FALSE]
  
  # Log transform for visualization
  abund_mat_plot <- log10(abund_mat + 1e-5)
  
  # Build sample metadata annotation
  annot <- sub %>%
    select(Sample, Organism, Host, Study, Project, Region_amplified) %>%
    distinct()
  
  # Align annotation to samples in matrix
  annot <- annot[match(colnames(abund_mat_plot), annot$Sample), , drop = FALSE]
  rownames(annot) <- annot$Sample
  annot$Sample <- NULL
  
  # Draw heatmap
  pheatmap(
    abund_mat_plot,
    annotation_col = annot,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    fontsize = 9,
    main = paste0("Taxa ≥ 1% (", tax_label, ") — ", species_name)
  )
}

## 2. Generate heatmaps for all species
species_vec <- sort(unique(comp_df_heat$Organism))

heatmaps_over1_by_species <- list()

for (sp in species_vec) {
  heatmaps_over1_by_species[[sp]] <- make_species_heatmap_over1(
    comp_df = comp_df_heat,
    species_name = sp,
    tax_label = heat_rank,
    cutoff = abundance_cutoff
  )
}

## 3. Example: show heatmap for L. vannamei
# heatmaps_over1_by_species$`L.vannamei gut metagenome`



## 3. Function to make ONE heatmap for a given species
make_species_heatmap_from_df <- function(comp_df, species_name,
                                         top_n = 30,
                                         tax_label = heat_rank) {
  
  message("\n=== Building heatmap for species: ", species_name, " ===")
  
  # Subset to this species
  sub <- comp_df %>%
    dplyr::filter(Organism == species_name)
  
  if (nrow(sub) == 0) {
    warning("No rows for species: ", species_name, " — skipping.")
    return(NULL)
  }
  
  # 3a. Build Taxon x Sample matrix of abundances
  mat_df <- sub %>%
    dplyr::select(Sample, Taxon, Abundance) %>%
    dplyr::group_by(Taxon, Sample) %>%
    dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = Sample,
      values_from = Abundance,
      values_fill = 0
    )
  
  # rownames = taxa
  tax_labels <- mat_df$Taxon
  mat_df$Taxon <- NULL
  abund_mat <- as.matrix(mat_df)
  rownames(abund_mat) <- tax_labels
  
  # 3b. Pick top N taxa by total abundance
  if (nrow(abund_mat) == 0) {
    warning("No taxa left for species: ", species_name)
    return(NULL)
  }
  
  n_taxa_available <- nrow(abund_mat)
  n_to_take <- min(top_n, n_taxa_available)
  top_taxa <- names(sort(rowSums(abund_mat), decreasing = TRUE))[1:n_to_take]
  abund_mat <- abund_mat[top_taxa, , drop = FALSE]
  
  # 3c. Optional log transform
  abund_mat_plot <- log10(abund_mat + 1e-5)
  
  # 3d. Build sample annotation from comp_df
  annot <- sub %>%
    dplyr::select(
      Sample,
      Organism,
      Host,
      Study,
      Project,
      Region_amplified
    ) %>%
    dplyr::distinct()
  
  # Align annotation rows to matrix columns
  annot <- annot %>%
    dplyr::filter(Sample %in% colnames(abund_mat_plot))
  
  rownames(annot) <- annot$Sample
  annot$Sample <- NULL
  
  # Final safety check
  annot <- annot[colnames(abund_mat_plot), , drop = FALSE]
  
  # 3e. Draw heatmap
  pheatmap(
    abund_mat_plot,
    annotation_col = annot,
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    scale          = "row",
    fontsize       = 9,
    main           = paste0(
      "Top ", n_to_take, " ", tax_label,
      " — ", species_name
    )
  )
}

## 4. Build heatmaps for ALL species (stored in a list, but also plotted)
species_vec <- sort(unique(comp_df_heat$Organism))

heatmaps_by_species <- list()

for (sp in species_vec) {
  heatmaps_by_species[[sp]] <- make_species_heatmap_from_df(
    comp_df = comp_df_heat,
    species_name = sp,
    top_n = top_n,
    tax_label = heat_rank
  )
}

## 5. To re-draw a specific species heatmap, e.g. vannamei:
# heatmaps_by_species$`L.vannamei gut metagenome`
