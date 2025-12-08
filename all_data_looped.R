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




