### ===== AD DATASETS INTEGRATION SCRIPT - SEQUENTIAL VERSION =====
# # Load required libraries
cat("\nLoading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# Verify Seurat version
cat(paste("Seurat version:", packageVersion("Seurat"), "\n"))

# NO PARALLEL PROCESSING - Sequential only
options(future.globals.maxSize = 50 * 1024^3)  # 50GB limit just in case

# Function for log progress
log_progress <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
  flush.console()
}

# Function to report memory usage
report_memory <- function(step_name) {
  gc_result <- gc()
  used_mb <- sum(gc_result[,2])
  log_progress(paste("Memory after", step_name, ":", round(used_mb/1024, 1), "GB used"))
}


### ===== LOAD MATHYS_2020 ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork) # For combining plots
library(celldex)
library(SingleR)
library(plyr)
library(pheatmap)
library(HGNChelper)
library(Matrix)
data_dir <- "/pathway/to/data/dir"

# Read files 
matrix_path <- paste0(data_dir, "CellRangerOutput_matrix.mtx")
genes_path <- paste0(data_dir, "CellRangerOutput_genes.tsv")
barcodes_path <- paste0(data_dir, "CellRangerOutput_barcodes.tsv")

# Load the matrix
mat <- readMM(matrix_path)

# Load genes and barcodes
genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

# Set dimnames
rownames(mat) <- genes$V1
colnames(mat) <- barcodes$V1

# Create Seurat object
mathys <- CreateSeuratObject(counts = mat, 
                             project = "Mathys_AD",
                             min.cells = 3,
                             min.features = 200)
#QUALITY CONTROL
# Add mitochondrial percentage
mathys[["percent.mt"]] <- PercentageFeatureSet(mathys, pattern = "^MT-")
mathys <- subset(mathys, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

#LOAD AND ADD METADATA
sample_key <- read.csv(paste0(data_dir, "snRNAseqPFC_BA10_Sample_key.csv"))
clinical <- read.csv(paste0(data_dir, "ROSMAP_clinical.csv"))

#NORMALIZATION AND SCALING 
mathys <- NormalizeData(mathys)
mathys <- FindVariableFeatures(mathys, selection.method = "vst", nfeatures = 2000)
mathys <- ScaleData(mathys)
mathys <- RunPCA(mathys, features = VariableFeatures(object = mathys))
ElbowPlot(mathys, ndims = 50)

#CLUSTERING
mathys <- FindNeighbors(mathys, dims = 1:30)
mathys <- FindClusters(mathys, resolution = 0.8)
mathys <- RunUMAP(mathys, dims = 1:30)

saveRDS(mathys, file = "PFC_BA10_Mathys_AD_processed.rds")


### ===== LOAD MATHYS_2024####
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(celldex)
library(SingleR)
library(plyr)
library(pheatmap)
library(HGNChelper)
library(Matrix)

# Read the data
entorhinal <- readRDS("/Users/amoghchangavi/R_projects/Human Single Cell Integration/AD_snRNA_Hippocampus_Entorhinal_Mathys/Entorhinal_cortex.rds")
hippocampus <- readRDS("/Users/amoghchangavi/R_projects/Human Single Cell Integration/AD_snRNA_Hippocampus_Entorhinal_Mathys/Hippocampus.rds")

# Get count matrices
counts_ent <- GetAssayData(entorhinal, slot = "counts")
counts_hip <- GetAssayData(hippocampus, slot = "counts")

# Get metadata
meta_ent <- entorhinal@meta.data
meta_hip <- hippocampus@meta.data

# Clear original objects
rm(entorhinal, hippocampus)
gc()

# Add prefixes to cell names
colnames(counts_ent) <- paste0("EC_", colnames(counts_ent))
colnames(counts_hip) <- paste0("HC_", colnames(counts_hip))
rownames(meta_ent) <- paste0("EC_", rownames(meta_ent))
rownames(meta_hip) <- paste0("HC_", rownames(meta_hip))

# Add brain region to metadata
meta_ent$brain_region <- "Entorhinal_cortex"
meta_hip$brain_region <- "Hippocampus"

all_genes <- union(rownames(counts_ent), rownames(counts_hip))
genes_ent_only <- setdiff(rownames(counts_ent), rownames(counts_hip))
genes_hip_only <- setdiff(rownames(counts_hip), rownames(counts_ent))

# Create small zero matrices only for missing genes
if(length(genes_hip_only) > 0) {
  zeros_for_ent <- Matrix(0, 
                          nrow = length(genes_hip_only), 
                          ncol = ncol(counts_ent),
                          sparse = TRUE)
  rownames(zeros_for_ent) <- genes_hip_only
  colnames(zeros_for_ent) <- colnames(counts_ent)
  counts_ent <- rbind(counts_ent, zeros_for_ent)
  rm(zeros_for_ent)
}

if(length(genes_ent_only) > 0) {
  zeros_for_hip <- Matrix(0, 
                          nrow = length(genes_ent_only), 
                          ncol = ncol(counts_hip),
                          sparse = TRUE)
  rownames(zeros_for_hip) <- genes_ent_only
  colnames(zeros_for_hip) <- colnames(counts_hip)
  counts_hip <- rbind(counts_hip, zeros_for_hip)
  rm(zeros_for_hip)
}

# Reorder to ensure same gene order
counts_ent <- counts_ent[all_genes, ]
counts_hip <- counts_hip[all_genes, ]
combined_counts <- cbind(counts_ent, counts_hip)

# Clear intermediate objects
rm(counts_ent, counts_hip, all_genes, genes_ent_only, genes_hip_only)
gc()

# Combine metadata
combined_meta <- rbind(meta_ent, meta_hip)
rm(meta_ent, meta_hip)
gc()

# Create Seurat object
mathys_2024_core <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = combined_meta,
  project = "Mathys_2024_Core_AD_Regions",
  min.cells = 3,
  min.features = 200
)

rm(combined_counts, combined_meta)
gc()

mathys_2024_core[["percent.mt"]] <- PercentageFeatureSet(mathys_2024_core, pattern = "^MT-")
mathys_2024_core <- subset(mathys_2024_core, 
                           subset = nFeature_RNA > 200 & 
                             nFeature_RNA < 6000 & 
                             percent.mt < 20)




mathys_2024_core <- NormalizeData(mathys_2024_core)
mathys_2024_core <- FindVariableFeatures(mathys_2024_core, 
                                         selection.method = "vst", 
                                         nfeatures = 3000)
mathys_2024_core <- ScaleData(mathys_2024_core)
mathys_2024_core <- RunPCA(mathys_2024_core, features = VariableFeatures(mathys_2024_core))

#INTEGRATION
library(harmony)
mathys_2024_core <- RunHarmony(mathys_2024_core, 
                                 group.by.vars = "brain_region",
                                 plot_convergence = TRUE)
  
# Use harmony embeddings for downstream
reduction_to_use <- "harmony"
mathys_2024_core <- FindNeighbors(mathys_2024_core, 
                                  reduction = reduction_to_use, 
                                  dims = 1:30)
mathys_2024_core <- FindClusters(mathys_2024_core, resolution = 0.8)
mathys_2024_core <- RunUMAP(mathys_2024_core, 
                            reduction = reduction_to_use, 
                            dims = 1:30)

saveRDS(mathys_2024_core, file = "mathys_2024_core_processed.rds")


### ===== LOAD MORABITO ####
data_dir <- "path/to/data/dir"

# Read files 
matrix_path <- paste0(data_dir, "snRNA_counts.mtx")
genes_path <- paste0(data_dir, "genes.csv")
barcodes_path <- paste0(data_dir, "barcodes_rna.csv")

# Load the matrix
mat <- readMM(matrix_path)

# Load genes and barcodes
genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

# Set dimnames
rownames(mat) <- genes$V1
colnames(mat) <- barcodes$V1

# Create Seurat object
uci_seurat <- CreateSeuratObject(counts = mat, 
                                 project = "UCI_Multiomics",
                                 min.cells = 3,
                                 min.features = 200)
metadata <- read.csv("snRNA_metadta.csv", row.names = 1) 


# If the cell names match perfectly, add metadata directly
if(all(rownames(metadata) %in% colnames(uci_seurat))) {
  uci_seurat <- AddMetaData(uci_seurat, metadata = metadata)
} else {
  print("Cell name mismatch detected. Checking patterns...")
  common_cells <- intersect(rownames(metadata), colnames(uci_seurat))
  print(paste("Common cells found:", length(common_cells)))
  
  # Subset to common cells
  metadata_matched <- metadata[common_cells, ]
  uci_seurat <- uci_seurat[, common_cells]
  
  # Add metadata
  uci_seurat <- AddMetaData(uci_seurat, metadata = metadata_matched)
}

# Add the pre-computed UMAP coordinates
umap_coords <- as.matrix(metadata[, c("UMAP_1", "UMAP_2")])
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Create a DimReduc object for the UMAP
uci_seurat[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, 
                                             key = "UMAP_", 
                                             assay = DefaultAssay(uci_seurat))

# Check the data
table(uci_seurat$Diagnosis)
table(uci_seurat$celltype)
table(uci_seurat$cluster)


# Pre-computed annotations exist, but we still run standard workflow for integration with other datasets
uci_seurat <- NormalizeData(uci_seurat)
uci_seurat <- FindVariableFeatures(uci_seurat, selection.method = "vst", nfeatures = 2000)
uci_seurat <- ScaleData(uci_seurat)
uci_seurat <- RunPCA(uci_seurat)
ElbowPlot(uci_seurat, ndims = 50)

# Add dataset identifier
uci_seurat$dataset <- "UCI_Multiomics"

saveRDS(uci_seurat, file = "uci_multiomics_processed.rds")


### ===== LOAD FUJITA ####
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(patchwork)

#LOAD H5SEURAT OBJECT 
log_progress("Loading microglia.h5Seurat object...")
microglia <- LoadH5Seurat("microglia.h5Seurat")


print("\nFirst few rows of metadata:")
print(head(microglia@meta.data))

#ADD FUJITA-STANDARD METADATA
log_progress("Adding Fujita-standard metadata fields...")

# Dataset identifier
microglia$dataset_name <- "microglia_dataset"

# Add DLPFC brain region 
microglia$brain_region <- "DLPFC"

# Add basic QC metrics following Fujita standards
microglia$nCount_RNA <- microglia$nUMI  
microglia$nFeature_RNA <- microglia$nGene  

# Calculate mitochondrial gene percentage
microglia$percent.mt <- PercentageFeatureSet(microglia, pattern = "^MT-")


#QUALITY CONTROL 
log_progress("Applying QC filters...")
microglia <- subset(microglia, subset = nFeature_RNA > 200 & 
                      nFeature_RNA < 5000 & 
                      percent.mt < 20)

log_progress(paste("After QC filtering:", ncol(microglia), "cells remain"))

#NORMALIZATION/PCA
microglia <- SCTransform(microglia, 
                         variable.features.n = 2000,
                         conserve.memory = TRUE,
                         verbose = TRUE)
microglia <- RunPCA(microglia, npcs = 30, verbose = TRUE)
optimal_pcs <- 15

microglia <- FindNeighbors(microglia, dims = 1:optimal_pcs)
microglia <- FindClusters(microglia, 
                          resolution = 0.2,  
                          algorithm = 1)# Louvain algorithm

#UMAP GENERATION 
log_progress("Generating UMAP...")
microglia <- RunUMAP(microglia, dims = 1:optimal_pcs)
saveRDS(microglia, file = file.path(output_dir, "microglia_fujita_processed.rds"))


### ===== LOAD BLANCHARD ####
# Load required libraries
library(Matrix)
library(Seurat)
library(dplyr)

# Set data path
data_path <- "/Users/amoghchangavi/SavioMount/AD_snRNAseq_Prefrontal_Blanchard"

# Load the count matrix/ gene names
counts <- readMM(file.path(data_path, "qc_counts.mtx"))
gene_names <- read.table(file.path(data_path, "qc_gene_names.txt"), 
                         header = FALSE, 
                         stringsAsFactors = FALSE)
colnames(gene_names) <- "gene_name"

# Load column metadata
cell_metadata <- read.csv(file.path(data_path, "qc_column_metadata.csv"), 
                          row.names = 1, 
                          stringsAsFactors = FALSE)

# Assign row and column names to count matrix
rownames(counts) <- gene_names$gene_name
colnames(counts) <- rownames(cell_metadata)

# Create Seurat object
blanchard_seurat <- CreateSeuratObject(
  counts = counts,
  meta.data = cell_metadata,
  project = "Blanchard_APOE4"
)
blanchard_seurat$dataset <- "blanchard"

# Save as RDS file
output_file <- file.path(data_path, "blanchard_apoe4_seurat.rds")
cat(sprintf("Saving Seurat object to: %s\n", output_file))
saveRDS(blanchard_seurat, file = output_file)


### ===== LOAD LI (Slurm Script for HPC cluster submission) ####
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

#FUNCTION: PROCESS LI TEMPORAL WITH REAL METADATA
process_li_temporal_with_real_metadata <- function(data_dir, metadata_csv_path) {
  
  # Load metadata from CSV file
  metadata <- read.csv(metadata_csv_path, stringsAsFactors = FALSE)
  
  # Save metadata for reference
  metadata_path <- file.path(dirname(data_dir), "li_temporal_metadata_loaded.csv")
  cat("Metadata saved to:", metadata_path, "\n")
  
  # Find all tar.gz files
  tar_files <- list.files(data_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)
  cat("Found", length(tar_files), "tar.gz files\n")
  
  # Extract all files
  extract_dir <- file.path(data_dir, "extracted")
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  
  seurat_list <- list()
  failed_samples <- character()
  
  for (i in seq_along(tar_files)) {
    tar_file <- tar_files[i]
    cat("Processing", i, "of", length(tar_files), ":", basename(tar_file), "\n")
    sample_id <- gsub(".*_(GEX_\\d+R)\\.tar\\.gz$", "\\1", basename(tar_file))
    
    tryCatch({
      # Create sample directory
      sample_dir <- file.path(extract_dir, sample_id)
      dir.create(sample_dir, showWarnings = FALSE)
      
      # Extract tar.gz file
      untar(tar_file, exdir = sample_dir)
      
      # Find the 10X files
      matrix_file <- list.files(sample_dir, pattern = "matrix", full.names = TRUE, recursive = TRUE)[1]
      barcodes_file <- list.files(sample_dir, pattern = "barcodes", full.names = TRUE, recursive = TRUE)[1]
      features_file <- list.files(sample_dir, pattern = "features", full.names = TRUE, recursive = TRUE)[1]
      
      if (is.na(matrix_file) || is.na(barcodes_file) || is.na(features_file)) {
        cat("Warning: Missing files for sample", sample_id, "\n")
        failed_samples <- c(failed_samples, sample_id)
        next
      }
      
      # Load 10X data
      sample_data <- Read10X(data.dir = dirname(matrix_file))
      
      # Create Seurat object
      seurat_obj <- CreateSeuratObject(
        counts = sample_data,
        project = "Li_Temporal",
        min.cells = 3,
        min.features = 200
      )
      
      # Get metadata for this sample
      sample_meta <- metadata[metadata$sample_id == sample_id, ]
      
      if (nrow(sample_meta) == 0) {
        cat("Warning: No metadata found for sample", sample_id, "\n")
        failed_samples <- c(failed_samples, sample_id)
        next
      }
      
      # Add sample-level metadata
      seurat_obj$sample_id <- sample_id
      seurat_obj$apoe_genotype <- sample_meta$apoe_genotype[1]
      seurat_obj$apoe_category <- sample_meta$apoe_category[1]
      seurat_obj$pathology_group <- sample_meta$pathology_group[1]
      seurat_obj$sex <- sample_meta$sex[1]
      seurat_obj$brain_region <- "Temporal_Cortex"
      seurat_obj$dataset_name <- "li"
      
      # Add mitochondrial and ribosomal gene percentages
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
      
      # Make cell names unique by adding sample prefix
      seurat_obj <- RenameCells(seurat_obj, add.cell.id = sample_id)
      seurat_list[[sample_id]] <- seurat_obj
  
      cat("Sample", sample_id, "processed:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes\n")
      
    }, error = function(e) {
      cat("ERROR processing sample", sample_id, ":", e$message, "\n")
      failed_samples <- c(failed_samples, sample_id)
    })
    
    # Memory cleanup after each sample
    gc()
  }
  
  if (length(failed_samples) > 0) {
    cat("Failed samples:", paste(failed_samples, collapse = ", "), "\n")
  }
  if (length(seurat_list) == 0) {
    stop("No samples processed successfully!")
  }
  
  # Merge all samples
  cat("Merging", length(seurat_list), "samples...\n")
  li_temporal_merged <- merge(seurat_list[[1]], y = seurat_list[-1])
  
  # Clear individual objects to save memory
  rm(seurat_list)
  gc()
  
  # Basic QC filtering
  cat("Applying QC filters...\n")
  cat("Before filtering:", ncol(li_temporal_merged), "cells\n")
  
  li_temporal_filtered <- subset(li_temporal_merged, 
                                 subset = nFeature_RNA > 500 & 
                                   nFeature_RNA < 7000 & 
                                   percent.mt < 20)
  
  cat("After QC filtering:", ncol(li_temporal_filtered), "cells\n")
  
  # Clear merged object to save memory
  rm(li_temporal_merged)
  gc()
  
  # Standard processing
  cat("Running normalization and scaling...\n")
  li_temporal_filtered <- NormalizeData(li_temporal_filtered, verbose = FALSE)
  li_temporal_filtered <- FindVariableFeatures(li_temporal_filtered, nfeatures = 3000, verbose = FALSE)
  li_temporal_filtered <- ScaleData(li_temporal_filtered, verbose = FALSE)
  li_temporal_filtered <- RunPCA(li_temporal_filtered, npcs = 50, verbose = FALSE)
  
 
  # Save processed object
  output_path <- file.path(output_dir, "li_temporal_cortex_processed_real_metadata.rds")
  cat("Saving processed object to:", output_path, "\n")
  saveRDS(li_temporal_filtered, output_path)
  
  summary_path <- file.path(output_dir, "li_temporal_processing_summary_real_metadata.txt")
  writeLines(summary_text, summary_path)
  cat(summary_text)
  
  return(li_temporal_filtered)
}

#MAIN EXECUTION
main <- function() {
  cat("=== STARTING LI TEMPORAL CORTEX PROCESSING ===\n")
  
  # Set paths
  data_dir <- "/path/to/data/dir"
  metadata_csv_path <- "/path/to/metadata"
  
  # Check if paths exist
  if (!dir.exists(data_dir)) {
    stop("Data directory not found: ", data_dir)
  }
  
  if (!file.exists(metadata_csv_path)) {
    stop("Metadata CSV file not found: ", metadata_csv_path)
  }
  
  # Process the dataset
  li_temporal <- process_li_temporal_with_real_metadata(data_dir, metadata_csv_path)
  
  cat("=== PROCESSING COMPLETED SUCCESSFULLY ===\n")
  
  return(li_temporal)
}

# Execute main function
if (!interactive()) {
  main()
}


### ===== LOAD PROCESSED DATASETS =====
log_progress("Loading processed datasets...")

# Define paths - looking for any .rds file in each directory
find_rds <- function(dir_path) {
  if(!dir.exists(dir_path)) {
    stop(paste("Directory not found:", dir_path))
  }

  rds_files <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if(length(rds_files) == 0) {
    stop(paste("No .rds file found in", dir_path))
  }

  if(length(rds_files) > 1) {
    processed <- grep("processed", rds_files, value = TRUE)
    if(length(processed) > 0) {
      return(processed[1])
    }
  }

  return(rds_files[1])
}

# Find .rds files in each directory
datasets <- list(
  mathys_2024 = find_rds("/pathway/to/rds/object"),
  morabito = find_rds("/pathway/to/rds/object"),
  mathys_2020 = find_rds("/pathway/to/rds/object"),
  li=find_rds("/pathway/to/rds/object"),
  fujita=find_rds("/pathway/to/rds/object"),
  blanchard=find_rds("/pathway/to/rds/object")
)

# Load each dataset with error handling
seurat_list <- list()
for(name in names(datasets)) {
  log_progress(paste("Loading", name, "from", datasets[[name]]))

  tryCatch({
    obj <- readRDS(datasets[[name]])
    obj$dataset_name <- name
    log_progress(paste(name, "loaded:", ncol(obj), "cells,", nrow(obj), "genes"))
    seurat_list[[name]] <- obj

    # Clean up
    rm(obj)
    gc()
  }, error = function(e) {
    log_progress(paste("ERROR loading", name, ":", e$message))
    stop(e)
  })
}


### ===== CONVERT MATHYS ENSEMBL IDS =====
library(biomaRt)

# Get rownames (Ensembl IDs)
ensembl_ids <- rownames(seurat_list[["mathys_2020"]])

# Use Ensembl BioMart to map
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Clean up
rm(ensembl_ids)
gc()

# Remove entries with empty symbols
mapping <- mapping[mapping$hgnc_symbol != "", ]

# Build lookup table
id_map <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
rm(mapping)
gc()

# Subset Seurat object to mappable genes
obj <- seurat_list[["mathys_2020"]]
keep_genes <- intersect(rownames(obj), names(id_map))
obj <- subset(obj, features = keep_genes)

# Rename genes to gene symbols
rownames(obj) <- id_map[rownames(obj)]

# Save back to list
seurat_list[["mathys_2020"]] <- obj

# Clean up
rm(obj, id_map, keep_genes)
gc()

log_progress("Mathys_2020 gene names mapped to HGNC symbols.")


### ===== ADDING PATIENT METADATA FOR MATHYS2020 =====
log_progress("Applying special pre-processing for mathys_2020...")

# Map the numeric sample suffix to a projid
mathys_sample_key_path <- read.csv("snRNAseqPFC_BA10_Sample_key.csv")
if (file.exists(mathys_sample_key_path)) {
  tryCatch({
    mathys_sample_key <- read.csv(mathys_sample_key_path)
    required_key_cols <- c("projid", "sample")
    if (!all(required_key_cols %in% colnames(mathys_sample_key))) {
      missing_cols <- required_key_cols[!required_key_cols %in% colnames(mathys_sample_key)]
      stop(sprintf("The mathys_sample_key.csv file is missing required column(s): %s. Found columns: %s",
                   paste(missing_cols, collapse = ", "),
                   paste(colnames(mathys_sample_key), collapse = ", ")))
    }

    mathys_obj <- seurat_list[['mathys_2020']]

    # Create a data frame for the new metadata
    metadata_to_add <- data.frame(row.names = colnames(mathys_obj))
    sample_numbers <- gsub(".*-", "", rownames(metadata_to_add))

    # Create a lookup table: sample number -> projid
    sample_to_projid <- setNames(as.character(mathys_sample_key$projid), as.character(mathys_sample_key$sample))

    # Map the sample numbers to projids
    metadata_to_add$projid <- sample_to_projid[sample_numbers]

    # Use AddMetaData to safely add the new column to the object
    seurat_list[['mathys_2020']] <- AddMetaData(mathys_obj, metadata = metadata_to_add)

    rm(mathys_obj, mathys_sample_key, metadata_to_add)

  }, error = function(e) {
    log_progress(paste("ERROR during Mathys 2020 pre-processing:", e$message))
    stop(e)
  })
} else {
  log_progress(paste("WARNING: Mathys 2020 sample key file not found at:", mathys_sample_key_path, ". Patient IDs will be missing."))
}


### ===== SEX METADATA MAPPING =====
log_progress("Establishing ID mapping...")
log_progress(paste(rep("-", 30), collapse = ""))

seurat_meta_cols <- colnames(integrated@meta.data)

if("projid" %in% seurat_meta_cols) {
  id_col <- "projid"
} else if("individualID" %in% seurat_meta_cols) {
  id_col <- "individualID"
} else if("donor" %in% seurat_meta_cols) {
  id_col <- "donor"
} else {
  log_progress("No direct ID match found. Examining sample metadata values...")
  for(col in head(seurat_meta_cols, 10)) {
    sample_vals <- head(unique(integrated@meta.data[[col]]), 5)
    log_progress(sprintf("Column '%s': %s", col, paste(sample_vals, collapse = ", ")))
  }
  stop("Could not find ID mapping column. Please check metadata structure and update script.")
}

log_progress(sprintf("Using ID column: %s", id_col))

seurat_ids <- unique(integrated@meta.data[[id_col]])
clinical_ids <- unique(clinical_data$projid)
overlapping_ids <- intersect(seurat_ids, clinical_ids)


#COMPREHENSIVE SEX DETECTION STRATEGY 

#Check existing Sex column
existing_sex_count <- 0
if("Sex" %in% colnames(integrated@meta.data)) {
  log_progress("✓ Found existing Sex column")

  sex_values <- table(integrated$Sex, useNA = "ifany")
  log_progress("Current Sex column values:")
  print(sex_values)

  integrated$Sex_Standardized <- case_when(
    integrated$Sex %in% c("M", "Male", "male", "MALE", "m") ~ "Male",
    integrated$Sex %in% c("F", "Female", "female", "FEMALE", "f") ~ "Female",
    integrated$Sex == 1 ~ "Male",
    integrated$Sex == 0 ~ "Female",
    TRUE ~ "Unknown"
  )

  existing_sex_count <- sum(integrated$Sex_Standardized %in% c("Male", "Female"))
  log_progress(sprintf("Existing Sex column covers %d cells after standardization", existing_sex_count))

} else {
  log_progress("No existing Sex column found")
  integrated$Sex_Standardized <- "Unknown"
}

#Check for OTHER sex-related columns
log_progress("Checking for additional sex-related columns...")

possible_sex_cols <- c("sex", "gender", "Gender", "GENDER", "Sex_standardized", "sex_standardized","msex","SEX")
found_additional_sex <- FALSE

for(col in possible_sex_cols) {
  if(col %in% colnames(integrated@meta.data) && col != "Sex") {
    log_progress(sprintf("Found additional sex column: %s", col))

    col_values <- table(integrated@meta.data[[col]], useNA = "ifany")
    log_progress(sprintf("Values in %s:", col))
    print(col_values)

    missing_indices <- which(integrated$Sex_Standardized == "Unknown")

    # Standardize the additional sex column first
    additional_sex <- case_when(
      integrated@meta.data[[col]] %in% c("M", "Male", "male", "MALE", "m", 1) ~ "Male",
      integrated@meta.data[[col]] %in% c("F", "Female", "female", "FEMALE", "f", 0) ~ "Female",
      TRUE ~ "Unknown"
    )

    # Use a single, vectorized ifelse to fill in missing values
    integrated$Sex_Standardized <- ifelse(
      integrated$Sex_Standardized == "Unknown" & additional_sex != "Unknown",
      additional_sex,
      integrated$Sex_Standardized
    )
    found_additional_sex <- TRUE
    new_count <- sum(integrated$Sex_Standardized %in% c("Male", "Female"))
    additional_contribution <- new_count - existing_sex_count
    log_progress(sprintf("✓ Column %s contributed %d additional cells", col, additional_contribution))
    existing_sex_count <- new_count
  }
}

if(!found_additional_sex) {
  log_progress("No additional sex columns found")
}

#MATHYS 2020 BARCODE SUFFIX MAPPING
log_progress("=== APPLYING MATHYS 2020 BARCODE SUFFIX MAPPING ===")

mathys_sex_count <- 0
if(!is.null(mathys_sample_key) && "dataset_name" %in% colnames(integrated@meta.data)) {

  # Identify Mathys 2020 cells
  mathys_identifiers <- c("mathys_2020", "Mathys_2020", "mathys2020", "AD_snRNAseqPFC_BA10_Mathys", "Mathys_AD")
  mathys_cells_indices <- NULL
  mathys_dataset_name <- NULL

  for(identifier in mathys_identifiers) {
    if(identifier %in% integrated$dataset_name) {
      mathys_cells_indices <- which(integrated$dataset_name == identifier)
      mathys_dataset_name <- identifier
      break
    }
  }

  if(!is.null(mathys_cells_indices) && length(mathys_cells_indices) > 0) {
    log_progress(sprintf("Found %d Mathys 2020 cells with identifier '%s'",
                         length(mathys_cells_indices), mathys_dataset_name))

    # Extract cell barcodes and sample numbers
    mathys_cell_barcodes <- colnames(integrated)[mathys_cells_indices]
    sample_numbers <- as.numeric(gsub(".*-", "", mathys_cell_barcodes))
    unique_sample_numbers <- unique(sample_numbers)

    log_progress(sprintf("Extracted sample numbers: %s (total: %d unique)",
                         paste(head(sort(unique_sample_numbers), 10), collapse = ", "),
                         length(unique_sample_numbers)))

    # Verify mapping exists
    sample_key_samples <- unique(mathys_sample_key$sample)
    matching_samples <- intersect(unique_sample_numbers, sample_key_samples)

    if(length(matching_samples) > 0) {
      log_progress(sprintf("Sample mapping available: %d/%d samples match",
                           length(matching_samples), length(unique_sample_numbers)))

      # Create mappings
      sample_to_projid <- setNames(mathys_sample_key$projid, mathys_sample_key$sample)
      mathys_projids <- sample_to_projid[as.character(sample_numbers)]

      # Get clinical data
      clinical_subset <- clinical_data[clinical_data$projid %in% mathys_projids, ]
      projid_to_clinical <- clinical_subset
      rownames(projid_to_clinical) <- as.character(clinical_subset$projid)

      # Apply sex mapping to Mathys cells
      for(i in seq_along(mathys_cells_indices)) {
        cell_idx <- mathys_cells_indices[i]
        projid <- mathys_projids[i]

        if(integrated$Sex_Standardized[cell_idx] == "Unknown" &&
           !is.na(projid) && as.character(projid) %in% rownames(projid_to_clinical)) {

          msex_val <- projid_to_clinical[as.character(projid), "msex"]
          if(!is.na(msex_val)) {
            if(msex_val == 1) {
              integrated$Sex_Standardized[cell_idx] <- "Male"
              mathys_sex_count <- mathys_sex_count + 1
            } else if(msex_val == 0) {
              integrated$Sex_Standardized[cell_idx] <- "Female"
              mathys_sex_count <- mathys_sex_count + 1
            }
          }
        }
      }

      log_progress(sprintf("MATHYS BARCODE MAPPING: Added sex data for %d cells", mathys_sex_count))

      existing_sex_count <- sum(integrated$Sex_Standardized %in% c("Male", "Female"))
    }
  }
}


### ===== ROSMAP METADATA MAPPING FOR REMAINING UNKNOWNS ====
log_progress("Adding standard ROSMAP mapping for remaining unknowns...")

cell_metadata <- integrated@meta.data
cell_metadata$cell_barcode <- rownames(cell_metadata)

clinical_subset <- clinical_data[clinical_data$projid %in% overlapping_ids, ]
merged_metadata <- merge(cell_metadata, clinical_subset,
                         by.x = id_col, by.y = "projid", all.x = TRUE)

log_progress(sprintf("Standard ROSMAP mapping available for %d cells", sum(!is.na(merged_metadata$msex))))

rosmap_sex <- ifelse(merged_metadata$msex == 1, "Male",
                     ifelse(merged_metadata$msex == 0, "Female", "Unknown"))

merged_metadata_ordered <- merged_metadata[match(rownames(integrated@meta.data),
                                                 merged_metadata$cell_barcode), ]
rosmap_sex_ordered <- ifelse(merged_metadata_ordered$msex == 1, "Male",
                             ifelse(merged_metadata_ordered$msex == 0, "Female", "Unknown"))

missing_sex_indices <- which(integrated$Sex_Standardized == "Unknown")
log_progress(sprintf("Filling %d remaining missing sex values using standard ROSMAP data...", length(missing_sex_indices)))

rosmap_contribution <- 0
# Identify which cells are currently missing sex data
is_missing <- integrated$Sex_Standardized == "Unknown"
can_be_filled <- !is.na(rosmap_sex_ordered) & rosmap_sex_ordered != "Unknown"

# Calculate contribution before the update
indices_to_update <- which(is_missing & can_be_filled)

# Apply the update in a single vectorized operation
if (length(indices_to_update) > 0) {
  integrated$Sex_Standardized[indices_to_update] <- rosmap_sex_ordered[indices_to_update]
}

# Now, create the APOE category using a new, reliable column.
merged_metadata_ordered$APOE_Final <- case_when(
  merged_metadata_ordered$apoe_genotype_final %in% c(22, 23, 33) ~ "No_E4",
  merged_metadata_ordered$apoe_genotype_final %in% c(24, 34) ~ "One_E4",
  merged_metadata_ordered$apoe_genotype_final == 44 ~ "Two_E4",
  TRUE ~ "Unknown"
)
integrated$APOE_Final <- merged_metadata_ordered$APOE_Final
log_progress("Created the 'APOE_Category' column.")


### ===== REMOVING DUPLICATES BEFORE INTEGRATION =====
log_progress("De-duplication Step 1: Standardizing patient IDs...")

seurat_list <- lapply(seurat_list, function(obj) {
  potential_id_cols <- c("individualID", "projid", "specimenID", "Sample.ID", "sample_id")
  cols_to_use <- intersect(potential_id_cols, colnames(obj@meta.data))
  
  if (length(cols_to_use) == 0) {
    obj$unified_patient_id <- NA_character_
    log_progress(paste("WARNING: No potential patient ID columns found in",
                       obj$dataset_name[1], "- cells cannot be de-duplicated."))
    return(obj)
  }
  
  # 1. Isolate only the potential ID columns and add unified ID back to the Seurat object. 
  id_df <- obj@meta.data[, cols_to_use, drop = FALSE]
  id_df[] <- lapply(id_df, as.character)
  unified_id_vector <- do.call(dplyr::coalesce, id_df)
  obj$unified_patient_id <- unified_id_vector
  
  return(obj)
})

#Creating a master metadata table
all_metadata <- lapply(names(seurat_list), function(name) {
  meta <- seurat_list[[name]]@meta.data
  return(meta[, c("unified_patient_id", "dataset_name")])
}) %>%
  bind_rows() %>%
  filter(!is.na(unified_patient_id))

dataset_priority <- c("mathys_2024", "li", "fujita", "morabito", "mathys_2020", "blanchard")

duplicate_patients <- all_metadata %>%
  group_by(unified_patient_id) %>%
  filter(n_distinct(dataset_name) > 1) %>%
  ungroup()

if (nrow(duplicate_patients) > 0) {
  datasets_to_keep <- duplicate_patients %>%
    mutate(dataset_name = factor(dataset_name, levels = dataset_priority)) %>%
    group_by(unified_patient_id) %>%
    arrange(dataset_name) %>%
    slice(1) %>%
    ungroup()
  
  removals <- duplicate_patients %>%
    anti_join(datasets_to_keep, by = c("unified_patient_id", "dataset_name"))
  
  if (nrow(removals) > 0) {
    log_progress(paste("Found", n_distinct(removals$unified_patient_id), "duplicate patients to handle."))
    print(removals %>% count(dataset_name, unified_patient_id))
  }
} else {
  removals <- data.frame(unified_patient_id = character(), dataset_name = character()) # Create empty frame if no duplicates
}


#Filtering Seurat objects 
for (name in names(seurat_list)) {
  patients_to_remove_from_this_dataset <- removals %>%
    filter(dataset_name == name) %>%
    pull(unified_patient_id)
  
  if (length(patients_to_remove_from_this_dataset) > 0) {
    log_progress(paste("In", name, ", flagging cells from",
                       length(patients_to_remove_from_this_dataset), "duplicate patients for removal."))
    
    patient_ids_in_object <- seurat_list[[name]]$unified_patient_id
    cells_to_remove_mask <- patient_ids_in_object %in% patients_to_remove_from_this_dataset
    cells_to_remove_mask[is.na(cells_to_remove_mask)] <- FALSE
    cells_to_keep <- colnames(seurat_list[[name]])[!cells_to_remove_mask]
    
    original_cell_count <- ncol(seurat_list[[name]])
    seurat_list[[name]] <- subset(seurat_list[[name]], cells = cells_to_keep)
    
    log_progress(paste("Filtered", name, ": removed",
                       original_cell_count - ncol(seurat_list[[name]]), "cells. New count:",
                       ncol(seurat_list[[name]])))
  } else {
    log_progress(paste("No duplicate patients to remove from", name))
  }
}


### ===== NORMALIZE/VARIABLE FEATURES/MERGE INTEGRATED OBJECT =====
for (name in names(seurat_list)) {
  log_progress(paste("Processing", name))
  obj <- seurat_list[[name]]
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  seurat_list[[name]] <- obj
  
  # Clean up
  rm(obj)
  gc()
}

#Add dataset identifiers
log_progress("Adding dataset identifiers to each object...")

for(name in names(seurat_list)) {
  if("dataset_name" %in% colnames(seurat_list[[name]]@meta.data)) {
    log_progress(paste(name, "already has dataset_name"))
  } else {
    log_progress(paste(name, "missing dataset_name - adding it"))
    seurat_list[[name]]$dataset_name <- name
  }
}


integrated <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "AD_Integration"
)

log_progress(paste("Merged object:", nrow(integrated), "genes,", ncol(integrated), "cells"))

# Clean up
rm(seurat_list)


### ===== CREATE BRAIN REGION VARIABLE =====
log_progress("Checking existing metadata columns:")
region_cols <- grep("region|area|tissue|brain|cortex|hippocampus|entorhinal",
                    colnames(integrated@meta.data),
                    value = TRUE, ignore.case = TRUE)

# Initialize brain_region column
integrated$brain_region <- "Unknown"

# Assign regions based on dataset using base R indexing (bulletproof approach)
integrated$brain_region[integrated$dataset_name == "morabito"] <- "Prefrontal_Cortex"
integrated$brain_region[integrated$dataset_name == "mathys_2020"] <- "Prefrontal_Cortex"
integrated$brain_region[integrated$dataset_name == "mathys_2024"] <- "Hippocampus_Entorhinal"
integrated$brain_region[integrated$dataset_name == "li"] <- "Temporal_Cortex"
integrated$brain_region[integrated$dataset_name == "fujita"] <- "Prefrontal_Cortex"
integrated$brain_region[integrated$dataset_name == "blanchard"] <- "Prefrontal_Cortex"

#MATHYS brain region assignment
meta_cols <- colnames(integrated@meta.data)

# Look for any preserved orig.ident or brain_region from original object
if("orig.ident" %in% meta_cols) {
  mathys_orig_ident <- integrated$orig.ident[mathys_cells]
  orig_table <- table(mathys_orig_ident)

  if(!"brain_region" %in% colnames(integrated@meta.data)) {
    integrated$brain_region <- "Unknown"
    log_progress("Created new brain_region column")
  } else {
    log_progress("brain_region column already exists, updating...")
  }


  # Find mathys_2024 cells with HC orig.ident
  mathys_hc_indices <- which(integrated$dataset_name == "mathys_2024" &
                               integrated$orig.ident == "HC")

  # Find mathys_2024 cells with EC orig.ident
  mathys_ec_indices <- which(integrated$dataset_name == "mathys_2024" &
                               integrated$orig.ident == "EC")


  # Assign brain regions
  if(length(mathys_hc_indices) > 0) {
    integrated$brain_region[mathys_hc_indices] <- "Hippocampus"
    log_progress(paste("✓ Assigned", length(mathys_hc_indices), "HC cells to Hippocampus"))
  }
  if(length(mathys_ec_indices) > 0) {
    integrated$brain_region[mathys_ec_indices] <- "Entorhinal_Cortex"
    log_progress(paste("✓ Assigned", length(mathys_ec_indices), "EC cells to Entorhinal_Cortex"))
  }

  # Clean up indices
  rm(mathys_hc_indices, mathys_ec_indices)
  gc()

} else {
  log_progress("orig.ident column not found - cannot assign brain regions")
}


### ===== STANDARD PROCESSING FOR INTEGRATED OBJECT =====
integrated <- NormalizeData(integrated, verbose = FALSE)
integrated <- FindVariableFeatures(integrated, nfeatures = 3000, verbose = FALSE)
integrated <- ScaleData(integrated, features = VariableFeatures(integrated), verbose = TRUE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

#HARMONY INTEGRATION (Correct for dataset batch effects)
integrated <- RunHarmony(
  object = integrated,
  group.by = "dataset_name",
  reduction.use = "pca",
  dims = 1:30,
  plot_convergence = TRUE,
  verbose = TRUE
)

#DOWNSTREAM ANALYSIS 
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.6, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30, verbose = FALSE)


### ===== IDENTIFY MYELOID CLUSTERS =====
log_progress("Identifying microglia/macrophage clusters (1, 19, 28, 35 at res 0.6)...")

res_06_cols <- colnames(integrated@meta.data)[grepl("0.6|res.0.6|resolution.0.6|snn_res.0.6|RNA_snn_res.0.6|integrated_snn_res.0.6",
                                                    colnames(integrated@meta.data))]
log_progress(sprintf("Potential res 0.6 columns: %s", paste(res_06_cols, collapse = ", ")))

if(length(res_06_cols) == 0) {
  cluster_cols <- colnames(integrated@meta.data)[grepl("cluster|res|resolution", colnames(integrated@meta.data), ignore.case = TRUE)]
  log_progress(sprintf("Available clustering columns: %s", paste(cluster_cols, collapse = ", ")))
  stop("No resolution 0.6 clustering column found!")
}

cluster_col <- res_06_cols[1]
log_progress(sprintf("Using cluster column: %s", cluster_col))

available_clusters <- unique(integrated@meta.data[[cluster_col]])
log_progress(sprintf("Available clusters: %s", paste(sort(as.character(available_clusters)), collapse = ", ")))

target_clusters <- c("1", "19","28","35")
is_myeloid_cluster <- integrated@meta.data[[cluster_col]] %in% target_clusters
myeloid_indices <- which(is_myeloid_cluster)

if(length(myeloid_indices) == 0) {
  target_clusters_num <- c(1,19,28,35)
  is_myeloid_cluster <- integrated@meta.data[[cluster_col]] %in% target_clusters_num
  myeloid_indices <- which(is_myeloid_cluster)
}

log_progress(sprintf("Found %d cells in clusters 1, 19, 28,35", length(myeloid_indices)))

if(length(myeloid_indices) == 0) {
  log_progress("Cluster breakdown:")
  cluster_table <- table(integrated@meta.data[[cluster_col]])
  print(head(cluster_table, 20))
  stop("No cells found in target clusters 1, 19, 28, 35")
}

myeloid_subset <- integrated[, myeloid_indices]
log_progress(sprintf("Microglia/macrophage cluster subset: %d cells", ncol(myeloid_subset)))

myeloid_subset$Target_Cluster <- as.character(myeloid_subset@meta.data[[cluster_col]])

saveRDS(myeloid_subset, file.path(output_dir, "myeloid_subset_comprehensive_annotated_SEXONLY.rds"))


### ===== LOAD DATA AND FIX MYELOID SUBSET METADATA =====
log_progress("Loading myeloid_subset data...")

myeloid_subset <- readRDS("myeloid_subset_comprehensive_annotated_SEXONLY.rds")
ROSMAP <- read.csv("ROSMAP_clinical.csv")
output_dir <- "/path/to/output/dir"

cell_count <- ncol(myeloid_subset)
gene_count <- nrow(myeloid_subset)
log_progress(sprintf("Dataset: %d genes × %d cells", gene_count, cell_count))
colnames(myeloid_subset@meta.data)

# --- Step 1: Prepare the Metadata ---
original_metadata <- myeloid_subset@meta.data
original_metadata$cell_barcode <- rownames(original_metadata)

# --- Step 2: Perform a Safe Join ---
ROSMAP$projid <- as.character(ROSMAP$projid)
new_metadata <- left_join(original_metadata, 
                          select(ROSMAP, -any_of("individualID")), 
                          by = "projid")


# --- Step 3: Intelligently Combine Duplicate Columns ---
conflicting_cols <- sub("\\.x$", "", colnames(new_metadata)[grepl("\\.x$", colnames(new_metadata))])

log_progress(paste("Found", length(conflicting_cols), "conflicting columns to merge:", paste(conflicting_cols, collapse=", ")))

for (col in conflicting_cols) {
  col_x <- paste0(col, ".x")
  col_y <- paste0(col, ".y")
  new_metadata[[col_x]] <- as.character(new_metadata[[col_x]])
  new_metadata[[col_y]] <- as.character(new_metadata[[col_y]])
  
  # Use coalesce() to take the first non-NA value.
  new_metadata[[col]] <- coalesce(new_metadata[[col_x]], new_metadata[[col_y]])
  new_metadata[[col_x]] <- NULL
  new_metadata[[col_y]] <- NULL
}

existing_sex_count <- sum(new_metadata$Sex_Standardized %in% c("Male", "Female"))
possible_sex_cols <- c("sex", "gender", "Gender", "GENDER", "Sex_standardized", "sex_standardized","msex","SEX")
found_additional_sex <- FALSE

for(col in possible_sex_cols) {
  if(col %in% colnames(new_metadata) && col != "Sex") {
    log_progress(sprintf("Found additional sex column: %s", col))
    
    col_values <- table(new_metadata[[col]], useNA = "ifany")
    log_progress(sprintf("Values in %s:", col))
    print(col_values)
    
    missing_indices <- which(new_metadata$Sex_Standardized == "Unknown")
    
    
    additional_sex <- case_when(
      new_metadata[[col]] %in% c("M", "Male", "male", "MALE", "m", 1) ~ "Male",
      new_metadata[[col]] %in% c("F", "Female", "female", "FEMALE", "f", 0) ~ "Female",
      TRUE ~ "Unknown"
    )
    
    # Use a single, vectorized ifelse to fill in missing values
    new_metadata$Sex_Standardized <- ifelse(
      new_metadata$Sex_Standardized == "Unknown" & additional_sex != "Unknown",
      additional_sex,
      new_metadata$Sex_Standardized
    )
    
    found_additional_sex <- TRUE
    new_count <- sum(new_metadata$Sex_Standardized %in% c("Male", "Female"))
    additional_contribution <- new_count - existing_sex_count
    log_progress(sprintf("Column %s contributed %d additional cells", col, additional_contribution))
    existing_sex_count <- new_count
  }
}

if(!found_additional_sex) {
  log_progress("No additional sex columns found")
}


#Finalize and Re-assign to Seurat Object
rownames(new_metadata) <- new_metadata$cell_barcode
new_metadata$cell_barcode <- NULL
myeloid_subset@meta.data <- new_metadata

# --- Step 4: Fixing Li APOE_Final column ---

myeloid_subset_fix$APOE_Final <- case_when(
  # If apoe_category has meaningful data (i.e., not 0 or NA), use it
  !is.na(myeloid_subset_fix$apoe_category) & myeloid_subset_fix$apoe_category != "0" ~ myeloid_subset_fix$apoe_category,
  
  # Otherwise, keep the existing value from APOE_Final
  TRUE ~ myeloid_subset_fix$APOE_Final
)

# Set the now-unused 'apoe_category' column to NULL to delete it
myeloid_subset_fix$apoe_category <- NULL
unique(myeloid_subset_fix$sample_id[myeloid_subset_fix$dataset_name == "li"])


### ===== MYELOID SUBCLUSTERING ####
counts <- LayerData(myeloid_subset, assay = "RNA", layer = "counts")

# 2. Extract the metadata 
metadata <- myeloid_subset@meta.data

# 3. Create a new, clean Seurat object from the extracted data
myeloid_subset_fix <- CreateSeuratObject(counts = counts, meta.data = metadata)
rm(myeloid_subset)

# Process for subclustering
DefaultAssay(myeloid_subset_fix) <- "RNA"
myeloid_subset_fix <- NormalizeData(myeloid_subset_fix)
myeloid_subset_fix <- FindVariableFeatures(myeloid_subset_fix, nfeatures = 2000)
myeloid_subset_fix <- ScaleData(myeloid_subset_fix)
myeloid_subset_fix <- RunPCA(myeloid_subset_fix, npcs = 30)


# Now, run Harmony as you did (using your original, correct parameters)
log_progress("Running Harmony integration...")

myeloid_subset_fix <- RunHarmony(
  object = myeloid_subset_fix,
  group.by = "dataset_name",
  reduction.use = "pca",  
  dims = 1:30,           
  plot_convergence = TRUE,
  verbose = TRUE
)

# **IMPORTANT: Now re-run FindNeighbors and RunUMAP using the 'harmony' reduction**
log_progress("Running clustering/UMAP *with* integration...")
myeloid_subset_fix <- FindNeighbors(myeloid_subset_fix, reduction = "harmony", dims = 1:20)
myeloid_subset_fix <- FindClusters(myeloid_subset_fix, resolution = 0.1)
myeloid_subset_fix <- RunUMAP(myeloid_subset_fix, reduction = "harmony", dims = 1:20, reduction.name = "umap.with_harmony")

# **Store the 'with harmony' clusters for clarity**
myeloid_subset_fix$clusters_with_harmony <- myeloid_subset_fix$seurat_clusters


# Keep the original cluster identities
Idents(myeloid_subset_fix) <- "clusters_with_harmony" 

p_micro_macro <- Seurat::DimPlot(
  myeloid_subset_fix, 
  reduction = "umap.with_harmony", 
  group.by = "clusters_with_harmony",
  raster = FALSE , label = TRUE
) +
  labs(title = "Myeloid UMAP") +
  theme_minimal()
p_micro_macro


### ===== MYELOID SUBSET ANALYSIS (UMAP/VIOLIN PLOTS) =====
myeloid_subset_fix$Final_Classification <-NULL
cluster_to_celltype <- c(
  "0" = "Microglia 1",
  "1" = "Microglia 2",
  "2" = "Microglia 3",
  "3" = "BAMs",
  "4" = "Microglia 4",
  "5" = "Peripheral Leukocytes",
  "6" = "Microglia 5",
  "7" = "Microglia 6"
)

# Assign to metadata
myeloid_subset_fix@meta.data$Final_Classification <- cluster_to_celltype[
  as.character(myeloid_subset_fix@meta.data$clusters_with_harmony)
]



##UMAP/Violin Plots

cell_type_order <- c(
  "Peripheral Leukocytes",
  "BAMs",
  "Microglia 1",
  "Microglia 2",
  "Microglia 3",
  "Microglia 4",
  "Microglia 5",
  "Microglia 6"
)

# Reorder the factor levels
myeloid_subset_fix$Final_Classification <- factor(
  myeloid_subset_fix$Final_Classification,
  levels = cell_type_order
)

#Extract UMAP coordinates and create ggplot UMAP
umap_data <- data.frame(
  UMAP_1 = myeloid_subset_fix@reductions$umap.with_harmony@cell.embeddings[, 1],
  UMAP_2 = myeloid_subset_fix@reductions$umap.with_harmony@cell.embeddings[, 2],
  Cell_Type = myeloid_subset_fix$Final_Classification
)

# Calculate centroids for labels
umap_centroids <- umap_data %>%
  group_by(Cell_Type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

# Create custom color palette
cell_colors <- c(
  "Peripheral Leukocytes" = "#FF7F00",
  "BAMs" = "#33A02C",
  "Microglia 1" = "#1F78B4",
  "Microglia 2" = "lavender",
  "Microglia 3" = "#B15928",
  "Microglia 4" = "#A6CEE3",
  "Microglia 5" = "#FB9A99",
  "Microglia 6" = "#FDBF6F"
)

p_umap <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Cell_Type)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  labs(title = "Myeloid Cell Type Annotations",
       x = "UMAP 1",
       y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))+
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

#Extract expression data for marker genes
marker_genes <- c("PTPRC", "CD3E", "CD2", "IL7R",           
                  "ITGAM", "FCN1", "S100A8", "S100A9",
                  "MRC1", "LYVE1", "CD163",                  
                  "P2RY12", "SPP1", "CX3CR1")              

# Check which genes are present
genes_present <- marker_genes[marker_genes %in% rownames(myeloid_subset_fix)]
genes_missing <- marker_genes[!marker_genes %in% rownames(myeloid_subset_fix)]

if(length(genes_missing) > 0) {
  cat("\nWarning: Missing genes:", paste(genes_missing, collapse = ", "), "\n")
}

# Extract expression data
expression_data <- GetAssayData(myeloid_subset_fix, slot = "data")[genes_present, , drop = FALSE]
expression_df <- as.data.frame(t(as.matrix(expression_data)))
expression_df$Cell_Type <- myeloid_subset_fix$Final_Classification
expression_df$Cell_ID <- rownames(expression_df)

# Convert to long format for ggplot
expression_long <- expression_df %>%
  pivot_longer(
    cols = all_of(genes_present),
    names_to = "Gene",
    values_to = "Expression"
  )

# Ensure gene order
expression_long$Gene <- factor(expression_long$Gene, levels = genes_present)

#Violin Plots
violin_plots <- list()

for(gene in genes_present) {
  gene_data <- expression_long %>% filter(Gene == gene)
  
  p <- ggplot(gene_data, aes(x = Cell_Type, y = Expression, fill = Cell_Type)) +
    geom_violin(scale = "width", trim = TRUE) +
    scale_fill_manual(values = cell_colors) +
    labs(title = gene, x = "", y = "Expression") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11,face = "bold"),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      panel.grid.major.x = element_blank()
    )
  
  violin_plots[[gene]] <- p
}

# Combine
combined_violins <- wrap_plots(violin_plots, ncol = 4)

figure_grid <- p_umap / combined_violins + 
  plot_layout(heights = c(1, 2.5))

ggsave(file.path(output_dir, "Cell_Type_Annotation_Violins.png"),
       figure_grid,
       width = 18,
       height = 16,
       dpi = 300,
       bg = "white")
ggsave(file.path(output_dir, "Cell_Type_UMAP.pdf"),
       p_umap,
       width = 24,
       height = 16,
       bg = "white")


##Stacked Violin Plot (genes are rows)
#Load Required Libraries ---
library(tidyr)  

#Define Genes and Colors ---
genes_to_plot <- c("PTPRC", "CD3E", "CD2", "IL7R", "ITGAM", "FCN1", "S100A8", "S100A9",
                   "MRC1", "LYVE1", "CD163", "P2RY12", "SPP1","CX3CR1")


cell_colors <- c(
  "Peripheral Leukocytes" = "#FF7F00",
  "BAMs" = "#33A02C",
  "Microglia 1" = "#1F78B4",
  "Microglia 2" = "lavender",
  "Microglia 3" = "#B15928",
  "Microglia 4" = "#A6CEE3",
  "Microglia 5" = "#FB9A99",
  "Microglia 6" = "#FDBF6F"
)

#Create the Long-Format Data Frame ---

expression_matrix <- GetAssayData(myeloid_subset_fix, 
                                  slot = "data")[genes_to_plot, ]
plot_data <- as.data.frame(t(as.matrix(expression_matrix)))
plot_data$Cell_Type <- myeloid_subset_fix$Final_Classification

expression_long <- pivot_longer(plot_data,
                                cols = -Cell_Type,
                                names_to = "Gene",
                                values_to = "Expression")

#Create the Stacked Plot ---

#Factor the genes to plot in your specified order (top to bottom)
expression_long$Gene <- factor(expression_long$Gene, levels = rev(genes_to_plot))
expression_long$Cell_Type <- factor(expression_long$Cell_Type, levels = names(cell_colors))

p_stacked_violins <- ggplot(expression_long, 
                            aes(x = Cell_Type, y = Expression, fill = Cell_Type)) +
  
  #Add the violins
  geom_violin(scale = "width", trim = TRUE) +
  
  #Add the color scale and legend title
  scale_fill_manual(values = cell_colors, name = "Cluster annotation") +
  facet_grid(Gene ~ ., scales = "free_y", switch = "y") +
  labs(title = "Gene Expression Levels", 
       x = "Cell Type", 
       y = "Expression") +
  theme_minimal() +
  theme(
    #General plot text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(), 
    
    #Y-axis 
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "italic", size = 10),
    
    #Grid lines and spacing
    panel.grid.major.x = element_blank(), 
    panel.spacing.y = unit(0.2, "lines")  
  )

ggsave(file.path(output_dir, "Myeloid_Faceted_Violin_Plot.pdf"),
       p_stacked_violins, 
       width = 10,  # Adjust width as needed
       height = 8,  # Adjust height as needed
       dpi = 300, 
       bg = "white")


### ===== MYELOID SUBSET ANALYSIS(SEX COMPARISON PER CLUSTER) ======

sex_de_results <- list()
sex_summary <- data.frame()
myeloid_subset_fix$Final_Classification <- as.character(myeloid_subset_fix$Final_Classification)

# Get unique clusters
clusters <- unique(myeloid_subset_fix$Final_Classification)
clusters <- clusters[!is.na(clusters)]  


for(cluster in clusters) {
  cat(sprintf("\n=== Analyzing %s: Male vs Female ===\n", cluster))
  cluster_idx <- which(myeloid_subset_fix$Final_Classification == cluster)
  
  if(length(cluster_idx) == 0) {
    cat("Skipping - no cells found for this cluster\n")
    next
  }
  
  cluster_cells <- myeloid_subset_fix[, cluster_idx]
  
  # Remove Unknown sex
  sex_idx <- which(cluster_cells$Sex_Standardized %in% c("Male", "Female"))
  
  if(length(sex_idx) == 0) {
    cat("Skipping - no cells with known sex\n")
    next
  }
  
  cluster_cells <- cluster_cells[, sex_idx]
  
  # Check cell counts
  sex_counts <- table(cluster_cells$Sex_Standardized)
  
  if(length(sex_counts) < 2) {
    cat("Skipping - only one sex present\n")
    next
  }
  cat(sprintf("Female: %d cells, Male: %d cells\n", 
              sex_counts["Female"], sex_counts["Male"]))
  # Skip if too few cells
  if(any(sex_counts < 10)) {
    cat("Skipping - insufficient cells (need >=10 per group)\n")
    next
  }
  
  # Run DE
  Idents(cluster_cells) <- "Sex_Standardized"
  
  tryCatch({
    de_results <- FindMarkers(cluster_cells,
                              ident.1 = "Female",
                              ident.2 = "Male",
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              verbose = FALSE)
    
    de_results$gene <- rownames(de_results)
    de_results$cluster <- cluster
    de_results$comparison <- "Female_vs_Male"
    
    # Filter significant
    de_sig <- de_results %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.25)
    
    cat(sprintf("Total DEGs (p_adj < 0.05, |log2FC| > 0.25): %d\n", nrow(de_sig)))
    cat(sprintf("  Female-biased: %d\n", sum(de_sig$avg_log2FC > 0)))
    cat(sprintf("  Male-biased: %d\n", sum(de_sig$avg_log2FC < 0)))
    
    # Save results
    sex_de_results[[cluster]] <- de_results
    
    write.csv(de_sig,
              file.path(de_output_dir, "Supplemental_Sex", 
                        sprintf("%s_Female_vs_Male_Significant.csv", gsub("/", "_", cluster))),
              row.names = FALSE)
    
  }, error = function(e) {
    cat(sprintf("Error: %s\n", e$message))
  })
}


### ===== MYELOID SUBSET ANALYSIS (APOE E4 COMPARISONS PER CLUSTER) =====

apoe_de_results <- list()
apoe_summary <- data.frame()

for(cluster in clusters) {
  cat(sprintf("\n=== Analyzing %s: E4-expressing vs Non-E4 ===\n", cluster))
  
  #Subsetting
  cluster_idx <- which(myeloid_subset_fix$Final_Classification == cluster)
  
  if(length(cluster_idx) == 0) {
    cat("Skipping - no cells found for this cluster\n")
    next
  }
  cluster_cells <- myeloid_subset_fix[, cluster_idx]
  
  # Remove Unknown APOE
  e4_idx <- which(cluster_cells$E4_Status %in% c("E4_Expressing", "Non_E4"))
  if(length(e4_idx) == 0) {
    cat("Skipping - no cells with known APOE status\n")
    next
  }
  cluster_cells <- cluster_cells[, e4_idx]
  
  # Check cell counts
  e4_counts <- table(cluster_cells$E4_Status)
  if(length(e4_counts) < 2) {
    cat("Skipping - only one E4 status present\n")
    next
  }
  
  cat(sprintf("E4-expressing: %d cells, Non-E4: %d cells\n", 
              e4_counts["E4_Expressing"], e4_counts["Non_E4"]))
  
  # Skip if too few cells or if either group has 0 cells
  if(any(e4_counts == 0) || any(e4_counts < 10)) {
    cat("Skipping - insufficient cells (need >=10 per group)\n")
    next
  }
  
  # Run DE
  Idents(cluster_cells) <- "E4_Status"
  
  tryCatch({
    de_results <- FindMarkers(cluster_cells,
                              ident.1 = "E4_Expressing",
                              ident.2 = "Non_E4",
                              min.pct = 0.1,
                              logfc.threshold = 0.25,
                              test.use = "wilcox",
                              verbose = FALSE)
    
    de_results$gene <- rownames(de_results)
    de_results$cluster <- cluster
    de_results$comparison <- "E4_vs_NonE4"
    
    # Filter significant
    de_sig <- de_results %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
    
    cat(sprintf("Total DEGs (p_adj < 0.05, |log2FC| > 0.25): %d\n", nrow(de_sig)))
    cat(sprintf("  E4-upregulated: %d\n", sum(de_sig$avg_log2FC > 0)))
    cat(sprintf("  E4-downregulated: %d\n", sum(de_sig$avg_log2FC < 0)))
    
    # Save results
    apoe_de_results[[cluster]] <- de_results
    
    write.csv(de_sig,
              file.path(de_output_dir, "Supplemental_APOE", 
                        sprintf("%s_E4_vs_NonE4_Significant.csv", gsub("/", "_", cluster))),
              row.names = FALSE)
    
  }, error = function(e) {
    cat(sprintf("Error: %s\n", e$message))
  })
}


### ===== MYELOID SUBSET ANALYSIS (APOE E4 COMPARISONS SPLIT BY SEX PER CLUSTER) =====
sex_apoe_de_results <- list()
sex_apoe_summary <- data.frame()

for(cluster in clusters) {
  cat(sprintf("\n=== Analyzing %s ===\n", cluster))
  
  #Subsetting
  cluster_idx <- which(myeloid_subset_fix$Final_Classification == cluster)
  
  if(length(cluster_idx) == 0) {
    cat("Skipping - no cells found for this cluster\n")
    next
  }
  
  cluster_cells <- myeloid_subset_fix[, cluster_idx]
  
  # Remove Unknown sex and Unknown APOE
  good_cells_idx <- which(cluster_cells$Sex_Standardized %in% c("Male", "Female") & 
                            cluster_cells$E4_Status %in% c("E4_Expressing", "Non_E4"))
  
  if(length(good_cells_idx) == 0) {
    cat("Skipping - no cells with both known sex and APOE status\n")
    next
  }
  
  cluster_cells <- cluster_cells[, good_cells_idx]
  
  for(sex in c("Female", "Male")) {
    cat(sprintf("\n--- %s: E4 vs Non-E4 in %s ---\n", cluster, sex))
    
    # Subset to this sex
    sex_idx <- which(cluster_cells$Sex_Standardized == sex)
    
    if(length(sex_idx) == 0) {
      cat(sprintf("Skipping - no %s cells\n", sex))
      next
    }
    
    sex_cells <- cluster_cells[, sex_idx]
    
    # Check counts
    e4_counts <- table(sex_cells$E4_Status)
    
    if(length(e4_counts) < 2) {
      cat("Skipping - only one E4 status present\n")
      next
    }
    
    cat(sprintf("E4-expressing: %d cells, Non-E4: %d cells\n", 
                e4_counts["E4_Expressing"], e4_counts["Non_E4"]))
    
    # Skip if 0 cells or too few cells
    if(any(e4_counts == 0) || any(e4_counts < 10)) {
      cat("Skipping - insufficient cells (need >=10 per group)\n")
      next
    }
    
    # Run DE
    Idents(sex_cells) <- "E4_Status"
    
    tryCatch({
      de_results <- FindMarkers(sex_cells,
                                ident.1 = "E4_Expressing",
                                ident.2 = "Non_E4",
                                min.pct = 0.1,
                                logfc.threshold = 0.25,
                                test.use = "wilcox",
                                verbose = FALSE)
      
      de_results$gene <- rownames(de_results)
      de_results$cluster <- cluster
      de_results$sex <- sex
      de_results$comparison <- sprintf("E4_vs_NonE4_%s", sex)
      
      # Filter significant
      de_sig <- de_results %>%
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)
      
      cat(sprintf("Total DEGs: %d\n", nrow(de_sig)))
      cat(sprintf("  E4-upregulated: %d\n", sum(de_sig$avg_log2FC > 0)))
      cat(sprintf("  E4-downregulated: %d\n", sum(de_sig$avg_log2FC < 0)))
      
      # Save results
      comparison_name <- sprintf("%s_%s", cluster, sex)
      sex_apoe_de_results[[comparison_name]] <- de_results
      
      
      write.csv(de_sig,
                file.path(de_output_dir, "Main_Sex_APOE", 
                          sprintf("%s_%s_E4_vs_NonE4_Significant.csv", 
                                  gsub("/", "_", cluster), sex)),
                row.names = FALSE)
      
    }, error = function(e) {
      cat(sprintf("Error: %s\n", e$message))
    })
  }
}


### ===== MYELOID SUBSET ANALYSIS (HLA GENES) =====
library(ComplexHeatmap)
library(circlize) 
output_dir <- "/path/to/output/dir"

#Create the Custom Sex/APOE Metadata ---
group_levels <- c("Female non-E4", "Female E4", "Male non-E4", "Male E4")

if (!"Sex_APOE_Group" %in% colnames(myeloid_subset_fix@meta.data)) {
  myeloid_subset_fix <- myeloid_subset_fix %>%
    AddMetaData(., col.name = "Sex_APOE_Group", 
                metadata = case_when(
                  .$Sex_Standardized == "Female" & .$APOE_Final == "No_E4" ~ "Female non-E4",
                  .$Sex_Standardized == "Female" & .$APOE_Final %in% c("One_E4", "Two_E4") ~ "Female E4",
                  .$Sex_Standardized == "Male" & .$APOE_Final == "No_E4" ~ "Male non-E4",
                  .$Sex_Standardized == "Male" & .$APOE_Final %in% c("One_E4", "Two_E4") ~ "Male E4",
                  TRUE ~ "Other"
                ))
}

#Find HLA-D Genes (Now correctly filtered) ---
available_genes <- rownames(myeloid_subset_fix)
hla_d_genes_all <- grep("^HLA-D", available_genes, value = TRUE)
hla_d_genes_all <- intersect(hla_d_genes_all, available_genes)

myeloid_subset_filtered <- subset(myeloid_subset_fix, 
                                  subset = Sex_APOE_Group %in% group_levels)

hla_expr_sums <- rowSums(GetAssayData(myeloid_subset_filtered, slot = "data")[hla_d_genes_all, ])
hla_d_genes <- names(hla_expr_sums[hla_expr_sums > 0])

if(length(hla_d_genes) == 0) {
  stop("No HLA-D genes *with expression in the 4 target groups* found.")
}
message(sprintf("Found %d total HLA-D genes, %d of which have expression in target groups:", 
                length(hla_d_genes_all), length(hla_d_genes)))
print(hla_d_genes)


#Define Colors and Get Cluster List ---
heatmap_col <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
all_clusters <- sort(unique(myeloid_subset_fix$Final_Classification))
message(sprintf("Found %d clusters to process...", length(all_clusters)))

#Establish a FIXED GENE ORDER (using the filtered object) ---

# Calculate average expression from the filtered object
overall_avg_expr <- AverageExpression(
  myeloid_subset_filtered,
  features = hla_d_genes, # Use the new, filtered gene list
  group.by = "Sex_APOE_Group",
  assay = "RNA"
)$RNA %>% as.matrix()

# Scale this for clustering
overall_scaled_matrix <- t(scale(t(overall_avg_expr)))
overall_scaled_matrix[is.nan(overall_scaled_matrix)] <- 0

# Perform clustering to get a fixed order
gene_hclust <- hclust(dist(overall_scaled_matrix))
fixed_gene_order <- gene_hclust$labels[gene_hclust$order]

#Loop, Calculate Averages, and Save Heatmaps ---
for (cluster_name in all_clusters) {
  message(sprintf("Processing heatmap for cluster: %s", cluster_name))
  
  #Subset the *original* Seurat object
  cluster_subset <- subset(myeloid_subset_fix, 
                           subset = Final_Classification == cluster_name)
  
  #Check for cells in our 4 groups
  cells_in_groups <- sum(cluster_subset$Sex_APOE_Group %in% group_levels)
  if(cells_in_groups == 0) {
    message("Skipping, no cells found in target Sex/APOE groups.")
    next
  }
  
  #Calculate Average Expression
  avg_expr_data <- AverageExpression(
    cluster_subset,
    features = hla_d_genes, #Use the correctly filtered gene list
    group.by = "Sex_APOE_Group",
    assay = "RNA"
  )$RNA
  
  avg_expr <- as.matrix(avg_expr_data)
  
  #Create a full matrix to ensure all 4 columns are present in order
  cluster_expr_matrix <- matrix(
    NA, 
    nrow = length(hla_d_genes),
    ncol = length(group_levels),
    dimnames = list(hla_d_genes, group_levels) 
  )
  
  #Find common genes AND groups
  present_groups <- intersect(group_levels, colnames(avg_expr))
  present_genes <- intersect(hla_d_genes, rownames(avg_expr)) 
  
  #Fill the template matrix *by name*
  if (length(present_genes) > 0 && length(present_groups) > 0) {
    cluster_expr_matrix[present_genes, present_groups] <- avg_expr[present_genes, present_groups, drop = FALSE]
  }
  
  #Scale the data row-wise
  cluster_scaled_matrix <- t(scale(t(cluster_expr_matrix)))
  cluster_scaled_matrix[is.nan(cluster_scaled_matrix)] <- 0
  
  #Create the Heatmap
  ht <- Heatmap(
    cluster_scaled_matrix,
    name = "Scaled Avg. Expr.",
    col = heatmap_col,
    na_col = "grey90", 
    
    # Columns
    cluster_columns = FALSE,     
    column_order = group_levels, 
    column_names_gp = gpar(fontsize = 10),
    column_title = cluster_name,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    # Rows
    cluster_rows = FALSE, 
    row_order = fixed_gene_order, # Use the fixed order
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8)
  )
  
  #Save to PDF
  safe_cluster_name <- gsub("[^a-zA-Z0-9_]", "_", cluster_name)
  pdf_filename <- file.path(output_dir, 
                            sprintf("Heatmap_HLA-D_AvgExpr_%s.pdf", safe_cluster_name))
  
  pdf(pdf_filename, width = 7, height = 8)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
}