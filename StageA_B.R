# ===================================================
# 1. Required Libraries
# ===================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

libs <- c("GEOquery", "affy", "oligo", "Biobase", "AnnotationDbi", 
          "hgu133plus2.db", "hgu133a.db", "hgu133plus2cdf", "hgu133acdf", 
          "hugene10sttranscriptcluster.db", "openxlsx", "R.utils")

for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) BiocManager::install(lib)
  library(lib, character.only = TRUE)
}

# ===================================================
# 2. Configuration & Directory Setup
# ===================================================
base_dir   <- "C:/UBILab_Research"
raw_dir    <- file.path(base_dir, "Data/GEO_Raw")
cel_dir    <- file.path(base_dir, "Data/CEL_Extracted")
norm_dir   <- file.path(base_dir, "Data/Normalized")
gene_dir   <- file.path(base_dir, "Data/Gene_Mapped")
merged_dir <- file.path(base_dir, "Data/Merged")

dirs <- c(raw_dir, cel_dir, norm_dir, gene_dir, merged_dir)
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

gse_full_list <- c("GSE55235", "GSE55584", "GSE77298", "GSE48780", 
                   "GSE206848", "GSE55457", "GSE93720", "GSE97779", 
                   "GSE15061", "GSE62872")

st_array_sets <- c("GSE62872") 

# ===================================================
# 3. Robust Platform Handling Functions
# ===================================================

# Standardizes platform names from CEL headers to match Bioconductor packages
ensure_cdf_loaded <- function(chip_type) {
  # Remove all special characters and lowercase (e.g., "HG-U133_Plus_2" -> "hgu133plus2")
  clean_type <- gsub("[^[:alnum:]]", "", tolower(chip_type))
  message("Attempting to match platform: ", clean_type)
  
  if (grepl("hgu133plus2", clean_type) || grepl("hgu133plus", clean_type)) {
    library(hgu133plus2cdf)
    return("hgu133plus2")
  }
  if (grepl("hgu133a", clean_type)) {
    library(hgu133acdf)
    return("hgu133a")
  }
  # Fallback for HuGene ST arrays processed via oligo
  if (grepl("hugene10st", clean_type)) {
    return("hugene10sttranscriptcluster")
  }
  
  stop("Platform not recognized: ", chip_type)
}

get_db_from_platform <- function(platform_name) {
  p_clean <- gsub("[^[:alnum:]]", "", tolower(platform_name))
  if (grepl("hugene10st", p_clean)) return(hugene10sttranscriptcluster.db)
  if (grepl("hgu133plus2", p_clean)) return(hgu133plus2.db)
  if (grepl("hgu133a", p_clean)) return(hgu133a.db)
  stop("No database found for: ", platform_name)
}

write_table_xlsx <- function(df, path) {
  write.xlsx(df, path, rowNames = FALSE, overwrite = TRUE)
}

# ===================================================
# 4. Processing Pipeline
# ===================================================

# List to keep track of successfully normalized platforms
success_platforms <- list()

for (gse_id in gse_full_list) {
  cel_files <- list.files(file.path(cel_dir, gse_id), pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)
  if (length(cel_files) == 0) next
  
  out_norm <- file.path(norm_dir, paste0(gse_id, "_expression_data.xlsx"))
  message("\n--- Processing: ", gse_id, " ---")
  
  tryCatch({
    # A. Normalization Stage
    if (gse_id %in% st_array_sets) {
      raw_data  <- oligo::read.celfiles(cel_files)
      norm_data <- oligo::rma(raw_data, target = "core")
      platform_name <- "hugene10sttranscriptcluster"
    } else {
      header <- affyio::read.celfile.header(cel_files[1])
      platform_name <- ensure_cdf_loaded(header[["cdfName"]])
      raw_data  <- affy::ReadAffy(filenames = cel_files)
      norm_data <- affy::rma(raw_data)
    }
    
    expr_mat <- Biobase::exprs(norm_data)
    expr_df  <- data.frame(probe_id = rownames(expr_mat), expr_mat, check.names = FALSE)
    write_table_xlsx(expr_df, out_norm)
    
    # B. Mapping Stage (Immediate mapping to ensure Gene_Mapped folder fills up)
    message("Mapping Probes to Genes for: ", gse_id)
    db_obj <- get_db_from_platform(platform_name)
    gene_symbols <- mapIds(db_obj, keys = as.character(expr_df$probe_id), 
                           column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
    
    expr_df$gene <- unname(gene_symbols)
    gene_mapped <- expr_df[!is.na(expr_df$gene) & expr_df$gene != "", ]
    gene_mapped <- gene_mapped[!duplicated(gene_mapped$gene), ]
    
    # Clean up column names and reorder
    final_df <- gene_mapped[, c("gene", setdiff(colnames(gene_mapped), c("probe_id", "gene")))]
    colnames(final_df) <- gsub("^(GSM\\d+).*", "\\1", colnames(final_df))
    
    write_table_xlsx(final_df, file.path(gene_dir, paste0(gse_id, "_gene_unique.xlsx")))
    success_platforms[[gse_id]] <- platform_name
    
    message("SUCCESS: ", gse_id, " is fully processed and mapped.")
    
  }, error = function(e) message("CRITICAL ERROR for ", gse_id, ": ", e$message))
}

# ===================================================
# 5. Final Compendium Merging
# ===================================================
gene_files <- list.files(gene_dir, pattern = "_gene_unique\\.xlsx$", full.names = TRUE)

if (length(gene_files) > 1) {
  message("\n--- Creating Final Compendium ---")
  all_data <- lapply(gene_files, function(f) {
    df <- read.xlsx(f)
    rownames(df) <- df$gene
    df$gene <- NULL
    return(df)
  })
  
  common_genes <- Reduce(intersect, lapply(all_data, rownames))
  
  if (length(common_genes) > 0) {
    merged_mat <- do.call(cbind, lapply(all_data, function(df) df[common_genes, , drop=FALSE]))
    merged_df  <- data.frame(gene = common_genes, merged_mat, check.names = FALSE)
    write_table_xlsx(merged_df, file.path(merged_dir, "ra_merged_source_compendium.xlsx"))
    message("ALL PROCESSES FINISHED! Check 'Merged' folder for output.")
  }
}