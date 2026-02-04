#!/usr/bin/env Rscript

###Run this script in the same subdirectory where the GDC data was downloaded
###Run separate ancestry groups in separate subdirectories

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

cat("MAF Somatic Mutation Frequency Analysis\n\n")

###(1) LOCATE MAF FILES###

# Find all MAF files one subdirectory down
find_maf_files <- function() {
  cat("(1) Searching for MAF files in subdirectories...\n")
  
  # Find all .maf files in immediate subdirectories
  maf_files <- list.files(
    path = ".",
    pattern = "\\.maf$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Filter to only include files that are exactly one subdirectory down
  maf_files <- maf_files[grepl("^\\./[^/]+/[^/]+\\.maf$", maf_files)]
  
  if (length(maf_files) == 0) {
    cat("No MAF files found in subdirectories.\n")
    
    stop("No valid MAF files found.")
  }
  
  cat(sprintf("Found %d MAF files:\n", length(maf_files)))
  for (i in seq_along(maf_files)) {
    cat(sprintf("  %d. %s\n", i, maf_files[i]))
  }
  cat("\n")
  
  return(maf_files)
}

###(2) READ MAF FILES###

cat("(2) Reading MAF files...\n")

# Function to read MAF file

read_maf_safe <- function(file_path) {
  cat(sprintf("Reading: %s\n", basename(file_path)))
  
  tryCatch({
    # First, check if file exists and is readable
    if (!file.exists(file_path)) {
      cat(sprintf("File does not exist: %s\n", file_path))
      return(NULL)
    }
    
    if (file.size(file_path) == 0) {
      cat(sprintf("File is empty: %s\n", basename(file_path)))
      return(NULL)
    }
    
    # Read with data.table
    maf_data <- fread(file_path, 
                      sep = "\t", 
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      showProgress = FALSE)
    
    if (nrow(maf_data) == 0) {
      cat(sprintf("File contains no data rows: %s\n", basename(file_path)))
      return(NULL)
    }
    
    # Check for required columns
    required_cols <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome")
    missing_cols <- setdiff(required_cols, colnames(maf_data))
    
    if (length(missing_cols) > 0) {
      cat(sprintf("Missing required columns in %s: %s\n", 
                  basename(file_path), 
                  paste(missing_cols, collapse = ", ")))
      cat(sprintf("Available columns: %s\n", 
                  paste(head(colnames(maf_data), 10), collapse = ", ")))
      return(NULL)
    }
    
    # Filter for somatic mutations (exclude germline if present)
    if ("Mutation_Status" %in% colnames(maf_data)) {
      before_count <- nrow(maf_data)
      maf_data <- maf_data[Mutation_Status == "Somatic" | is.na(Mutation_Status)]
      cat(sprintf("Filtered to somatic mutations: %d -> %d\n", before_count, nrow(maf_data)))
    }
    
    # Remove silent mutations, if present
    if ("Variant_Classification" %in% colnames(maf_data)) {
      before_count <- nrow(maf_data)
      maf_data <- maf_data[Variant_Classification != "Silent"]
      cat(sprintf("Removed silent mutations: %d -> %d\n", before_count, nrow(maf_data)))
    }
    
    if (nrow(maf_data) == 0) {
      cat(sprintf("No mutations remaining after filtering in %s\n", basename(file_path)))
      return(NULL)
    }
    
    cat(sprintf("Loaded %d mutations from %d unique samples\n", 
                nrow(maf_data), 
                length(unique(maf_data$Tumor_Sample_Barcode))))
    return(maf_data)
    
  }, error = function(e) {
    cat(sprintf("Error reading %s: %s\n", basename(file_path), e$message))
    return(NULL)
  })
}

###(3) PROCESS MAF FILES###

# Function to process all MAF files
process_all_mafs <- function(maf_files) {
  cat("(3) Processing MAF files...\n")
  
  all_mutations <- list()
  total_samples <- character(0)
  
  for (i in seq_along(maf_files)) {
    maf_data <- read_maf_safe(maf_files[i])
    
    if (!is.null(maf_data) && nrow(maf_data) > 0) {
      # Store mutations with source file info
      maf_data$source_file <- basename(maf_files[i])
      all_mutations[[i]] <- maf_data
      
      # Collect unique sample IDs
      total_samples <- c(total_samples, unique(maf_data$Tumor_Sample_Barcode))
    } else {
      cat(sprintf("Skipping file with no valid data: %s\n", basename(maf_files[i])))
    }
  }
  
  if (length(all_mutations) == 0) {
    stop("No valid MAF data could be loaded from any file!")
  }
  
  # Combine all mutations
  combined_mutations <- rbindlist(all_mutations, fill = TRUE)
  total_unique_samples <- unique(total_samples)
  
  cat(sprintf("Combined data: %d mutations across %d unique samples from %d files\n", 
              nrow(combined_mutations), 
              length(total_unique_samples),
              length(all_mutations)))
  cat("\n")
  
  return(list(
    mutations = combined_mutations,
    total_samples = length(total_unique_samples)
  ))
}

###(4) CALCULATE MUTATION FREQUENCIES###

# Function to calculate mutation frequencies
calculate_mutation_frequencies <- function(mutations_data, total_samples) {
  cat("(4) Calculating mutation frequencies per gene...\n")
  
  mutations <- mutations_data$mutations
  
  # Calculate per-gene statistics
  gene_stats <- mutations %>%
    group_by(Hugo_Symbol, Chromosome) %>%
    summarise(
      mutated_samples = n_distinct(Tumor_Sample_Barcode),
      total_mutations = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      mutation_percentage = round((mutated_samples / total_samples) * 100, 2)
    ) %>%
    arrange(desc(mutated_samples)) %>%
    slice_head(n = 25)
  
  # Clean chromosome names
  gene_stats$Chromosome <- gsub("^chr", "", gene_stats$Chromosome)
  gene_stats$Chromosome <- gsub("^23$", "X", gene_stats$Chromosome)
  gene_stats$Chromosome <- gsub("^24$", "Y", gene_stats$Chromosome)
  
  cat(sprintf("Calculated frequencies for %d genes\n", nrow(gene_stats)))
  cat("Top 5 most frequently mutated genes in this cohort:\n")
  
  for (i in 1:min(5, nrow(gene_stats))) {
    cat(sprintf("  %d. %s (chr%s): %d/%d samples (%.2f%%)\n",
                i,
                gene_stats$Hugo_Symbol[i],
                gene_stats$Chromosome[i],
                gene_stats$mutated_samples[i],
                total_samples,
                gene_stats$mutation_percentage[i]))
  }
  cat("\n")
  
  return(gene_stats)
}

###(5) FORMAT RESULTS AS TSV###

# Function to write formatted TSV results
write_results <- function(gene_stats, total_samples) {
  output_file <- "top_25_mutated_genes.tsv"
  
  cat("(5) Writing results to", output_file, "...\n")
  
  # Prepare final output with proper column names
  final_results <- gene_stats %>%
    select(
      Gene_Symbol = Hugo_Symbol,
      Chromosome = Chromosome,
      Samples_With_Mutations = mutated_samples,
      Total_Mutations = total_mutations,
      Mutation_Percentage = mutation_percentage
    )
  
  # Create file connection
  con <- file(output_file, "w")
  
  # Write header comments
  writeLines(c(
    "# MAF Somatic Mutation Frequency Analysis",
    paste0("# Generated: ", Sys.time()),
    paste0("# Total samples analyzed: ", total_samples),
    "# Top 25 most frequently mutated genes",
    "#"
  ), con)
  
  # Write column headers with tab separation
  header_line <- paste(colnames(final_results), collapse = "\t")
  writeLines(header_line, con)
  
  # Write data rows with tabs
  for (i in 1:nrow(final_results)) {
    data_line <- paste(
      final_results$Gene_Symbol[i],
      final_results$Chromosome[i], 
      final_results$Samples_With_Mutations[i],
      final_results$Total_Mutations[i],
      final_results$Mutation_Percentage[i],
      sep = "\t"
    )
    writeLines(data_line, con)
  }
  
  # Close the file connection
  close(con)
  
  cat("Results written successfully!\n")
  cat(sprintf("Analysis Summary:\n"))
  cat(sprintf("   - Total samples: %d\n", total_samples))
  cat(sprintf("   - Genes analyzed: %d\n", nrow(final_results)))
  cat(sprintf("   - Output file: %s\n", output_file))
  
  return(output_file)
}

# Main function
main <- function() {
  start_time <- Sys.time()
  
  tryCatch({
    cat("Starting MAF analysis...\n\n")
    
    # Step 1: Find & read MAF files
    maf_files <- find_maf_files()
    
    # Step 2: Process all MAF files
    mutations_data <- process_all_mafs(maf_files)
    
    # Step 3: Calculate mutation frequencies
    gene_stats <- calculate_mutation_frequencies(mutations_data, mutations_data$total_samples)
    
    # Step 4: Write results
    output_file <- write_results(gene_stats, mutations_data$total_samples)
    
    end_time <- Sys.time()
    cat(sprintf("\nAnalysis completed successfully in %.2f seconds.", 
                as.numeric(difftime(end_time, start_time, units = "secs"))))
    
  }, error = function(e) {
    cat(sprintf("Fatal error: %s\n", e$message))
    cat("Debugging information:\n")
    cat(sprintf("  - Current directory: %s\n", getwd()))
    cat(sprintf("  - R version: %s\n", R.version.string))
    cat("  - Working directory contents:\n")
    files <- list.files(".", recursive = FALSE)
    for (file in head(files, 10)) {
      cat(sprintf("    %s\n", file))
    }
    if (length(files) > 10) {
      cat(sprintf("    ... and %d more files\n", length(files) - 10))
    }
    
    # Return error status if running as script
    if (!interactive()) {
      quit(status = 1)
    }
  })
}

cat("Running in", ifelse(interactive(), "interactive", "script"), "mode\n\n")

# Run the analysis immediately
main()
