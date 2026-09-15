# Diet vs Urine Data Correlation Analysis in R
# Author: Harsha Gouda
# Description: Pearson correlations between food-diary intake and urine food composition scores
#              in the Crohn's food re-introduction cohort (Figure 5; MassIVE MSV000096125)
# Place the input files (see README) in the working directory, or set the path below.
# setwd("path/to/input/files")
# Load required libraries
library(tidyverse)
library(corrr)
library(Hmisc)
library(writexl)
library(openxlsx)

# Function to load and clean data
load_and_prepare_data <- function(diet_file, urine_file) {
  cat("Loading datasets...\n")
  
  # Read diet data
  diet_data <- read.csv(diet_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Diet data loaded:", nrow(diet_data), "rows,", ncol(diet_data), "columns\n")
  
  # Read urine data
  urine_data <- read.csv(urine_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Urine data loaded:", nrow(urine_data), "rows,", ncol(urine_data), "columns\n")
  
  # Merge datasets on filename
  merged_data <- merge(diet_data, urine_data, by = "filename", all = FALSE)
  cat("Merged data:", nrow(merged_data), "matched samples\n")
  
  # Identify numeric columns
  diet_numeric_cols <- setdiff(names(diet_data), c("Patient_ID", "filename", "Group", "TimePoint"))
  diet_numeric_cols <- diet_numeric_cols[sapply(diet_data[diet_numeric_cols], is.numeric)]
  
  urine_numeric_cols <- setdiff(names(urine_data), "filename")
  urine_numeric_cols <- urine_numeric_cols[sapply(urine_data[urine_numeric_cols], is.numeric)]
  
  cat("Diet numeric variables:", length(diet_numeric_cols), "\n")
  cat("Urine numeric variables:", length(urine_numeric_cols), "\n")
  
  return(list(
    merged_data = merged_data,
    diet_numeric_cols = diet_numeric_cols,
    urine_numeric_cols = urine_numeric_cols
  ))
}

# Function for selective correlation analysis
selective_correlation_analysis <- function(merged_data, diet_cols, urine_cols, 
                                           significance_level = 0.05, method = "pearson") {
  cat("Running selective correlation analysis...\n")
  cat("Diet columns:", length(diet_cols), "\n")
  cat("Urine columns:", length(urine_cols), "\n")
  cat("Total correlations to compute:", length(diet_cols) * length(urine_cols), "\n")
  
  results <- data.frame()
  
  for (diet_col in diet_cols) {
    for (urine_col in urine_cols) {
      # Extract values and remove missing data
      diet_values <- merged_data[[diet_col]]
      urine_values <- merged_data[[urine_col]]
      
      # Create complete cases
      complete_cases <- complete.cases(diet_values, urine_values)
      diet_clean <- diet_values[complete_cases]
      urine_clean <- urine_values[complete_cases]
      
      if (length(diet_clean) >= 3 && sd(diet_clean) > 0 && sd(urine_clean) > 0) {
        # Calculate correlation
        cor_test <- cor.test(diet_clean, urine_clean, method = method)
        
        # Store results
        results <- rbind(results, data.frame(
          diet_column = diet_col,
          urine_column = urine_col,
          correlation = cor_test$estimate,
          p_value = cor_test$p.value,
          confidence_low = cor_test$conf.int[1],
          confidence_high = cor_test$conf.int[2],
          sample_size = length(diet_clean),
          significant = cor_test$p.value < significance_level,
          method = method,
          stringsAsFactors = FALSE
        ))
      } else {
        # Store NA results for insufficient data
        results <- rbind(results, data.frame(
          diet_column = diet_col,
          urine_column = urine_col,
          correlation = NA,
          p_value = NA,
          confidence_low = NA,
          confidence_high = NA,
          sample_size = sum(complete_cases),
          significant = FALSE,
          method = method,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Sort by absolute correlation value
  results <- results[order(-abs(results$correlation), na.last = TRUE), ]
  row.names(results) <- NULL
  
  return(results)
}

# Function for comprehensive correlation analysis
all_correlation_analysis <- function(merged_data, diet_numeric_cols, urine_numeric_cols,
                                     significance_level = 0.05, min_correlation = 0.1, 
                                     method = "pearson") {
  cat("Running comprehensive correlation analysis...\n")
  cat("Total correlations to compute:", length(diet_numeric_cols) * length(urine_numeric_cols), "\n")
  
  results <- data.frame()
  total_correlations <- length(diet_numeric_cols) * length(urine_numeric_cols)
  completed <- 0
  
  for (diet_col in diet_numeric_cols) {
    for (urine_col in urine_numeric_cols) {
      completed <- completed + 1
      
      # Progress update
      if (completed %% 1000 == 0) {
        cat("Progress:", completed, "/", total_correlations, 
            "(", round(100 * completed / total_correlations, 1), "%)\n")
      }
      
      # Extract values and remove missing data
      diet_values <- merged_data[[diet_col]]
      urine_values <- merged_data[[urine_col]]
      
      # Create complete cases
      complete_cases <- complete.cases(diet_values, urine_values)
      diet_clean <- diet_values[complete_cases]
      urine_clean <- urine_values[complete_cases]
      
      if (length(diet_clean) >= 3 && sd(diet_clean) > 0 && sd(urine_clean) > 0) {
        # Calculate correlation
        cor_test <- cor.test(diet_clean, urine_clean, method = method)
        
        # Only store if meets minimum correlation threshold
        if (!is.na(cor_test$estimate) && abs(cor_test$estimate) >= min_correlation) {
          results <- rbind(results, data.frame(
            diet_column = diet_col,
            urine_column = urine_col,
            correlation = cor_test$estimate,
            p_value = cor_test$p.value,
            confidence_low = cor_test$conf.int[1],
            confidence_high = cor_test$conf.int[2],
            sample_size = length(diet_clean),
            significant = cor_test$p.value < significance_level,
            method = method,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Sort by absolute correlation value
  results <- results[order(-abs(results$correlation), na.last = TRUE), ]
  row.names(results) <- NULL
  
  cat("Completed! Found", nrow(results), "correlations >=", min_correlation, "\n")
  return(results)
}

# Function to generate summary statistics
generate_summary <- function(results, significance_level = 0.05) {
  summary_stats <- list(
    total_correlations = nrow(results),
    significant_correlations = sum(results$significant, na.rm = TRUE),
    positive_correlations = sum(results$correlation > 0, na.rm = TRUE),
    negative_correlations = sum(results$correlation < 0, na.rm = TRUE),
    strong_correlations = sum(abs(results$correlation) >= 0.5, na.rm = TRUE),
    moderate_correlations = sum(abs(results$correlation) >= 0.3 & abs(results$correlation) < 0.5, na.rm = TRUE),
    weak_correlations = sum(abs(results$correlation) >= 0.1 & abs(results$correlation) < 0.3, na.rm = TRUE),
    mean_correlation = mean(abs(results$correlation), na.rm = TRUE),
    median_correlation = median(abs(results$correlation), na.rm = TRUE)
  )
  
  return(summary_stats)
}

# Function to create visualization
create_correlation_plots <- function(results, top_n = 20) {
  library(ggplot2)
  library(gridExtra)
  
  # Top correlations plot
  top_results <- head(results[!is.na(results$correlation), ], top_n)
  
  p1 <- ggplot(top_results, aes(x = reorder(paste(diet_column, "-", urine_column), correlation), 
                                y = correlation, fill = significant)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red3")) +
    labs(title = paste("Top", top_n, "Correlations"),
         x = "Variable Pairs", y = "Pearson Correlation (r)",
         fill = "Significant") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # Distribution of correlations
  p2 <- ggplot(results[!is.na(results$correlation), ], aes(x = correlation)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Distribution of Correlations",
         x = "Pearson Correlation (r)", y = "Frequency") +
    theme_minimal()
  
  # P-value vs correlation
  p3 <- ggplot(results[!is.na(results$correlation) & !is.na(results$p_value), ], 
               aes(x = abs(correlation), y = -log10(p_value), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red3")) +
    labs(title = "Correlation vs Significance",
         x = "|Pearson Correlation (r)|", y = "-log10(p-value)",
         color = "Significant") +
    theme_minimal()
  
  return(list(top_correlations = p1, distribution = p2, volcano = p3))
}

# Function to export results
export_results <- function(results, summary_stats, filename_prefix = "correlation_analysis") {
  # Create timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Export to Excel with multiple sheets
  excel_filename <- paste0(filename_prefix, "_", timestamp, ".xlsx")
  
  wb <- createWorkbook()
  
  # Results sheet
  addWorksheet(wb, "Correlation_Results")
  writeData(wb, "Correlation_Results", results)
  
  # Summary sheet
  addWorksheet(wb, "Summary_Statistics")
  summary_df <- data.frame(
    Metric = names(summary_stats),
    Value = unlist(summary_stats)
  )
  writeData(wb, "Summary_Statistics", summary_df)
  
  # Save Excel file
  saveWorkbook(wb, excel_filename, overwrite = TRUE)
  cat("Results exported to:", excel_filename, "\n")
  
  # Also export CSV
  csv_filename <- paste0(filename_prefix, "_", timestamp, ".csv")
  write.csv(results, csv_filename, row.names = FALSE)
  cat("Results exported to:", csv_filename, "\n")
  
  return(list(excel = excel_filename, csv = csv_filename))
}

# Main analysis function
run_correlation_analysis <- function(diet_file = "IPENS_diet_metadata.csv",
                                     urine_file = "food_metadata_level3.csv",
                                     analysis_type = "both", # "selective", "all", or "both"
                                     selected_diet_cols = NULL,
                                     selected_urine_cols = NULL,
                                     significance_level = 0.05,
                                     min_correlation = 0.1,
                                     method = "pearson",
                                     export_results_flag = TRUE,
                                     create_plots = TRUE) {
  
  # Load and prepare data
  data_list <- load_and_prepare_data(diet_file, urine_file)
  merged_data <- data_list$merged_data
  diet_numeric_cols <- data_list$diet_numeric_cols
  urine_numeric_cols <- data_list$urine_numeric_cols
  
  results_list <- list()
  
  # Selective analysis
  if (analysis_type %in% c("selective", "both")) {
    if (is.null(selected_diet_cols)) {
      # Use first 10 diet columns as example
      selected_diet_cols <- head(diet_numeric_cols, 16)
      cat("No diet columns specified. Using first 10 as example.\n")
    }
    
    if (is.null(selected_urine_cols)) {
      # Use all urine columns
      selected_urine_cols <- urine_numeric_cols
      cat("No urine columns specified. Using all urine columns.\n")
    }
    
    selective_results <- selective_correlation_analysis(
      merged_data, selected_diet_cols, selected_urine_cols, significance_level, method
    )
    
    selective_summary <- generate_summary(selective_results, significance_level)
    
    cat("\n=== SELECTIVE ANALYSIS SUMMARY ===\n")
    cat("Total correlations:", selective_summary$total_correlations, "\n")
    cat("Significant correlations:", selective_summary$significant_correlations, "\n")
    cat("Strong correlations (|r| >= 0.5):", selective_summary$strong_correlations, "\n")
    cat("Mean |correlation|:", round(selective_summary$mean_correlation, 3), "\n")
    
    results_list$selective <- list(results = selective_results, summary = selective_summary)
    
    if (export_results_flag) {
      export_results(selective_results, selective_summary, "selective_correlation_analysis")
    }
  }
  
  # Comprehensive analysis
  if (analysis_type %in% c("all", "both")) {
    cat("\n=== RUNNING COMPREHENSIVE ANALYSIS ===\n")
    all_results <- all_correlation_analysis(
      merged_data, diet_numeric_cols, urine_numeric_cols, significance_level, min_correlation, method
    )
    
    all_summary <- generate_summary(all_results, significance_level)
    
    cat("\n=== COMPREHENSIVE ANALYSIS SUMMARY ===\n")
    cat("Total correlations (|r| >=", min_correlation, "):", all_summary$total_correlations, "\n")
    cat("Significant correlations:", all_summary$significant_correlations, "\n")
    cat("Strong correlations (|r| >= 0.5):", all_summary$strong_correlations, "\n")
    cat("Mean |correlation|:", round(all_summary$mean_correlation, 3), "\n")
    
    results_list$comprehensive <- list(results = all_results, summary = all_summary)
    
    if (export_results_flag) {
      export_results(all_results, all_summary, "comprehensive_correlation_analysis")
    }
  }
  
  # Create plots
  if (create_plots && length(results_list) > 0) {
    for (analysis_name in names(results_list)) {
      results <- results_list[[analysis_name]]$results
      if (nrow(results) > 0) {
        plots <- create_correlation_plots(results)
        
        # Save plots
        pdf(paste0(analysis_name, "_correlation_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"))
        print(plots$top_correlations)
        print(plots$distribution)
        print(plots$volcano)
        dev.off()
        
        cat("Plots saved for", analysis_name, "analysis\n")
      }
    }
  }
  
  return(results_list)
}

# ===============================================================================
# USAGE EXAMPLES
# ===============================================================================

# Example 1: Run selective analysis with specific columns
example_selective_analysis <- function() {
  # Define specific columns of interest
  selected_diet <- c("L1_Fruits_g", "L1_Vegetables_g", "L1_Meat & meat products_g")
  
  selected_urine <- c("Complex_Meat", "Fruit", "Vegetable", "Egg")
  
  results <- run_correlation_analysis(
    analysis_type = "selective",
    selected_diet_cols = selected_diet,
    selected_urine_cols = selected_urine,
    significance_level = 0.05
  )
  
  return(results)
}

# Example 2: Run comprehensive analysis
example_comprehensive_analysis <- function() {
  results <- run_correlation_analysis(
    analysis_type = "all",
    min_correlation = 0.2,  # Only show correlations >= 0.2
    significance_level = 0.01  # Stricter significance level
  )
  
  return(results)
}

# Example 3: Run both analyses
example_both_analyses <- function() {
  results <- run_correlation_analysis(
    analysis_type = "both",
    min_correlation = 0.15,
    significance_level = 0.05
  )
  
  return(results)
}

# ===============================================================================
# RUN THE ANALYSIS
# ===============================================================================

# Uncomment the line below to run the analysis
results <- run_correlation_analysis()

cat("
=== CORRELATION ANALYSIS SCRIPT LOADED ===

To run the analysis, use one of these commands:

1. Basic run with default settings:
   results <- run_correlation_analysis()

2. Selective analysis with specific columns:
   results <- example_selective_analysis()

3. Comprehensive analysis only:
   results <- example_comprehensive_analysis()

4. Both analyses:
   results <- example_both_analyses()

5. Custom analysis:
   results <- run_correlation_analysis(
     analysis_type = 'both',
     significance_level = 0.01,
     min_correlation = 0.2
   )

The script will automatically:
- Load and merge your datasets
- Calculate Pearson correlations
- Generate summary statistics
- Export results to Excel and CSV
- Create visualization plots
- Save everything with timestamps

Files expected:
- IPENS_diet_metadata.csv
- food_readout_urine_UKfoodgroup.csv
")



##############################################################################################################
# Simple Hierarchical Clustered Heatmap for Diet-Urine Correlations
# Uses only pheatmap - no ComplexHeatmap dependency
# Excludes "Days" variable and keeps original urine names

# Load required libraries (only what's available)
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)

# Function to load and prepare correlation data
load_and_prepare_heatmap_data <- function(excel_file = "selective_correlation_analysis_20260109_120619.xlsx") {
  
  # Read the correlation results
  correlation_data <- read_excel(excel_file, sheet = "Correlation_Results")
  
  # Remove "Days" variable and filter for significant correlations
  filtered_data <- correlation_data %>%
    filter(diet_column != "Days",  # Remove Days variable
           significant == TRUE) %>%   # Only significant correlations
    dplyr::select(diet_column, urine_column, correlation, p_value, significant)
  
  cat("Loaded correlation data:\n")
  cat("Total correlations (excluding Days):", nrow(filtered_data), "\n")
  cat("Correlation range:", round(range(filtered_data$correlation, na.rm = TRUE), 3), "\n")
  
  # Show what we have
  cat("\nDiet variables (", length(unique(filtered_data$diet_column)), "):\n")
  cat(paste(sort(unique(filtered_data$diet_column)), collapse = ", "), "\n")
  
  cat("\nUrine variables (", length(unique(filtered_data$urine_column)), "):\n")
  cat(paste(sort(unique(filtered_data$urine_column)), collapse = ", "), "\n")
  
  return(filtered_data)
}

# Function to create correlation matrix for heatmap
create_correlation_matrix <- function(correlation_data) {
  
  # Get unique variables
  diet_vars <- sort(unique(correlation_data$diet_column))
  urine_vars <- sort(unique(correlation_data$urine_column))
  
  # Create matrix
  corr_matrix <- matrix(NA, 
                        nrow = length(diet_vars), 
                        ncol = length(urine_vars),
                        dimnames = list(diet_vars, urine_vars))
  
  # Fill matrix with correlation values
  for(i in 1:nrow(correlation_data)) {
    diet_var <- correlation_data$diet_column[i]
    urine_var <- correlation_data$urine_column[i]
    corr_val <- correlation_data$correlation[i]
    
    corr_matrix[diet_var, urine_var] <- corr_val
  }
  
  return(corr_matrix)
}

# Function to create significance matrix for annotation
create_significance_matrix <- function(correlation_data) {
  
  # Get unique variables
  diet_vars <- sort(unique(correlation_data$diet_column))
  urine_vars <- sort(unique(correlation_data$urine_column))
  
  # Create matrix
  sig_matrix <- matrix(NA, 
                       nrow = length(diet_vars), 
                       ncol = length(urine_vars),
                       dimnames = list(diet_vars, urine_vars))
  
  # Fill matrix with significance levels
  for(i in 1:nrow(correlation_data)) {
    diet_var <- correlation_data$diet_column[i]
    urine_var <- correlation_data$urine_column[i]
    p_val <- correlation_data$p_value[i]
    
    # Convert p-values to significance levels
    if (!is.na(p_val)) {
      if (p_val < 0.001) {
        sig_level <- "***"
      } else if (p_val < 0.01) {
        sig_level <- "**"
      } else if (p_val < 0.05) {
        sig_level <- "*"
      } else {
        sig_level <- ""
      }
    } else {
      sig_level <- ""
    }
    
    sig_matrix[diet_var, urine_var] <- sig_level
  }
  
  return(sig_matrix)
}

# Function to create annotations for diet variables
create_diet_annotations <- function(diet_vars) {
  
  diet_groups <- case_when(
    diet_vars %in% c("Energy(kcal)", "Fat(g)", "Protein(g)","DF/1000kcal", 
                     "Carbohydrate(g)", "Sugars(g)", "Fibre(g)") ~ "Macronutrients",
    diet_vars %in% c("Fat%", "SFA%", "Prot%", "CHO%", "Sugrs%") ~ "Nutrient Percentages",
    diet_vars %in% c("DF/1000kcal", "kcal (%EAR)", "Prot (% RNI)", "Dietary_fibre_%RNI_p") ~ "Nutrient Ratios",
    TRUE ~ "Other Variables"
  )
  
  annotation_df <- data.frame(
    Diet_Category = diet_groups,
    row.names = diet_vars
  )
  
  return(annotation_df)
}

# Function to create basic heatmap with pheatmap
create_basic_heatmap <- function(correlation_data, 
                                 title = "Diet-Urine Correlations Heatmap",
                                 filename = "basic_correlation_heatmap.svg") {
  
  # Create correlation matrix
  corr_matrix <- create_correlation_matrix(correlation_data)
  
  # Remove rows/columns with all NAs
  corr_matrix <- corr_matrix[rowSums(!is.na(corr_matrix)) > 0, 
                             colSums(!is.na(corr_matrix)) > 0]
  
  # Check if we have enough data for clustering
  if (nrow(corr_matrix) < 2 || ncol(corr_matrix) < 2) {
    cat("Not enough data for clustering. Creating simple heatmap without clustering.\n")
    
    # Create simple heatmap without clustering
    svg(filename, width = 12, height = 8)
    pheatmap(
      corr_matrix,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-0.35, 0.3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 10,
      fontsize_col = 10,
      na_col = "grey90",
      main = title,
      border_color = "white",
      cellwidth = 25,
      cellheight = 15
    )
    dev.off()
    return(corr_matrix)
  }
  
  # For clustering, we need to handle NAs by replacing them with 0 or removing incomplete cases
  # Option 1: Replace NAs with 0 (assuming no correlation)
  corr_matrix_clean <- corr_matrix
  corr_matrix_clean[is.na(corr_matrix_clean)] <- 0
  
  # Create annotations
  diet_annotations <- create_diet_annotations(rownames(corr_matrix_clean))
  
  # Define colors
  annotation_colors <- list(
    Diet_Category = c(
      "Macronutrients" = "#e74c3c",
      "Nutrient Percentages" = "#e74c3c",
      "Nutrient Ratios" = "#e74c3c",
      "Other Variables" = "#e74c3c"
    )
  )
  
  # Create the heatmap
  svg(filename, width = 12, height = 8)
  
  tryCatch({
    pheatmap(
      corr_matrix_clean,
      color = colorRampPalette(c("#2c3e50", "white", "#e74c3c"))(100),
      breaks = seq(-0.35, 0.3, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      annotation_row = diet_annotations,
      annotation_colors = annotation_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 10,
      fontsize_col = 10,
      na_col = "grey90",
      main = title,
      border_color = "white",
      cellwidth = 25,
      cellheight = 15
    )
  }, error = function(e) {
    cat("Clustering failed, creating heatmap without clustering\n")
    # Fallback: create heatmap without clustering
    pheatmap(
      corr_matrix,  # Use original matrix with NAs
      color = colorRampPalette(c("#2c3e50", "white", "#e74c3c"))(100),
      breaks = seq(-0.35, 0.3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = diet_annotations,
      annotation_colors = annotation_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 10,
      fontsize_col = 10,
      na_col = "grey90",
      main = paste(title, "(No Clustering)"),
      border_color = "white",
      cellwidth = 25,
      cellheight = 15
    )
  })
  
  dev.off()
  
  cat("Basic heatmap saved to:", filename, "\n")
  return(corr_matrix)
}

# Function to create heatmap with significance annotations (text overlay)
create_annotated_heatmap <- function(correlation_data,
                                     title = "Diet-Urine Correlations with Significance",
                                     filename = "annotated_correlation_heatmap.svg") {
  
  # Create matrices
  corr_matrix <- create_correlation_matrix(correlation_data)
  sig_matrix <- create_significance_matrix(correlation_data)
  
  # Remove rows/columns with all NAs
  valid_rows <- rowSums(!is.na(corr_matrix)) > 0
  valid_cols <- colSums(!is.na(corr_matrix)) > 0
  
  corr_matrix <- corr_matrix[valid_rows, valid_cols]
  sig_matrix <- sig_matrix[valid_rows, valid_cols]
  
  # For clustering, replace NAs with 0
  corr_matrix_clean <- corr_matrix
  corr_matrix_clean[is.na(corr_matrix_clean)] <- 0
  
  # Create annotations
  diet_annotations <- create_diet_annotations(rownames(corr_matrix_clean))
  
  # Define colors
  annotation_colors <- list(
    Diet_Category = c(
      "Macronutrients" = "#e74c3c",
      "Nutrient Percentages" = "#e74c3c",
      "Nutrient Ratios" = "#e74c3c",
      "Other Variables" = "#e74c3c"
    )
  )
  
  # Create the heatmap with text annotations
  svg(filename, width = 25, height = 10)
  
  tryCatch({
    pheatmap(
      corr_matrix_clean,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-0.35, 0.3, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      annotation_row = diet_annotations,
      annotation_colors = annotation_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 12,
      fontsize_col = 12,
      na_col = "grey90",
      main = paste(title, "\n(* p<0.05, ** p<0.01, *** p<0.001)"),
      border_color = "white",
      cellwidth = 30,
      cellheight = 18,
      display_numbers = sig_matrix,  # This adds the significance stars
      number_color = "white",
      fontsize_number = 10
    )
  }, error = function(e) {
    cat("Clustering failed, creating annotated heatmap without clustering\n")
    pheatmap(
      corr_matrix,  # Use original matrix
      color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
      breaks = seq(-0.35, 0.3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = diet_annotations,
      annotation_colors = annotation_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 12,
      fontsize_col = 12,
      na_col = "grey90",
      main = paste(title, "\n(* p<0.05, ** p<0.01, *** p<0.001) - No Clustering"),
      border_color = "white",
      cellwidth = 30,
      cellheight = 18,
      display_numbers = sig_matrix,
      number_color = "white",
      fontsize_number = 10
    )
  })
  
  dev.off()
  
  cat("Annotated heatmap saved to:", filename, "\n")
  return(list(correlation_matrix = corr_matrix, significance_matrix = sig_matrix))
}

# Function to create a simple numerical heatmap
create_numerical_heatmap <- function(correlation_data,
                                     title = "Diet-Urine Correlations (Values Displayed)",
                                     filename = "numerical_correlation_heatmap.pdf") {
  
  # Create correlation matrix
  corr_matrix <- create_correlation_matrix(correlation_data)
  
  # Remove rows/columns with all NAs
  corr_matrix <- corr_matrix[rowSums(!is.na(corr_matrix)) > 0, 
                             colSums(!is.na(corr_matrix)) > 0]
  
  # Check if we have enough data
  if (nrow(corr_matrix) < 2 || ncol(corr_matrix) < 2) {
    cat("Not enough data for clustering in numerical heatmap. Creating simple version.\n")
    
    # Simple version without clustering
    pdf(filename, width = 25, height = 10)
    pheatmap(
      corr_matrix,
      color = colorRampPalette(c("#2c3e50", "white", "#e74c3c"))(100),
      breaks = seq(-0.3, 0.3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 9,
      fontsize_col = 9,
      na_col = "grey90",
      main = paste(title, "(No Clustering)"),
      border_color = "white",
      cellwidth = 35,
      cellheight = 20,
      display_numbers = TRUE,
      number_color = "black",
      number_format = "%.2f",
      fontsize_number = 8
    )
    dev.off()
    return(corr_matrix)
  }
  
  # For clustering, replace NAs with 0
  corr_matrix_clean <- corr_matrix
  corr_matrix_clean[is.na(corr_matrix_clean)] <- 0
  
  # Additional check: ensure no infinite values
  corr_matrix_clean[is.infinite(corr_matrix_clean)] <- 0
  
  # Check for valid correlations (all values should be between -1 and 1)
  corr_matrix_clean[corr_matrix_clean > 1] <- 1
  corr_matrix_clean[corr_matrix_clean < -1] <- -1
  
  # Create annotations
  diet_annotations <- create_diet_annotations(rownames(corr_matrix_clean))
  
  # Define colors
  annotation_colors <- list(
    Diet_Category = c(
      "Macronutrients" = "#e74c3c",
      "Nutrient Percentages" = "#3498db",
      "Nutrient Ratios" = "#f39c12",
      "Other Variables" = "#95a5a6"
    )
  )
  
  # Create the heatmap with correlation values displayed
  pdf(filename, width = 14, height = 10)
  
  tryCatch({
    # Try with clustering first
    pheatmap(
      corr_matrix_clean,
      color = colorRampPalette(c("#2c3e50", "white", "#e74c3c"))(100),
      breaks = seq(-0.3, 0.3, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      annotation_row = diet_annotations,
      annotation_colors = annotation_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 9,
      fontsize_col = 9,
      na_col = "grey90",
      main = title,
      border_color = "white",
      cellwidth = 35,
      cellheight = 20,
      display_numbers = corr_matrix,  # Show original correlation values (with NAs)
      number_color = "black",
      number_format = "%.2f",
      fontsize_number = 8
    )
  }, error = function(e) {
    cat("Clustering failed in numerical heatmap, creating version without clustering\n")
    cat("Error was:", e$message, "\n")
    
    # Fallback: create heatmap without clustering
    pheatmap(
      corr_matrix,  # Use original matrix
      color = colorRampPalette(c("#2c3e50", "white", "#e74c3c"))(100),
      breaks = seq(-0.35, 0.3, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = diet_annotations,
      annotation_colors = annotation_colors,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 9,
      fontsize_col = 9,
      na_col = "grey90",
      main = paste(title, "(No Clustering)"),
      border_color = "white",
      cellwidth = 35,
      cellheight = 20,
      display_numbers = TRUE,
      number_color = "black",
      number_format = "%.2f",
      fontsize_number = 8
    )
  })
  
  dev.off()
  
  cat("Numerical heatmap saved to:", filename, "\n")
  return(corr_matrix)
}

# Function to create summary statistics
create_heatmap_summary <- function(correlation_data) {
  
  # Summary by diet variable
  diet_summary <- correlation_data %>%
    group_by(diet_column) %>%
    summarise(
      n_correlations = n(),
      mean_abs_correlation = mean(abs(correlation)),
      max_abs_correlation = max(abs(correlation)),
      min_p_value = min(p_value),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_abs_correlation))
  
  # Summary by urine variable
  urine_summary <- correlation_data %>%
    group_by(urine_column) %>%
    summarise(
      n_correlations = n(),
      mean_abs_correlation = mean(abs(correlation)),
      max_abs_correlation = max(abs(correlation)),
      min_p_value = min(p_value),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_abs_correlation))
  
  cat("\n=== HEATMAP SUMMARY STATISTICS ===\n")
  cat("Diet variables with strongest correlations:\n")
  print(head(diet_summary, 5))
  
  cat("\nUrine variables with strongest correlations:\n")
  print(head(urine_summary, 5))
  
  return(list(diet_summary = diet_summary, urine_summary = urine_summary))
}

# Main function to create all heatmaps (simplified version)
create_simple_heatmaps <- function(excel_file = "selective_correlation_analysis_20260109_120619.xlsx") {
  
  cat("=== CREATING SIMPLE HIERARCHICAL CLUSTERED HEATMAPS ===\n\n")
  
  # Load and prepare data
  correlation_data <- load_and_prepare_heatmap_data(excel_file)
  
  if (nrow(correlation_data) == 0) {
    cat("No significant correlations found (excluding Days variable)\n")
    return(invisible())
  }
  
  # Create basic heatmap
  cat("\nCreating basic heatmap...\n")
  corr_matrix <- create_basic_heatmap(correlation_data)
  
  # Create annotated heatmap with significance stars
  cat("Creating annotated heatmap with significance stars...\n")
  annotated_results <- create_annotated_heatmap(correlation_data)
  
  # Create numerical heatmap
  cat("Creating numerical heatmap...\n")
  numerical_matrix <- create_numerical_heatmap(correlation_data)
  
  # Create summary statistics
  summary_stats <- create_heatmap_summary(correlation_data)
  
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Files created:\n")
  cat("- Basic heatmap: basic_correlation_heatmap.pdf\n")
  cat("- Annotated heatmap: annotated_correlation_heatmap.pdf\n")
  cat("- Numerical heatmap: numerical_correlation_heatmap.pdf\n")
  
  return(list(
    correlation_data = correlation_data,
    correlation_matrix = annotated_results$correlation_matrix,
    significance_matrix = annotated_results$significance_matrix,
    summary_stats = summary_stats
  ))
}

# ===============================================================================
# EXECUTION
# ===============================================================================

cat("
=== SIMPLE HIERARCHICAL CLUSTERED HEATMAP SCRIPT ===

This simplified script works with your R installation and creates:
- Basic heatmap with clustering and annotations
- Annotated heatmap with significance stars overlay  
- Numerical heatmap showing correlation values

Features:
- Uses only pheatmap (no ComplexHeatmap dependency)
- Excludes 'Days' variable as requested
- Keeps original urine variable names
- Hierarchical clustering on rows and columns
- Diet category annotations (color-coded)
- Significance levels as text overlay
- UPDATED P-VALUE RANGES: * p<0.05, ** p<0.01, *** p<0.001

To run:
results <- create_simple_heatmaps()

This creates 3 PDF files:
- basic_correlation_heatmap.pdf
- annotated_correlation_heatmap.pdf  
- numerical_correlation_heatmap.pdf
")

# Uncomment to run:
results <- create_simple_heatmaps()
# Simple Scatter Plot for Specific Variables
# Just plots one diet variable vs one urine variable

library(ggplot2)

# Function to create a single scatter plot with auto-named output files
create_scatter_plot <- function(diet_file = "IPENS_diet_metadata.csv",
                                urine_file = "food_metadata_level3.csv",
                                diet_var = "Fibre(g)",
                                urine_var = "Complex_Meat",
                                output_file = NULL) {  # Now defaults to NULL for auto-naming
  
  # Create automatic filename if not provided
  if (is.null(output_file)) {
    # Clean variable names for filename (remove special characters)
    diet_clean <- gsub("[^A-Za-z0-9]", "_", diet_var)
    urine_clean <- gsub("[^A-Za-z0-9]", "_", urine_var)
    output_file <- paste0(diet_clean, "_vs_", urine_clean, "_scatter.svg")
  }
  
  # Load data
  diet_data <- read.csv(diet_file, stringsAsFactors = FALSE, check.names = FALSE)
  urine_data <- read.csv(urine_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Merge datasets
  merged_data <- merge(diet_data, urine_data, by = "filename", all = FALSE)
  
  # Extract variables
  x_data <- merged_data[[diet_var]]
  y_data <- merged_data[[urine_var]]
  
  # Remove missing values
  complete_cases <- complete.cases(x_data, y_data)
  x_clean <- x_data[complete_cases]
  y_clean <- y_data[complete_cases]
  
  # Calculate correlation
  cor_test <- cor.test(x_clean, y_clean, method = "pearson")
  
  # Create plot data
  plot_data <- data.frame(x = x_clean, y = y_clean)
  
  # Create significance indicator
  sig_level <- if (cor_test$p.value < 0.001) "***" else 
    if (cor_test$p.value < 0.01) "**" else 
      if (cor_test$p.value < 0.05) "*" else ""
  
  # Statistics text
  stats_text <- paste0(
    "r = ", round(cor_test$estimate, 3), sig_level, "\\n",
    "p = ", format(cor_test$p.value, scientific = TRUE, digits = 2), "\\n",
    "n = ", length(x_clean)
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.5, size = 3, color = "#3498db") +
    geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", fill = "#e74c3c", alpha = 0.2) +
    labs(
      x = diet_var,
      y = urine_var,
      title = paste(diet_var, "vs", urine_var)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16)
    ) +
    annotate("text", 
             x = Inf, y = Inf, 
             label = stats_text,
             hjust = 1.1, vjust = 1.1,
             size = 6,
             fontface = "italic",
             color = "#2c3e50")
  
  # Save plot
  ggsave(output_file, plot = p, width = 8, height = 6, device = "svg")
  
  cat("Scatter plot saved to:", output_file, "\\n")
  cat("Correlation:", round(cor_test$estimate, 3), "\\n")
  cat("P-value:", cor_test$p.value, "\\n")
  
  return(p)
}

# Now your usage will automatically create named files:
create_scatter_plot(diet_var = "Prot%", urine_var = "fleshy fruit")

# Example 2: Energy vs Fruit
# create_scatter_plot(diet_var = "Energy(kcal)", urine_var = "Fruit")

# Example 3: SFA% vs Complex_Meat
# create_scatter_plot(diet_var = "SFA%", urine_var = "Complex_Meat")

cat("
To create a scatter plot, use:

create_scatter_plot(diet_var = 'YOUR_DIET_VARIABLE', 
                   urine_var = 'YOUR_URINE_VARIABLE')

Examples:
create_scatter_plot(diet_var = 'Fibre(g)', urine_var = 'Complex_Meat')
create_scatter_plot(diet_var = 'Energy(kcal)', urine_var = 'Fruit')
create_scatter_plot(diet_var = 'SFA%', urine_var = 'Complex_Meat')
")

