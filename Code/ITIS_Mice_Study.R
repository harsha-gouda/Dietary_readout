
setwd("C:/Users/harsh/OneDrive - University of California, San Diego Health/Projects/Food_readouts/RA_MSV000097333/")

library(tidyverse)
library(tidyr)
library(stringr)
library(mixOmics)
library(ggplot2)
library(tibble)
library(dplyr)
library(pheatmap)
library(ggrepel)

#Please obtain food_metadata.tsv file using the GNPS supported food readout app with taskID given in the supporting information
metadata <- read_csv("metadata_RA_mice_2025.csv")
food_read <- read_tsv("food_metadata.tsv")
food_read <- food_read %>% filter(str_detect(filename, "sample")) %>% arrange(filename)
metadata <- metadata %>% arrange(filename)
m_filter <- metadata %>% 
  filter(str_detect(Timepoint, "Week|Baseline"))
food_read <- food_read %>% filter(filename %in% m_filter$filename)

#PCA analysis#
PCA_whole <- mixOmics::pca(food_read, ncomp = 4, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X, m_filter)
PCA_whole_scores <- PCA_whole_scores %>% filter(str_detect(Timepoint, "Week|Baseline"))
PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "Diet_intervention", alpha = 0.7, size = 3.5,
            title = paste("PCA", "RA_ITIS Diet", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 2),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 2),"%)", sep = ""),
            ggtheme = theme_bw()) +
  stat_ellipse(aes(colour = Diet_intervention), alpha = 1) + 
  theme(plot.title = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20, face = "bold"),
        panel.border = element_rect(color = "black", size = 1.5)) +
  coord_fixed()

print(PCA_plot)

#Volcano Plot

food_metadata <- food_read %>% filter(str_detect(Timepoint, "Week")) # Filter for Week data only

# Get numeric variables (excluding ID columns)
numeric_vars <- names(food_metadata)[sapply(food_metadata, is.numeric)]
numeric_vars <- numeric_vars[!grepl("ID|id", numeric_vars, ignore.case = TRUE)]

# Define baseline and treatment group for comparison
baseline_group <- "Western Diet"
treatment_group <- "Antiinflammatory diet"  # Change this to compare different groups

# Calculate statistics for volcano plot
volcano_data <- data.frame()

for (var in numeric_vars) {
  # Get data for baseline and treatment
  baseline_values <- food_metadata[food_metadata$Diet == baseline_group, var]
  treatment_values <- food_metadata[food_metadata$Diet == treatment_group, var]
  
  # Remove NA values
  baseline_values <- baseline_values[!is.na(baseline_values)]
  treatment_values <- treatment_values[!is.na(treatment_values)]
  
  # Skip if not enough data
  if (length(baseline_values) < 2 || length(treatment_values) < 2) next
  
  # Calculate fold change (log2)
  mean_baseline <- mean(baseline_values)
  mean_treatment <- mean(treatment_values)
  fold_change <- log2(mean_treatment / mean_baseline)
  
  # Perform t-test
  t_test <- t.test(treatment_values, baseline_values)
  p_value <- t_test$p.value
  
  # Store results
  volcano_data <- rbind(volcano_data, data.frame(
    Variable = var,
    log2FC = fold_change,
    pvalue = p_value,
    neg_log10_pvalue = -log10(p_value),
    mean_baseline = mean_baseline,
    mean_treatment = mean_treatment
  ))
}

# Add significance categories
volcano_data$significance <- "Not Significant"
volcano_data$significance[volcano_data$pvalue < 0.05 & volcano_data$log2FC > 0.5] <- "Up ITIS"
volcano_data$significance[volcano_data$pvalue < 0.05 & volcano_data$log2FC < -0.5] <- "Down ITIS"
volcano_data$significance[volcano_data$pvalue < 0.05 & abs(volcano_data$log2FC) <= 0.5] <- "Significant"

# Create volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_pvalue, color = significance)) +
  geom_point(size = 3, alpha = 0.95) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.7) +
  scale_color_manual(values = c("Up ITIS" = "#A2D17D", 
                                "Down ITIS" = "#A75252", 
                                "Significant" = "orange",
                                "Not Significant" = "gray")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(title = paste("Volcano Plot:", treatment_group, "vs", baseline_group),
       x = "Log2 Fold Change",
       y = "-Log10 p-value",
       color = "Significance") +
  geom_text_repel(data = volcano_data[volcano_data$pvalue < 0.05, ], 
                  aes(label = Variable), 
                  size = 3, 
                  box.padding = 0.3,
                  point.padding = 0.3)

print(volcano_plot)

ggsave("ITISvsWD_volcano_plot.svg", plot = volcano_plot, dpi = 300, width = 8, height = 6)

# Print summary table
cat("Summary of significant variables:\n")
significant_vars <- volcano_data[volcano_data$pvalue < 0.05, ]
if(nrow(significant_vars) > 0) {
  print(significant_vars[order(significant_vars$pvalue), c("Variable", "log2FC", "pvalue", "significance")])
} else {
  cat("No significant variables found.\n")
}

