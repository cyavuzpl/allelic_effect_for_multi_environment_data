library(vcfR)
library(dplyr)
library(ggplot2)
library(readr)
library(writexl)
library(tidyr)
library(readxl)

# Read input list
input_list <- read_csv("trait_list.csv", show_col_types = FALSE)

# Read signal file
protein_signals <- read_csv("xxxx_combined_signals.csv", show_col_types = FALSE)

# Create SNP column if not present
if (!"SNP" %in% colnames(protein_signals)) {
  if (all(c("CHROM", "POS") %in% colnames(protein_signals))) {
    protein_signals <- protein_signals %>%
      mutate(SNP = paste0(CHROM, "_", POS))
  } else {
    stop("'SNP' or 'CHROM' and 'POS' columns not found in protein_combined_signals.csv.")
  }
}

# Create empty table for results
all_results <- tibble()

# Process each row in the input list
for (row_num in 1:nrow(input_list)) {
  input <- input_list[row_num, ]
  
  vcf_file <- input$vcf_file
  phe_file <- input$phenotype_file
  chr <- as.character(input$chromosome)
  pos <- as.numeric(input$position)
  trait <- input$trait
  outname <- input$Outname
  
  message("Processing SNP: ", chr, ":", pos)
  
  # Extract the specific SNP region to a temporary VCF
  bcf_cmd <- paste("bcftools view -r", paste0(chr, ":", pos, "-", pos), vcf_file, "-Ov -o temp_snp.vcf")
  system(bcf_cmd)
  
  # Skip if the VCF is missing or empty
  if (!file.exists("temp_snp.vcf") || file.info("temp_snp.vcf")$size == 0) {
    warning("temp_snp.vcf is missing or empty. Skipping: ", chr, ":", pos)
    next
  }
  
  # Read VCF
  vcf <- tryCatch({
    read.vcfR("temp_snp.vcf", verbose = FALSE)
  }, error = function(e) {
    warning("VCF read error: ", chr, ":", pos)
    return(NULL)
  })
  
  if (is.null(vcf) || nrow(vcf@fix) == 0) {
    warning("Empty VCF data: ", chr, ":", pos)
    next
  }
  
  vcf_filtered <- vcf[vcf@fix[, "CHROM"] == chr & as.numeric(vcf@fix[, "POS"]) == pos, ]
  if (length(vcf_filtered@fix) == 0) {
    warning("Filtered VCF is empty: ", chr, ":", pos)
    next
  }
  
  # Get REF and ALT alleles
  ref_allele <- vcf_filtered@fix[1, "REF"]
  alt_allele <- vcf_filtered@fix[1, "ALT"]
  
  tidy_vcf <- vcfR2tidy(vcf_filtered, info_only = FALSE)
  gt_df <- tidy_vcf$gt
  
  if (!"Indiv" %in% colnames(gt_df)) {
    if ("sample" %in% colnames(gt_df)) {
      colnames(gt_df)[colnames(gt_df) == "sample"] <- "Indiv"
    } else {
      warning("No 'Indiv' or 'sample' column found in genotype data.")
      next
    }
  }
  
  if (!"gt_GT_alleles" %in% colnames(gt_df)) {
    warning("Column 'gt_GT_alleles' not found.")
    next
  }
  
  snp_tidy <- gt_df
  phenotypes <- read.table(phe_file, header = TRUE)
  colnames(phenotypes)[1] <- "Indiv"
  
  gwas_data <- left_join(snp_tidy, phenotypes, by = "Indiv")
  
  gwas_data <- gwas_data %>%
    mutate(Allele = case_when(
      gt_GT_alleles == paste0(ref_allele, "|", ref_allele) ~ "REF",
      gt_GT_alleles == paste0(alt_allele, "|", alt_allele) ~ "ALT",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Allele))
  
  # Check if there are at least two REF and ALT homozygotes
  allele_counts <- gwas_data %>%
    group_by(Allele) %>%
    summarise(count = n(), .groups = "drop")
  
  if (all(c("REF", "ALT") %in% allele_counts$Allele) &&
      all(allele_counts$count >= 2)) {
    
    trait_stats <- gwas_data %>%
      select(Allele, Value = all_of(trait)) %>%
      group_by(Allele) %>%
      summarise(avg_value = mean(Value, na.rm = TRUE),
                count = n(),
                .groups = "drop")
    
    reference_group <- gwas_data %>% filter(Allele == "REF")
    comparison_group <- gwas_data %>% filter(Allele == "ALT")
    
    percent_change <- (mean(comparison_group[[trait]], na.rm = TRUE) -
                         mean(reference_group[[trait]], na.rm = TRUE)) /
      mean(reference_group[[trait]], na.rm = TRUE) * 100
    
    snp_id <- paste0(chr, "_", pos)
    
    # Get Signal_P value matching both SNP and phenotype file (basename)
    matched_p <- protein_signals %>%
      filter(SNP == snp_id, basename == phe_file) %>%
      pull(P)
    
    matched_p <- ifelse(length(matched_p) == 0, NA, matched_p)
    
    out_data <- tibble(
      SNP = snp_id,
      Trait = trait,
      Reference = ref_allele,
      Comparison = alt_allele,
      Change = percent_change,
      Signal_P = matched_p,
      Outname = outname
    )
    
    all_results <- bind_rows(all_results, out_data)
    
  } else {
    message("Not enough REF and ALT homozygotes: ", chr, ":", pos)
  }
  
  # Clean up temporary VCF
  file.remove("temp_snp.vcf")
}

# Write results to Excel
write_xlsx(all_results, "V2_snp_new_results.xlsx")
