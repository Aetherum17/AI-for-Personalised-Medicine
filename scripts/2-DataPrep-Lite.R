library(rstudioapi) 
library(stringr)
library(dplyr)
library(progeny)
library(survival)
library(glmnet)
library(tidyverse)
library(survminer)

# Set Dynamic Working Directory
directory <- getSourceEditorContext()$path
directory <- str_replace(directory, "/[^/]*$", "")
directory <- str_replace(directory, "/[^/]*$", "")
setwd(directory)

# Find the best threshold to split
best_split_function <- function(values, survival) {
  # Define base threshold
  values_df <- data.frame(Values = values)
  threshold <- min(values_df$Values)
  
  # Create Data frame for storing thresholds and p-values
  threshold_df <- data.frame(Threshold = numeric(),pval = numeric(), stringsAsFactors = FALSE)
  # For each threshold with 0.01 step:
  while(threshold<max(values_df$Values)){
    # Put here to skip the min separation
    threshold = threshold+0.001
    # Calculate the threshold
    values_df$Thresholds <- cut(values_df$Values, breaks = c(min(values_df$Values), threshold, max(values_df$Values)), include.lowest = TRUE)
    # Check if there are at least 2 groups and check if groups are big enough
    if(length(unique(values_df$Thresholds))>1 & table(values_df$Thresholds)[1]>10 & table(values_df$Thresholds)[2]>10){
      # Calculate to obtain p-value
      survdiff_result <- survdiff(survival ~ Thresholds, data = values_df) 
      # Save data to repeat the process 
      threshold_df <- rbind(threshold_df, data.frame(Threshold = threshold, pval = survdiff_result$p)) 
    }
  }
  # Return the threshold with lowest p-value
  best_threshold <- (threshold_df$Threshold[which.min(threshold_df$pval)])
  return(best_threshold)
}

# Read data_mrna_seq_v2_rsem.txt
data_mrna_seq <- read.table(paste(directory, "data/ov_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Drop genes that do not have Hugo_Symbol
data_mrna_seq <- data_mrna_seq[data_mrna_seq$Hugo_Symbol != "" & !is.na(data_mrna_seq$Hugo_Symbol), ]
# Replace NA values with 0
data_mrna_seq[is.na(data_mrna_seq)] <- 0
# Save Row names
rownames(data_mrna_seq) <- make.unique(data_mrna_seq$Hugo_Symbol)
# Drop Entrez_Gene_Id and Hugo_Symbol
data_mrna_seq <- subset(data_mrna_seq, select = -c(Entrez_Gene_Id, Hugo_Symbol))
# Rename columns to match samples
colnames(data_mrna_seq) <- gsub("\\.", "-", colnames(data_mrna_seq))
# Transpose the data_mrna_seq data to match the patient/sample data
data_mrna_seq <- as.data.frame(t(data_mrna_seq))

# Add 1
# Save the data we need
genes_mrna <- c("ATM","CHEK2","HIPK2","MDM2","PPM1D","SIAH1","TP53","WSB1")
print(grep("TP53", colnames(data_mrna_seq), value = TRUE))

data_mrna_seq <- subset(data_mrna_seq, select = genes_mrna)

while(any(data_mrna_seq<0)){
  data_mrna_seq <- data_mrna_seq + 1 
}
data_mrna_seq[data_mrna_seq<0]

# Process Mutations data #######################################################

# Read Data Mutations file
data_mutations <- read.table(paste(directory, "data/ov_tcga_pan_can_atlas_2018/data_mutations.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Subset the columns we might need
data_mutations <- subset(data_mutations, select = c(Hugo_Symbol, Consequence, Tumor_Sample_Barcode))
# Drop genes that do not have Tumor_Sample_Barcode
data_mutations <- data_mutations[data_mutations$Tumor_Sample_Barcode != "" & !is.na(data_mutations$Tumor_Sample_Barcode), ]
# Merge mutations with their Consequence
data_mutations$Hugo_Symbol <- paste(data_mutations$Hugo_Symbol, data_mutations$Consequence, sep = "_")
# Drop Consequence column
data_mutations <- subset(data_mutations, select = -c(Consequence))
# Spread the data set to create a matrix, where columns will be samples, rows mutations and value would indicate if a mutation is present in sample
data_mutations <- data_mutations %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  mutate(value = 1) %>%
  spread(Tumor_Sample_Barcode, value, fill = 0)
# Name the row names
rownames(data_mutations) <- data_mutations$Hugo_Symbol
# Drop the Gene names from column
data_mutations <- subset(data_mutations, select = -c(Hugo_Symbol))
# Transpose the data_mutations data to match the patient/sample data
data_mutations <- as.data.frame(t(data_mutations))
# Save the data we need
tp53_mutations <- c("TP53_frameshift_variant","TP53_frameshift_variant,splice_region_variant","TP53_inframe_deletion","TP53_inframe_deletion,splice_region_variant","TP53_missense_variant","TP53_missense_variant,splice_region_variant","TP53_splice_acceptor_variant","TP53_splice_donor_variant","TP53_splice_region_variant,synonymous_variant","TP53_stop_gained","TP53_stop_gained,frameshift_variant","TP53_stop_gained,splice_region_variant")
print(grep("TP53", colnames(data_mutations), value = TRUE))

# Keep only the Mutations related to p53
data_mutations <- subset(data_mutations, select = tp53_mutations)
# Create a summary mutation column to split the samples
data_mutations$TP53_mut <- ifelse(rowSums(data_mutations) > 0, 1, 0)
# Leave only mutation column
data_mutations <- subset(data_mutations, select = c(TP53_mut))

# Process Patient & Sample Data ################################################

# Read data_clinical_patient.txt and ata_clinical_sample.txt
data_clinical_patient <- read.table(paste(directory, "data/ov_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
data_clinical_sample <- read.table(paste(directory, "data/ov_tcga_pan_can_atlas_2018/data_clinical_sample.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# Select the columns we might need
#data_clinical_patient <- subset(data_clinical_patient, select = c(PATIENT_ID, AGE, SEX, SMOKING_PACK_YEARS, STAGE, OS_MONTHS, OS_STATUS))
data_clinical_patient <- subset(data_clinical_patient, select = c(PATIENT_ID, OS_MONTHS, OS_STATUS))
data_clinical_sample <- subset(data_clinical_sample, select = c(PATIENT_ID, SAMPLE_ID, CANCER_TYPE_DETAILED, TMB_NONSYNONYMOUS))
# Merge Two data frmaes
data_clinical <- merge(data_clinical_sample, data_clinical_patient, by = "PATIENT_ID", all.x = TRUE)
# Drop patients that do not have OS_STATUS or OS_MONTHS
data_clinical <- data_clinical[data_clinical$OS_STATUS != "" & !is.na(data_clinical$OS_STATUS), ]
data_clinical <- data_clinical[data_clinical$OS_MONTHS != "" & !is.na(data_clinical$OS_MONTHS), ]
data_clinical <- data_clinical[((data_clinical$OS_MONTHS > 1 & data_clinical$OS_STATUS == "0:LIVING") | data_clinical$OS_STATUS == "1:DECEASED"), ]
# Surv() function in the {survival} package accepts by default TRUE/FALSE, where TRUE is event and FALSE is censored; 1/0 where 1 is event and 0 is censored
data_clinical$SURV_STATUS <- ifelse(data_clinical$OS_STATUS == "1:DECEASED", 1, 0)

# Select Samples who have mRNA seq and mutations data
data_mrna_seq$SAMPLE_ID <- rownames(data_mrna_seq)
rownames(data_mrna_seq) <- NULL
#data_mutations$SAMPLE_ID <- rownames(data_mutations)
#rownames(data_mutations) <- NULL

#full_sample_list <- Reduce(intersect, list(data_clinical$SAMPLE_ID, data_mrna_seq$SAMPLE_ID, data_mutations$SAMPLE_ID))
full_sample_list <- Reduce(intersect, list(data_clinical$SAMPLE_ID, data_mrna_seq$SAMPLE_ID))

data_clinical <- data_clinical[data_clinical$SAMPLE_ID %in% full_sample_list,]
data_mrna_seq <- data_mrna_seq[data_mrna_seq$SAMPLE_ID %in% full_sample_list,]
#data_mutations <- data_mutations[data_mutations$SAMPLE_ID %in% full_sample_list,]

# Unite Dataframes
data_full <- merge(data_clinical, data_mrna_seq, by = "SAMPLE_ID", all.x = TRUE)
#data_full <- merge(data_full, data_mutations, by = "SAMPLE_ID", all.x = TRUE)

#nb_params_ovarian_tp53_mut <- data_full[data_full$TP53_mut == 1,]
#nb_params_ovarian_tp53_wild <- data_full[data_full$TP53_mut == 0,]

#nb_params_ovarian_tp53_mut_write <- subset(nb_params_ovarian_tp53_mut, select = genes_mrna)
#nb_params_ovarian_tp53_wild_write <- subset(nb_params_ovarian_tp53_wild, select = genes_mrna)
nb_params_ovarian_full_write <- subset(data_full, select = genes_mrna)

#write.table(nb_params_ovarian_tp53_mut_write, file = paste(directory, "data/nb_params_ovarian_tp53_mut.csv", sep = "/"), row.names = FALSE, sep = ",", dec = ".")
#write.table(nb_params_ovarian_tp53_wild_write, file = paste(directory, "data/nb_params_ovarian_tp53_wild.csv", sep = "/"), row.names = FALSE, sep = ",", dec = ".")
write.table(nb_params_ovarian_full_write, file = paste(directory, "data/nb_params_ovarian_full.csv", sep = "/"), row.names = FALSE, sep = ",", dec = ".")

# Cox Regression #################################################################################################
data_clinical_surv <-  Surv(time = data_full$OS_MONTHS, event = data_full$SURV_STATUS)
#data_clinical_surv_tp53_mut <- Surv(time = nb_params_ovarian_tp53_mut$OS_MONTHS, event = nb_params_ovarian_tp53_mut$SURV_STATUS)
#data_clinical_surv_tp53_wild <- Surv(time = nb_params_ovarian_tp53_wild$OS_MONTHS, event = nb_params_ovarian_tp53_wild$SURV_STATUS)

##################################################################################################################
# p53s15DR
p53s15DR_data <- read.table(paste(directory, "results/p53s15DR_tp53.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

for(column in 1:length(colnames(p53s15DR_data)))
{
  formula <- as.formula(paste("data_clinical_surv ~", colnames(p53s15DR_data)[column]))
  fit.coxph <- coxph(formula,
                     data = p53s15DR_data)
  summary_fit <- summary(fit.coxph)
  print(paste("p-val: ", summary_fit$coefficients[, "Pr(>|z|)"], "HZ:", summary_fit$coefficients[, "exp(coef)"]))
  if(summary_fit$coefficients[, "Pr(>|z|)"]<0.05 & summary_fit$coefficients[, "exp(coef)"]!=1){
    print(colnames(p53s15DR_data)[column])
  }
}

fit.coxph <- coxph(data_clinical_surv ~ DDR_55.556,
                   data = p53s15DR_data)
ggforest(fit.coxph, data = p53s15DR_data)

# Copy DDR Data
gg_p53s15DR_data <- p53s15DR_data
# Cut the values of important genes for Hazard Ratio into two groups
gg_p53s15DR_data$DDR_89.000 <- cut(gg_p53s15DR_data$DDR_89.000, breaks = c(min(gg_p53s15DR_data$DDR_89.000), best_split_function(gg_p53s15DR_data$DDR_89.000, data_clinical_surv), max(gg_p53s15DR_data$DDR_89.000)), include.lowest = TRUE)
# Build survival curves for each gene
fit_DDR <- survfit(data_clinical_surv ~ DDR_89.000, data = gg_p53s15DR_data)
# Draw plots for each gene
ggsurvplot(fit_DDR, data = gg_p53s15DR_data, pval = T)

# p53s46DR #####################################################################
p53s46DR_data <- read.table(paste(directory, "results/p53s46DR_tp53.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

for(column in 1:length(colnames(p53s46DR_data)))
{
  formula <- as.formula(paste("data_clinical_surv ~", colnames(p53s46DR_data)[column]))
  fit.coxph <- coxph(formula,
                     data = p53s46DR_data)
  summary_fit <- summary(fit.coxph)
  print(paste("p-val: ", summary_fit$coefficients[, "Pr(>|z|)"], "HZ:", summary_fit$coefficients[, "exp(coef)"]))
  if(summary_fit$coefficients[, "Pr(>|z|)"]<0.05 & summary_fit$coefficients[, "exp(coef)"]!=1){
    print(colnames(p53s46DR_data)[column])
  }
}

fit.coxph <- coxph(data_clinical_surv ~ DDR_89.000,
                   data = p53s46DR_data)
ggforest(fit.coxph, data = p53s46DR_data)

# Copy DDR Data
gg_p53s46DR_data <- p53s46DR_data
# Cut the values of important genes for Hazard Ratio into two groups
gg_p53s46DR_data$DDR_89.000 <- cut(gg_p53s46DR_data$DDR_89.000, breaks = c(min(gg_p53s46DR_data$DDR_89.000), best_split_function(gg_p53s46DR_data$DDR_89.000, data_clinical_surv), max(gg_p53s46DR_data$DDR_89.000)), include.lowest = TRUE)
# Build survival curves for each gene
fit_DDR <- survfit(data_clinical_surv ~ DDR_89.000, data = gg_p53s46DR_data)
# Draw plots for each gene
ggsurvplot(fit_DDR, data = gg_p53s46DR_data, pval = T)

################################################################################
# Mutations split ##############################################################
###################################################################################################################
###################################################################################################################
# p53s15DR TP53 Mut #######################################################################################################
p53s15DR_data <- read.table(paste(directory, "results/p53s15DR_tp53_mut.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

for(column in 1:length(colnames(p53s15DR_data)))
{
  formula <- as.formula(paste("data_clinical_surv_tp53_mut ~", colnames(p53s15DR_data)[column]))
  fit.coxph <- coxph(formula,
                     data = p53s15DR_data)
  summary_fit <- summary(fit.coxph)
  print(paste("p-val: ", summary_fit$coefficients[, "Pr(>|z|)"], "HZ:", summary_fit$coefficients[, "exp(coef)"]))
  if(summary_fit$coefficients[, "Pr(>|z|)"]<0.05 & summary_fit$coefficients[, "exp(coef)"]!=1){
    print(colnames(p53s15DR_data)[column])
  }
}

fit.coxph <- coxph(data_clinical_surv_tp53_mut ~ DDR_1.000,
                        data = p53s15DR_data)
ggforest(fit.coxph, data = p53s15DR_data)

# Copy DDR Data
gg_p53s15DR_data <- p53s15DR_data
# Cut the values of important genes for Hazard Ratio into two groups
gg_p53s15DR_data$DDR_89.000 <- cut(gg_p53s15DR_data$DDR_89.000, breaks = c(min(gg_p53s15DR_data$DDR_89.000), best_split_function(gg_p53s15DR_data$DDR_89.000, data_clinical_surv_tp53_mut), max(gg_p53s15DR_data$DDR_89.000)), include.lowest = TRUE)
# Build survival curves for each gene
fit_DDR <- survfit(data_clinical_surv_tp53_mut ~ DDR_89.000, data = gg_p53s15DR_data)
# Draw plots for each gene
ggsurvplot(fit_DDR, data = gg_p53s15DR_data, pval = T)

# p53s46DR TP 53 Mut ####################################################################################################
p53s46DR_data <- read.table(paste(directory, "results/p53s46DR_tp53_mut.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

for(column in 1:length(colnames(p53s46DR_data)))
{
  formula <- as.formula(paste("data_clinical_surv_tp53_mut ~", colnames(p53s46DR_data)[column]))
  fit.coxph <- coxph(formula,
                     data = p53s46DR_data)
  summary_fit <- summary(fit.coxph)
  print(paste("p-val: ", summary_fit$coefficients[, "Pr(>|z|)"], "HZ:", summary_fit$coefficients[, "exp(coef)"]))
  if(summary_fit$coefficients[, "Pr(>|z|)"]<0.05 & summary_fit$coefficients[, "exp(coef)"]!=1){
    print(colnames(p53s46DR_data)[column])
  }
}

fit.coxph <- coxph(data_clinical_surv_tp53_mut ~ DDR_89.000,
                   data = p53s46DR_data)
ggforest(fit.coxph, data = p53s46DR_data)

# Copy DDR Data
gg_p53s46DR_data <- p53s46DR_data
# Cut the values of important genes for Hazard Ratio into two groups
gg_p53s46DR_data$DDR_89.000 <- cut(gg_p53s46DR_data$DDR_89.000, breaks = c(min(gg_p53s46DR_data$DDR_89.000), best_split_function(gg_p53s46DR_data$DDR_89.000, data_clinical_surv_tp53_mut), max(gg_p53s46DR_data$DDR_89.000)), include.lowest = TRUE)
# Build survival curves for each gene
fit_DDR <- survfit(data_clinical_surv_tp53_mut ~ DDR_89.000, data = gg_p53s46DR_data)
# Draw plots for each gene
ggsurvplot(fit_DDR, data = gg_p53s46DR_data, pval = T)

# p53s15DR TP53 Wild #######################################################################################################
p53s15DR_data <- read.table(paste(directory, "results/p53s15DR_tp53_wild.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

for(column in 1:length(colnames(p53s15DR_data)))
{
  formula <- as.formula(paste("data_clinical_surv_tp53_wild ~", colnames(p53s15DR_data)[column]))
  fit.coxph <- coxph(formula,
                     data = p53s15DR_data)
  summary_fit <- summary(fit.coxph)
  print(paste("p-val: ", summary_fit$coefficients[, "Pr(>|z|)"], "HZ:", summary_fit$coefficients[, "exp(coef)"]))
  if(summary_fit$coefficients[, "Pr(>|z|)"]<0.05 & summary_fit$coefficients[, "exp(coef)"]!=1){
    print(colnames(p53s15DR_data)[column])
  }
}

fit.coxph <- coxph(data_clinical_surv_tp53_wild ~ DDR_89.000,
                   data = p53s15DR_data)
ggforest(fit.coxph, data = p53s15DR_data)

# Copy DDR Data
gg_p53s15DR_data <- p53s15DR_data
# Cut the values of important genes for Hazard Ratio into two groups
gg_p53s15DR_data$DDR_89.000 <- cut(gg_p53s15DR_data$DDR_89.000, breaks = c(min(gg_p53s15DR_data$DDR_89.000), best_split_function(gg_p53s15DR_data$DDR_89.000, data_clinical_surv_tp53_wild), max(gg_p53s15DR_data$DDR_89.000)), include.lowest = TRUE)
# Build survival curves for each gene
fit_DDR <- survfit(data_clinical_surv_tp53_wild ~ DDR_89.000, data = gg_p53s15DR_data)
# Draw plots for each gene
ggsurvplot(fit_DDR, data = gg_p53s15DR_data, pval = T)

# p53s46DR TP 53 Wild ####################################################################################################
p53s46DR_data <- read.table(paste(directory, "results/p53s46DR_tp53_wild.csv", sep = "/"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

for(column in 1:length(colnames(p53s46DR_data)))
{
  formula <- as.formula(paste("data_clinical_surv_tp53_wild ~", colnames(p53s46DR_data)[column]))
  fit.coxph <- coxph(formula,
                     data = p53s46DR_data)
  summary_fit <- summary(fit.coxph)
  print(paste("p-val: ", summary_fit$coefficients[, "Pr(>|z|)"], "HZ:", summary_fit$coefficients[, "exp(coef)"]))
  if(summary_fit$coefficients[, "Pr(>|z|)"]<0.05 & summary_fit$coefficients[, "exp(coef)"]!=1){
    print(colnames(p53s46DR_data)[column])
  }
}

fit.coxph <- coxph(data_clinical_surv_tp53_wild ~ DDR_89.000,
                   data = p53s46DR_data)
ggforest(fit.coxph, data = p53s46DR_data)

# Copy DDR Data
gg_p53s46DR_data <- p53s46DR_data
# Cut the values of important genes for Hazard Ratio into two groups
gg_p53s46DR_data$DDR_89.000 <- cut(gg_p53s46DR_data$DDR_89.000, breaks = c(min(gg_p53s46DR_data$DDR_89.000), best_split_function(gg_p53s46DR_data$DDR_89.000, data_clinical_surv_tp53_wild), max(gg_p53s46DR_data$DDR_89.000)), include.lowest = TRUE)
# Build survival curves for each gene
fit_DDR <- survfit(data_clinical_surv_tp53_wild ~ DDR_89.000, data = gg_p53s46DR_data)
# Draw plots for each gene
ggsurvplot(fit_DDR, data = gg_p53s46DR_data, pval = T)


