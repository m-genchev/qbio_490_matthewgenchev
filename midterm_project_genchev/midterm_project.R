#Queing up necessary Packages and libraries
setwd("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/midterm_project_genchev")
dir.create("outputs")
setwd("outputs")
if (!require(BiocManager)){
library('BiocManager')
}
if (!require(TCGAbiolinks)){
library('TCGAbiolinks')
}
if (!require(maftools)){
library('maftools')
}
if (!require(survival)){
  install.packages("survival") 
}
library(survival)
if (!require(survminer)){
  install.packages("survminer") 
}
library(survminer)
if (!require(ggplot2)){
  install.packages("ggplot2") 
}
library(ggplot2)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

BiocManager::install("DESeq2")
library(DESeq2)

#Pulling GDC clinical data (Patient and radiation info)
clin_query <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Clinical", 
                       file.type = "xml")

GDCdownload(clin_query)

clinic <- GDCprepare_clinic(clin_query, 
                            clinical.info = "patient")
#Have to rename clinical data columns to a format that is readable by the MAF object
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

clinical_rad <- GDCprepare_clinic(query = clin_query,
                                  clinical.info = "radiation")
#Pulling MAF Object
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(maf_query)

maf <- GDCprepare(maf_query)

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinic,
                       isTCGA = TRUE)

#Pulling transcriptomic data
rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

#My question that I want to address has to do with radiation therapy; I want to look
#at whether certain genes are more mutated between radiation treatment patients that
#a complete response vs had diseases that weren't fully cured (Radioactive Progressive Disease, partial response, stable disease).
#I then want to look at the age distribution of these patients and see if there is any difference in age
#I want to look at the survival of the two groups of patients. Finally, I will look at a differential analysis of RNA
#based on the fully recovered vs still diseased patients.

#Data Cleaning and segmenting patients
#Masking for patients that have a measure of response
response_mask <- ifelse(clinical_rad$measure_of_response == "", F, T)
clinical_rad <- clinical_rad[response_mask, ]

#Segmenting Patients into fully recovered and not fully recovered
recovery_mask <- ifelse(clinical_rad$measure_of_response == "Complete Response",
                        "fully_recovered", "not_fully_recovered")
clinical_rad$recovery <- recovery_mask

#Finding Radiation Patients in clinical data
rad_patient <- clinic$Tumor_Sample_Barcode %in% clinical_rad$bcr_patient_barcode
rad_patient_clinic <- clinic [rad_patient, ]

#Combining clinic and clinical_rad by Barcode
clinic$bcr_patient_barcode <- clinic$Tumor_Sample_Barcode
clinic_rad_merge <- merge(clinic, clinical_rad, by = "bcr_patient_barcode", all.x = TRUE)

#Creating a seperate combined clinic and clinical_rad DataFrame that only has radiation patients
status_mask <- ifelse(is.na(clinic_rad_merge$measure_of_response), F, T)
clinic_rad_merge_recovery_only <- clinic_rad_merge[status_mask, ]

#Segmenting the original merged DF so I can compare no rad, fully recovered rad, and not fully recovered rad patient mortality
no_rad_mask <- ifelse(is.na(clinic_rad_merge$measure_of_response) | clinic_rad_merge$measure_of_response == "", "no_radiation",
                      ifelse(clinic_rad_merge$measure_of_response == "Complete Response", "fully_recovered", "not_fully_recovered"))
clinic_rad_merge$recovery <- no_rad_mask

#Creating a boxplot that compares ages of patients that recieved radiation therapy
boxplot(formula= clinic_rad_merge_recovery_only$age_at_initial_pathologic_diagnosis ~ clinic_rad_merge_recovery_only$recovery,
        xlab = 'Radiation Efficacy',
        ylab = 'Age at Diagnosis',
        main = 'Radiation Efficacy vs Age at Diagnosis',
        cex.axis = 0.5)

ggsave("Radiation_Efficacy_vs_Age.jpeg",
       device = "jpeg")

#Looks like the average age for fully recovered patients is higher, interesting
#Let's look at the survival rates for patients that recieved no radiation, those that fully recovered, and those that didn't fully recover

#Need to create a survival time column + recovery column for KM analysis
clinic_rad_merge_KM <- clinic_rad_merge

clinic_rad_merge_KM$survival_time <- ifelse(is.na(clinic_rad_merge_KM$days_to_death),
                                              clinic_rad_merge_KM$survival_time <- clinic_rad_merge_KM$days_to_last_followup,
                                              clinic_rad_merge_KM$survival_time <- clinic_rad_merge_KM$days_to_death)

#Removing -Inf values in survival time
inf_mask <- ifelse(clinic_rad_merge_KM$survival_time == "-Inf", F, T)
clinic_rad_merge_KM <- clinic_rad_merge_KM[inf_mask, ]

#Making a Death Event column + Cleaning NA values
clinic_rad_merge_KM$death_event <- ifelse(clinic_rad_merge_KM$vital_status == "Alive", 
                                          clinic_rad_merge_KM$death_event <- FALSE,
                                          clinic_rad_merge_KM$death_event <- TRUE)
na_mask <- !is.na(clinic_rad_merge_KM$death_event)
clinic_rad_merge_KM <- clinic_rad_merge_KM[na_mask, ]

#Creating a Survminer Object
survival_time <- c(ifelse(clinic_rad_merge_KM$vital_status == "Alive", clinic_rad_merge_KM$days_to_last_followup, clinic_rad_merge_KM$days_to_death))
survival_status <- c(ifelse(clinic_rad_merge_KM$vital_status == "Alive", F, T))
survival_object <- Surv(time = survival_time, event = survival_status)

fit_object <- survfit(survival_object~ recovery, data = clinic_rad_merge_KM)

#Creating the actual KM plot and saving it
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

KM_plot_radiation <- survplot$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16),
                                                         legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_radiation

ggsave("KM_Radiation_Treatment.jpeg",
       device = "jpeg")

#Now, let's look at the most mutated genes in radiation patients which requires the maf_object

#Adding a new column to the maf_object clinical data that looks at radiation therapy status
#Since my radiation data has more patients (including duplicates & unique) I need to get rid of these to create a new column

clinic_rad_merge_KM <- clinic_rad_merge[!duplicated(clinic_rad_merge$Tumor_Sample_Barcode), ]
counts_mask <- clinic_rad_merge_KM$Tumor_Sample_Barcode %in% maf_object@clinical.data$Tumor_Sample_Barcode
clinic_rad_merge_KM <- clinic_rad_merge_KM[counts_mask, ]
recovery_status <- clinic_rad_merge_KM$recovery

#making sure that the dataframes are lined up by Tumor Sample Barcode & Inputting radition treatment data
identical(clinic_rad_merge_KM$Tumor_Sample_Barcode, maf_object@clinical.data$Tumor_Sample_Barcode)

maf_object@clinical.data$recovery_status <- recovery_status

#As not_fully_recovered has a sample size of 10, I will make 3 oncoplots for the radiation treatments
#I need to subset the MAFs to see their unique oncoplots
full_recovery_mask <- ifelse(maf_object@clinical.data$recovery_status == "fully_recovered", T, F)
full_recovery_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[full_recovery_mask]
full_recovery_maf <- subsetMaf(maf = maf_object,
                               tsb = full_recovery_barcodes)

partial_recovery_mask <- ifelse(maf_object@clinical.data$recovery_status == "not_fully_recovered", T, F)
partial_recovery_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[partial_recovery_mask]
partial_recovery_maf <- subsetMaf(maf = maf_object,
                               tsb = partial_recovery_barcodes)

no_rad_treatment_mask <- ifelse(maf_object@clinical.data$recovery_status == "no_radiation", T, F)
no_rad_treatment_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[no_rad_treatment_mask]
no_rad_treatment_maf <- subsetMaf(maf = maf_object,
                               tsb = no_rad_treatment_barcodes)

#Creating the Oncoplots for the three subsetted MAF files
oncoplot(maf = full_recovery_maf,
         top = 5)
ggsave("full_recovery_oncoplot.jpeg",
       device = "jpeg")

oncoplot(maf = partial_recovery_maf,
         top = 5)
ggsave("partial_recovery_oncoplot.jpeg",
       device = "jpeg")

oncoplot(maf = no_rad_treatment_maf,
         top = 5)
ggsave("no_rad_treatment_oncoplot.jpeg",
       device = "jpeg")

#Now, let's look at the transcriptomic data, we need to subset the rna_se data into RNA counts, genes, and clinical dataframes
#Creating a DataFrame that contains clinical data associated with the RNA counts data
rna_clinical <- rna_se@colData
rna_clinical <- as.data.frame(rna_clinical)
#Need to get rid of the columns that aren't actually columns
treatments_mask <- !colnames(rna_clinical) %in% c("treatments", "primary_site", "disease_type")
rna_clinical <- rna_clinical[ , treatments_mask]
#Making the row names informative
row.names(rna_clinical) <- rna_clinical$barcode
#Need to add a radiation treatment column
radiation_patients_mask <- ifelse(clinic_rad_merge$has_radiations_information == "YES", T, F)
radiation_patients_barcode <- clinic_rad_merge$Tumor_Sample_Barcode[radiation_patients_mask]
rna_clinical$radiation_treatment <- ifelse(rna_clinical$patient %in% radiation_patients_barcode,
                                           rna_clinical$radiation_treatment <- "YES",
                                           rna_clinical$radiation_treatment <- "NO")
#Finally, let's clear out any non mutated tissue present
tissue_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_clinical <- rna_clinical[tissue_mask, ]

#Creating a DataFrame that contains the genes associated with the RNA counts
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
#Making the Rownames more informative
row.names(rna_genes) <- rna_genes$gene_id

#Creating a DataFrame that contains the actual RNA counts for patients
rna_counts <- rna_se@assays@data$unstranded
rna_counts <- as.data.frame(rna_counts)
#Making the row and column names more informative
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id
#Getting rid of non-cancerous tissue samples
rna_counts <- rna_counts[ , tissue_mask]

#Setting our radiation_treatment column as a factor for differential expression testing
rna_clinical$radiation_treatment <- factor(rna_clinical$radiation_treatment)

#Turning our confounding variables into factors, for this one it's the stage of breast cancer and gender
rna_clinical$ajcc_pathologic_stage <- factor(rna_clinical$ajcc_pathologic_stage)
rna_clinical$gender <- factor(rna_clinical$gender)

#For some reason my colnames and rownames get messed up so I rename them
colnames(rna_counts) <- rna_clinical$barcode
rownames(rna_counts) <- rna_genes$gene_id

#There are NA values in the ajcc stage and gender, so I will remove those from the dataset
na_mask <- ifelse(is.na(rna_clinical$ajcc_pathologic_stage)
                  | is.na(rna_clinical$gender)
                  , F, T)
rna_clinical <- rna_clinical[na_mask, ]
rna_counts <- rna_counts[ , na_mask]

#Removing any genes with less than 10 counts in total to make the DESeq process faster
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums < 10, F, T)
rna_counts <- rna_counts[low_counts_mask, ]
rna_genes <- rna_genes [low_counts_mask, ]

#Doing the DESeq process
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~ajcc_pathologic_stage + gender + radiation_treatment)

dds_obj <- DESeq(dds)

resultsNames(dds_obj)

results <- results(dds_obj, format = "DataFrame", contrast = c("radiation_treatment", "YES", "NO"))

#Creating a new DF so that we can make a differential expression analysis + volcano plot
# We need useful rows that include the gene name, log2fold change, the p value, the p adjusted value, and -log10 of padj
modified_results <- data.frame(rna_genes$gene_name, 
                               results@rownames, 
                               results@listData$log2FoldChange, 
                               results@listData$pvalue,
                               results@listData$padj,
                               -log10(results@listData$padj))

colnames(modified_results) <- c("gene_name", "ensembl", "log2_fold_change", "p_value", "p_adjusted", "-log10_p_adjusted") ## FIX
row.names(modified_results) <- modified_results$ensembl

#Creating the volcano plot
par(mar=c(1,1,1,1))
EnhancedVolcano(modified_results,
                title = "Gene Expression in Radiation vs Nonradiation Patients",
                lab = modified_results$gene_name,
                x = 'log2_fold_change',
                y = 'p_adjusted',
                labSize = 2,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                max.overlaps = 35,
                shape = 16,
                colAlpha = 0.5)
ggssave("Rad_Patient_Volcano.jpeg",
        device = "jpeg")