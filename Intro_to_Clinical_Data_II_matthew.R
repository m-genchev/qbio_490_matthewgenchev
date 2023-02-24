#Setting a working directory and checking to make sure that it's been set properly
setwd("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/analysis_data")
getwd()

#Loading in Data and libraries used in the project
library(BiocManager)
library(TCGAbiolinks)
library(maftools)
if (!require(survival)){
 install.packages("survival") 
}
if (!require(survminer)){
  install.packages("survminer") 
}
if (!require(ggplot2)){
  install.packages("ggplot2") 
}
library(ggplot2)
clinic <- read.csv("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/analysis_data/brca_clinical_data.csv")
clin_query <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Clinical", 
                       file.type = "xml")
clinical_drug <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query,
                                   clinical.info = "radiation")

#Exploring columns of the clinical data frame to find an interesting one and testing for NA values
colnames(clinic)
lymphnodeNaCount <- !is.na(clinic$number_of_lymphnodes_positive_by_he)
sum(lymphnodeNaCount)
#I ended up choosing the positive lymph node count by he.
#This is a discrete variable (countable amount of values in between values)

#Exploring clinical_drug data and making sure the selected data doesn't have many NAs
colnames(clinical_drug)
head(clinical_drug)
therapyTypeMask <- !is.na(clinical_drug$therapy_types)
sum(therapyTypeMask)
#This variable is categorical, and I chose the therapy type that the cancer was treated.
#This could include chemotherapy, hormone therapy, or more specialized therapy

#Exploring clinical_rad data and making sure the selected data doesn't have many NAs
colnames(clinical_rad)
radiationDosageMask <- !is.na(clinical_rad$radiation_dosage)
sum(radiationDosageMask)
#This variable is continuous, and I chose the amount of radiation a cancer received.
#This includes radiation values in cGy and Gy

#Hypothesis 1: Patients with higher positive lymph node counts received higher
#radiation values (if any at all) and were more likely to be treated with more
#aggressive forms of therapy (chemotherapy)

#Hypothesis 2: A higher positive lymph node count is related to a worse survival

#Hypothesis 3: Patients who had hormone therapy were more likely to survive than
#those who had chemotherapy

#Creating a merged DataFrame to create a scatter plot comparing lymph node count and
#the amount of radiation a patient was treated with NA values removed
rad_mask <- ifelse(is.na(clinical_rad$radiation_dosage), F, T)
cleaned_rad <- clinical_rad[rad_mask, ]

lymphnode_mask <- ifelse(is.na(clinic$number_of_lymphnodes_positive_by_he), F, T)
cleaned_clinic <- clinic[lymphnode_mask, ]

clinic_rad_merge <- merge(cleaned_clinic, cleaned_rad, by = "bcr_patient_barcode")

#Creating a scatterplot of lymph node count vs radiation received
plot(clinic_rad_merge$number_of_lymphnodes_positive_by_he,clinic_rad_merge$radiation_dosage,
     main = "Radiation Dosage vs Number of Positive Lymph Nodes",
     xlab = "Number of Positive Lymphnodes (By he)",
     ylab = "Radiation Dosage (Gy)")
#The plot is not too informative, as it seems that many times radiation is administered
#Without any lymph nodes being positive. Maybe next time we can remake this plot
#It does seem that the amount of radiation increases with the number of lymph nodes
#Found positive but the correlation is most likely not too strong

#Creating a Kaplan-Meier plot for positive lymph nodes
#Creating categorical variables for lymph nodes (1-10, 11-20, and 21+ possitive)
low_nodes_mask <- ifelse (cleaned_clinic$number_of_lymphnodes_positive_by_he <= 10, T, F)
mid_nodes_mask <- ifelse (cleaned_clinic$number_of_lymphnodes_positive_by_he >10 & cleaned_clinic$number_of_lymphnodes_positive_by_he <= 20, T, F)
high_nodes_masl <- ifelse (cleaned_clinic$number_of_lymphnodes_positive_by_he >10, T, F)

cleaned_clinic$positive_lymph_node_count <- ifelse(low_nodes_mask, "Low", ifelse(mid_nodes_mask, "Medium", "High"))

#Making a survival time column
cleaned_clinic$survival_time <- ifelse (is.na(cleaned_clinic$days_to_death),
                                        cleaned_clinic$survival_time <- cleaned_clinic$days_to_last_followup,
                                        cleaned_clinic$survival_time <- cleaned_clinic$days_to_death)
#Removing any -Inf values in the survival_time column
inf_mask1 <- ifelse(cleaned_clinic$survival_time == "-Inf", F, T)
cleaned_clinic <- cleaned_clinic[inf_mask1, ]

#making a death event (T/F) column
cleaned_clinic$death_event <- ifelse(cleaned_clinic$vital_status == "ALIVE", cleaned_clinic$death_event <- FALSE,
                                     cleaned_clinic$death_event <- TRUE)

#Creating a Survminer Object
survival_time <- c(ifelse(cleaned_clinic$vital_status == "ALIVE", cleaned_clinic$days_to_last_followup, cleaned_clinic$days_to_death))
survival_status <- c(ifelse(cleaned_clinic$vital_status == "ALIVE", F, T))
survival_object <- Surv(time = survival_time, event = survival_status)

fit_object <- survfit(survival_object~ positive_lymph_node_count, data = cleaned_clinic)

#Finally creating the Kaplan-Meier Plot
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

KM_plot_lymphnodes <- survplot$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16),
                                              legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_lymphnodes
#p-value of 0.0027. There isn't much data for medium and high lymph node counts so 
#maybe it would be worth it to further stratify low counts and see if there is a difference
#that is more significant.

#Creating a Kaplan-Meier plot for Therapy Types
#Creating a cleaned DataFrame that contains drug therapy data and clinical data without NA values
therapy_mask <- ifelse(is.na(clinical_drug$therapy_types), F, T)
cleaned_therapy <- clinical_drug[therapy_mask, ]

clinic_drug_merge <- merge(clinic, cleaned_therapy, by = "bcr_patient_barcode")
#Making a survival time column
clinic_drug_merge$survival_time <- ifelse (is.na(clinic_drug_merge$days_to_death),
                                           clinic_drug_merge$survival_time <- clinic_drug_merge$days_to_last_followup,
                                           clinic_drug_merge$survival_time <- clinic_drug_merge$days_to_death)
#Removing any -Inf values in the survival_time column
inf_mask2 <- ifelse(clinic_drug_merge$survival_time == "-Inf", F, T)
clinic_drug_merge <- clinic_drug_merge[inf_mask2, ]

#making a death event (T/F) column
clinic_drug_merge$death_event <- ifelse(clinic_drug_merge$vital_status == "ALIVE", clinic_drug_merge$death_event <- FALSE,
                                        clinic_drug_merge$death_event <- TRUE)

#Creating a Survminer Object
survival_time_drug <- c(ifelse(clinic_drug_merge$vital_status == "ALIVE", clinic_drug_merge$days_to_last_followup, clinic_drug_merge$days_to_death))
survival_status_drug <- c(ifelse(clinic_drug_merge$vital_status == "ALIVE", F, T))
survival_object_drug <- Surv(time = survival_time_drug, event = survival_status_drug)

fit_object_drug <- survfit(survival_object_drug~ therapy_types, data = clinic_drug_merge)

#Finally creating the Kaplan-Meier Plot
survplot <- ggsurvplot(fit_object_drug, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

KM_plot_therapy_types <- survplot$plot + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16),
                                                         legend.title = element_text(size=14), legend.text = element_text(size=12))
KM_plot_therapy_types
#p = 0.038, not really significant. This could be due to noise from the one-off
#therapies, so maybe rerunning with just chemotherapy vs hormone therapy might be significant