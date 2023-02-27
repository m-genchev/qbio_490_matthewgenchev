#Setting working directory and making sure that the correct working directory was set
setwd("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/analysis_data")
getwd()

#Installing libraries that I will need for the project
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

#Reading in the clinical data csv made in part 1
clinical <- read.csv("/Users/matthewgenchev/Desktop/USC/qbio_490_matthewgenchev/analysis_data/brca_clinical_data.csv")

#Initializing maf_object from part 1
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

#Looking at clinical variables to choose
colnames(maf_object@clinical.data)

#I decided to look at immunohistochemistry_level_result which tells if 
#a tumor has high amounts of HER2 receptors

#Checking the factors of immunohistochemistry_level_result
factor(maf_object@clinical.data$her2_immunohistochemistry_level_result)

#The factors are 0, 1+, 2+, 3+. A score of 0 or 1+ is HER2 negative and a score of
#3+ is HER2 positive. 2+ is borderline.
#Creating a new column HER2_Status which segments scores of 0 and 1+ to "NEGATIVE"
#and scores of 3+ to "POSITIVE"
maf_object@clinical.data$HER2_Status <- ifelse(maf_object@clinical.data$her2_immunohistochemistry_level_result == c("0", "1+"),
                                               "NEGATIVE", ifelse(maf_object@clinical.data$her2_immunohistochemistry_level_result == "3+",
                                                      "POSITIVE", NA))
maf_object@clinical.data$HER2_Status

#Creating a cooncoplot for HER2_Status
#Need to start by subsetting HER2 positive and HER2 negative patients in their own MAF objects
HER2_negative_mask <- ifelse(maf_object@clinical.data$HER2_Status == "POSITIVE" | is.na(maf_object@clinical.data$HER2_Status), F, T)

HER2_negative_pt_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[HER2_negative_mask]

HER2_negative_maf <- subsetMaf(maf = maf_object,
                               tsb = HER2_negative_pt_barcodes)

#Subsetting for a HER2 positive maf file
HER2_positive_mask <- ifelse(maf_object@clinical.data$HER2_Status == "NEGATIVE" | is.na(maf_object@clinical.data$HER2_Status), F, T)

HER2_positive_pt_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[HER2_positive_mask]

HER2_positive_maf <- subsetMaf(maf = maf_object,
                               tsb = HER2_positive_pt_barcodes)

#Creating the cooncoplot with the subsetted MAF files

#Have to manually find the top 10 genes that are mutated as the coOncoplot() function does not provide an
#automatoc way to do so
m1.genes = getGeneSummary(HER2_negative_maf)[1:10]
m2.genes = getGeneSummary(HER2_positive_maf)[1:10]
mdt = merge(m1.genes[,.(Hugo_Symbol, MutatedSamples)], m2.genes[,.(Hugo_Symbol, MutatedSamples)], by = 'Hugo_Symbol', all = TRUE)
mdt$MutatedSamples.x[is.na(mdt$MutatedSamples.x)] = 0
mdt$MutatedSamples.y[is.na(mdt$MutatedSamples.y)] = 0
mdt$max = apply(mdt[,.(MutatedSamples.x, MutatedSamples.y)], 1, max)
mdt = mdt[order(max, decreasing = TRUE)]

#Creating a vector that contains the top mutated genes
top_mutated <- c(mdt$Hugo_Symbol)
#This for some reason gave me 15 genes so I copied the top 10 and put them in a new vector
top_10_mutated_genes <- c("PIK3CA", "TP53", "TTN", "CDH1", "MAP3K1", "MUC16", "KMT2C", "GATA3", "HMCN1", "FAT3")

#Plotting the cooncoplot with the top 10 mutated genes
coOncoplot(m1 = HER2_negative_maf,
           m2 = HER2_positive_maf,
           genes = top_10_mutated_genes,
           m1Name = "Negative HER2 Receptor Status (IHC Scores of 0, 1+)",
           m2Name = "Positive HER2 Receptor Status (IHC Score of 3+)",
           borderCol = NA)

#I chose gene MAP3K1 to examine as it is mutated 2% of the time in positive HER2 receptors
#and 11% of the time in negative HER2 receptors.
#MAP3K1 is a kinase that is involved in the ERK pathway which is important in cell proliferation and growth
#The ERK pathway is also used as an alternative pathway for cells to grow that become resistant to HER2+ treatment.
#Possibly, these patients have been treated with an anti-HER2+ drug which caused the tumor cells to grow
#using the ERK pathway which requires an intact MAP3K1 protein which could explain why it's less mutated in
#cancers that are are positive for HER2 receptors as opposed to cancers which are negative.

#Creating a contingency table for MAP3K1 and HER2_status
#Need to create a new column with 2 levels for if a patient has a mutated or normal MAP3K1 gene
maf_object@data$MAP3K1_status <- ifelse(maf_object@data$Hugo_Symbol == "MAP3K1", "Mutated", "Normal")
colnames(maf_object@data)

#Creating new dataframes that only subset $MAP3K1_status or $HER2_Status alogn with Tumor_Sample_Barcode
#so that I can later combine them into a new DataFrame by $Tumor_Sample_Barcode so that I have the same number
#of rows to eventually create a contingency table and a fischer exact test
MAP3K1_df <- maf_object@data[ , c("MAP3K1_status", "Tumor_Sample_Barcode") ]
HER2_df <- maf_object@clinical.data[ , c("HER2_Status", "Tumor_Sample_Barcode")]

#Merged DataFrame by Tumor_Sample_Barcode
gene_receptor_merge <- merge(MAP3K1_df, HER2_df,
                             by = "Tumor_Sample_Barcode")

#Creating a contingency table and mosaic plotting
contig <- table(gene_receptor_merge$HER2_Status, gene_receptor_merge$MAP3K1_status)
mosaicplot(contig)

#Running a Fisher Exact Test on the contingency table
fisher_test <- fisher.test(contig)
fisher_test
#The fisher test yielded a p value of 0.001281 and an odds ratio of 5.381592
#This means that it is extremely statistically likely that HER2 status and MAP3K1 tumor mutations
#are correlated through the p-value. The odds ratio tells us that the odds of being HER2 negative
#and having a mutated MAP3K1 gene are 5.38x as likely as being HER2 positive and having a mutated MAP3K1 gene

#Creating a colollipop plot for differing HER2 status
lollipopPlot2(m1 = HER2_negative_maf,
              m2 = HER2_positive_maf,
              m1_name = "Negative HER2 Receptor Status (IHC Scores of 0, 1+)",
              m2_name = "Positive HER2 Receptor Status (IHC Score of 3+)",
              gene = "MAP3K1",
              )

#Creating a boolean vector for overall survival status
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$vital_status == 'Alive',
                                                           T, F)
#Creating a mafSurvival plot
mafSurvival(maf = maf_object,
            genes = "MAP3K1",
            time = "days_to_last_followup",
            Status = "Overall_Survival_Status",
            isTCGA = TRUE)
#There is a difference, a non-mutated MAP3K1 gene is able to control the cell cycle
#better than a mutated one, so it is likely that patients with a mutated MAP3K1 gene
#will not survive as long as those with a fully functioning one