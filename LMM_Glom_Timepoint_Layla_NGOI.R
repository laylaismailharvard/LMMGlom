
#LMM for Glom MTS 2,3,5 

# Load required package (stats is available by default in R)
library(stats)

# Read the data file


mts.2.3.5.glom <- read.csv("~/Desktop/summer25/mts.2.3.5.glom.quantilenorm.data.csv",
                           row.names = 1)
head(mts.2.3.5.glom)

treatment_factor <- sub("^[^_]+_([^_]+)_.*$", "\\1", colnames(mts.2.3.5.glom))

keep_columns <- treatment_factor =="None"

mts.2.3.5.glom <- mts.2.3.5.glom[, keep_columns]

data<-mts.2.3.5.glom

colnames(mts.2.3.5.glom)


############ For comparisions on Timepoint###########################################################


timepoint_factor <- factor(sub(".*_T(\\d+)_.*", "T\\1", colnames(mts.2.3.5.glom)),
                           levels = c("T0", "T16"))

# Set "None" as the reference level
timepoint_factor <- relevel(timepoint_factor, ref = "T0")  # Reference level set to "None"

print(levels(timepoint_factor))  # Verify levels are assigned correctly
# Expected output: [1] "T0"  "T16"

# Function to test linear hypothesis (contrast) and calculate linear fold change for one gene
perform_contrast_and_linear_fc <- function(gene_data, treatment_factor, contrast) {
  # Fit linear model: gene expression ~ treatment factor
  model <- lm(gene_data ~ treatment_factor)
  
  # Extract estimated coefficients
  beta <- coef(model)
  
  # Calculate linear fold change
  fold_change <- sum(contrast * beta)  # Linear fold change
  
  # Extract covariance matrix of coefficients
  V <- vcov(model)
  
  # Ensure the contrast is a column vector
  contrast <- as.matrix(contrast)

  # Compute estimated contrast: C %*% beta
  est <- sum(contrast * beta)

  # Compute standard error of the contrast: sqrt(C %*% V %*% t(C))
  se <- sqrt(t(contrast) %*% V %*% contrast)

  # Compute t statistic
  t_stat <- est / se

  # Degrees of freedom from residuals
  df <- model$df.residual

  # Compute two-sided p-value from t-distribution
  p_value <- 2 * pt(-abs(t_stat), df)

  return(list(linear_fc = fold_change, p_value = p_value))
}

# Define the contrasts based on the new levels
contrast_T16_vs_T0 <- c(0, 1)       # T16 vs. T0

# For each gene, calculate linear fold change and p-values for each contrast
results_T0_vs_T16 <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, timepoint_factor, contrast_T16_vs_T0))

# Extract p-values and linear fold changes for each contrast
pvals_T0_vs_T16 <- sapply(results_T0_vs_T16, function(x) x$p_value)
linear_fc_T0_vs_T16 <- sapply(results_T0_vs_T16, function(x) x$linear_fc)


# Apply FDR correction (Benjamini–Hochberg) for each set of p-values
fdr_T0_vs_T16 <- p.adjust(pvals_T0_vs_T16, method = "BH")


# Combine results into a data frame, including gene names and log2FC calculations
results <- data.frame(
  Gene = rownames(mts.2.3.5.glom),
  T16_vs_T0_FDR = fdr_T0_vs_T16,
  T16_vs_T0_linearFC = linear_fc_T0_vs_T16,
  T16_vs_T0_log2FC = log2(abs(linear_fc_T0_vs_T16)),
  
  row.names = NULL
)

# Inspect first few rows of the results
head(results)
#> head(results)
#Gene T16_vs_T0_FDR T16_vs_T0_linearFC T16_vs_T0_log2FC
#  A2M  9.257025e-05        -358.213748         8.484677
#  ACADM  4.605004e-02          -4.811989         2.266633
#  ACAT1  2.395027e-05         -10.094643         3.335518

# Set significance thresholds
fdr_threshold <- 0.05

# Add a "Significance" column based on FDR and log2FC thresholds
results$Significance <- with(results, 
                             (T16_vs_T0_FDR < fdr_threshold))

# View the updated results with the "Significance" column
head(results)
#> head(results)
#Gene T16_vs_T0_FDR T16_vs_T0_linearFC T16_vs_T0_log2FC Significance
#    A2M  9.257025e-05        -358.213748         8.484677         TRUE
#  ACADM  4.605004e-02          -4.811989         2.266633         TRUE
#  ACAT1  2.395027e-05         -10.094643         3.335518         TRUE

# Save results to an Excel file
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)
write.xlsx(results, file = "Results_with_T16_vs_T0_glom_noneonly.xlsx")
#results<-readxl::read_xlsx("C:/Users/Laptop Admin/Desktop/Slice culture data/Results_with_T16_vs_T0_glom_noneonly.xlsx")



ngoi <- results[results$Significance == FALSE, "Gene"]



#heatmap

glom.ngois<-mts.2.3.5.glom[ngoi,]

Annot_df <- data.frame(
  Timepoint = character(length = ncol(glom.ngois)),
  Treatment = character(length = ncol(glom.ngois)),
  DonorID = character(length = ncol(glom.ngois)),
  GlomLocation = character(length = ncol(glom.ngois)),
  Histology = character(length = ncol(glom.ngois))
)

Annot_df$Timepoint<-sub(".*_(T[0-9]+)_.*", "\\1",colnames(glom.ngois))
Annot_df$Treatment<-sub(".*_(.*?)_.*", "\\1", colnames(glom.ngois))
Annot_df$DonorID<-sub("^([A-Za-z0-9]+)_.*", "\\1",  colnames(glom.ngois))
Annot_df$GlomLocation <- sub("^[^_]+_[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(glom.ngois))
Annot_df$Histology <- sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(glom.ngois))



row.names(Annot_df) <- colnames(glom.ngois)


name.treatments <- function(timepoint, treatment) {
  if (timepoint == "T0" & treatment == "None") { 
    return("Naive")
  } else if (timepoint == "T16" & treatment == "None") {
    return("Vehicle")
  } else if (timepoint == "T16" & treatment == "IFN") {
    return("IFN")
  } else if (timepoint == "T16" & treatment == "LPS") {
    return("LPS")
  } else if (timepoint == "T16" & treatment == "PIC") {
    return("PIC")
  } else if (timepoint == "T16" & treatment == "Baricitinib") {
    return("Baricitinib")
  } else if (timepoint == "T16" & treatment == "Baricitinib+IFN") {
    return("Baricitinib+IFN")
  } else {
    return(NA)
  }
}



#Apply to data
Annot_df$Treatment <- mapply(name.treatments, Annot_df$Timepoint, Annot_df$Treatment)

desired_order <- c("Naive", "Vehicle", "IFN", "LPS", "PIC")




# Ensure that 'Treatment' column is a factor with the desired order
Annot_df$Treatment <- factor(Annot_df$Treatment, levels = desired_order)
Annot_df$GlomLocation <- factor(Annot_df$GlomLocation, levels = c("Edge","Central"))
Annot_df$DonorID <- factor(Annot_df$DonorID, levels = c("MTS002","MTS003","MTS005"))
# Sort Annot_df based on Treatment and GlomLocation
sorted_Annot_df <- Annot_df[order(Annot_df$Treatment, Annot_df$GlomLocation,Annot_df$DonorID), ]

# View the sorted annotated dataframe
head(sorted_Annot_df)



ordered.glom.data<-glom.ngois[,order(Annot_df$Treatment, Annot_df$GlomLocation)]

# Check if the columns match
identical(colnames(ordered.glom.data), rownames(sorted_Annot_df)) # Should return TRUE

my.breaks <- c(seq(-2, 2, by=0.1)) 

my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("white", "orange", "red"))(length(my.breaks)/2))


library(pheatmap)
#pheatmap(ordered.glom.data, color=my.colors, cellwidth=3, 
#        cellheight=0.02, lwd=0.5, border_color="darkgrey", fontsize_col=5, fontsize_row=5, scale="row", 
#        breaks= my.breaks, treeheight_row=10,  show_rownames = FALSE, show_colnames = TRUE, 
#        annotation_col = sorted_Annot_df,cluster_cols=FALSE, cluster_rows=TRUE, clustering_distance_rows="euclidean")



#save it PNG
png("Glom.T16vT0.heatmap.png", width = 1200, height = 1500, res = 150)

heatmap_plot <- pheatmap(ordered.glom.data,
                         color = my.colors,
                         cellwidth = 3,
                         cellheight = 0.02,
                         lwd = 0.5,
                         border_color = "darkgrey",
                         fontsize_col = 5,
                         fontsize_row = 5,
                         scale = "row",
                         breaks = my.breaks,
                         treeheight_row = 10,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         annotation_col = sorted_Annot_df,
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean")

# force the plot to render
print(heatmap_plot)

dev.off()



# ---- START: Identify Upregulated Genes in T16 vs T0 ----

results.test<-results

# Add signed log2 fold change (preserve direction)
results.test$Signed_log2FC <- log2(abs(results$T16_vs_T0_linearFC)) * sign(results$T16_vs_T0_linearFC)

# Set FDR threshold
fdr_threshold <- 0.05

# Add Significance and Direction columns
results.test$Sig <- results.test$T16_vs_T0_FDR < fdr_threshold
results.test$Direction <- ifelse(results.test$Significance,
                            ifelse(results.test$Signed_log2FC > 0, "Up in T16", "Down in T16"),
                            "Not Significant")

# Subset: Upregulated genes in T16
upregulated_in_T16 <- subset(results.test, Direction == "Up in T16")
upregulated_in_T0<-subset(results.test, Direction == "Down in T16")

# View top results (optional)
head(upregulated_in_T16)

# Save to CSV
write.csv(upregulated_in_T16, "Upregulated_GLom_Genes-T16_T16_vs_T0.csv", row.names = FALSE)
write.csv(upregulated_in_T0, "Upregulated_Glom_Genes-T0_T16_vs_T0.csv", row.names = FALSE)
write.csv(results.test, "results.T16vT0..gloms.allgene+sig.csv", row.names = FALSE)
#################################
#Volcanos 

##################################
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)


res1<-results.test
res1$T16_vs_T0_pvals<-pvals_T0_vs_T16
res1$Signed_log2FC<-results.test$Signed_log2FC



# Set output file first
pdf(file = "MTS2,3,5-Quantile-T0-vs-T16-.None.Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_log2FC', # log2FC
                         y = 'T16_vs_T0_pvals', # p-value
                         title = 'T0 vs. T16 Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)




print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()


glom.res1<-res1[,c(1:6,8:9)]
glom.upndown.sig.genes<-ngoi
#write.csv(upregulated_in_T16, "Upregulated_GLom_Genes-T16_T16_vs_T0.csv", row.names = FALSE)

###############################################################################################################
                                 #TUBULES    #didnt need to edit because there is already only "none" treatment 

################################################TUBULES################################################

# Read the data file

data <- read.csv("C:/Users/Laptop Admin/Desktop/Slice culture data/Quantile Normalized datasets/mts.2.3.5.tub.quantilenorm.data.csv",
                 row.names = 1)
mts.2.3.5.tub<-data 

############ For comparisions on Timepoint###########################################################


timepoint_factor <- factor(sub(".*_T(\\d+)_.*", "T\\1", colnames(mts.2.3.5.tub)),
                           levels = c("T0", "T16"))

# Set "None" as the reference level
timepoint_factor <- relevel(timepoint_factor, ref = "T0")  # Reference level set to "None"

print(levels(timepoint_factor))  # Verify levels are assigned correctly
# Expected output: [1] "T0"  "T16"

# Function to test linear hypothesis (contrast) and calculate linear fold change for one gene
perform_contrast_and_linear_fc <- function(gene_data, treatment_factor, contrast) {
  # Fit linear model: gene expression ~ treatment factor
  model <- lm(gene_data ~ treatment_factor)
  
  # Extract estimated coefficients
  beta <- coef(model)
  
  # Calculate linear fold change
  fold_change <- sum(contrast * beta)  # Linear fold change
  
  # Extract covariance matrix of coefficients
  V <- vcov(model)
  
  # Ensure the contrast is a column vector
  contrast <- as.matrix(contrast)
  
  # Compute estimated contrast: C %*% beta
  est <- sum(contrast * beta)
  
  # Compute standard error of the contrast: sqrt(C %*% V %*% t(C))
  se <- sqrt(t(contrast) %*% V %*% contrast)
  
  # Compute t statistic
  t_stat <- est / se
  
  # Degrees of freedom from residuals
  df <- model$df.residual
  
  # Compute two-sided p-value from t-distribution
  p_value <- 2 * pt(-abs(t_stat), df)
  
  return(list(linear_fc = fold_change, p_value = p_value))
}

# Define the contrasts based on the new levels
contrast_T16_vs_T0 <- c(0, 1)       # T16 vs. T0

# For each gene, calculate linear fold change and p-values for each contrast
results_T0_vs_T16 <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, timepoint_factor, contrast_T16_vs_T0))

# Extract p-values and linear fold changes for each contrast
pvals_T0_vs_T16 <- sapply(results_T0_vs_T16, function(x) x$p_value)
linear_fc_T0_vs_T16 <- sapply(results_T0_vs_T16, function(x) x$linear_fc)


# Apply FDR correction (Benjamini–Hochberg) for each set of p-values
fdr_T0_vs_T16 <- p.adjust(pvals_T0_vs_T16, method = "BH")


# Combine results into a data frame, including gene names and log2FC calculations
results <- data.frame(
  Gene = rownames(mts.2.3.5.tub),
  T16_vs_T0_FDR = fdr_T0_vs_T16,
  T16_vs_T0_linearFC = linear_fc_T0_vs_T16,
  T16_vs_T0_log2FC = log2(abs(linear_fc_T0_vs_T16)),
  
  row.names = NULL
)

# Inspect first few rows of the results
head(results)
#> head(results)
#Gene T16_vs_T0_FDR T16_vs_T0_linearFC T16_vs_T0_log2FC
#    A2M  1.242387e-06        -73.4637844         6.198961
#  ACADM  1.352359e-07        -37.0455546         5.211229
#  ACAT1  3.146101e-08        -66.2191847         6.049177


# Set significance thresholds
fdr_threshold <- 0.05

# Add a "Significance" column based on FDR and log2FC thresholds
results$Significance <- with(results, 
                             (T16_vs_T0_FDR < fdr_threshold))

# View the updated results with the "Significance" column
head(results)
#> head(results)
#Gene T16_vs_T0_FDR T16_vs_T0_linearFC T16_vs_T0_log2FC Significance
#    A2M  1.242387e-06        -73.4637844         6.198961         TRUE
#  ACADM  1.352359e-07        -37.0455546         5.211229         TRUE
#  ACAT1  3.146101e-08        -66.2191847         6.049177         TRUE

# Save results to an Excel file
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)
write.xlsx(results, file = "Results_with_T16_vs_T0_Tub_Levels.xlsx")


ngoi <- results[results$Significance == TRUE, "Gene"]



#heatmap

tub.ngois<-mts.2.3.5.tub[ngoi,]

Annot_df <- data.frame(
  Timepoint = character(length = ncol(tub.ngois)),
  Treatment = character(length = ncol(tub.ngois)),
  DonorID = character(length = ncol(tub.ngois)),
  Location = character(length = ncol(tub.ngois)),
  Histology = character(length = ncol(tub.ngois))
)

Annot_df$Timepoint<-sub(".*_(T[0-9]+)_.*", "\\1",colnames(tub.ngois))
Annot_df$Treatment<-sub(".*_(.*?)_.*", "\\1", colnames(tub.ngois))
Annot_df$DonorID<-sub("^([A-Za-z0-9]+)_.*", "\\1",  colnames(tub.ngois))
Annot_df$Location <- sub("^[^_]+_[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(tub.ngois))
Annot_df$Histology <- sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(tub.ngois))



row.names(Annot_df) <- colnames(tub.ngois)


name.treatments <- function(timepoint, treatment) {
  if (timepoint == "T0" & treatment == "None") { 
    return("Naive")
  } else if (timepoint == "T16" & treatment == "None") {
    return("Vehicle")
  } else if (timepoint == "T16" & treatment == "IFN") {
    return("IFN")
  } else if (timepoint == "T16" & treatment == "LPS") {
    return("LPS")
  } else if (timepoint == "T16" & treatment == "PIC") {
    return("PIC")
  } else if (timepoint == "T16" & treatment == "Baricitinib") {
    return("Baricitinib")
  } else if (timepoint == "T16" & treatment == "Baricitinib+IFN") {
    return("Baricitinib+IFN")
  } else {
    return(NA)
  }
}



#Apply to data
Annot_df$Treatment <- mapply(name.treatments, Annot_df$Timepoint, Annot_df$Treatment)

desired_order <- c("Naive", "Vehicle", "IFN", "LPS", "PIC")




# Ensure that 'Treatment' column is a factor with the desired order
Annot_df$Treatment <- factor(Annot_df$Treatment, levels = desired_order)
Annot_df$Location <- factor(Annot_df$Location, levels = c("Edge","Central"))
Annot_df$DonorID <- factor(Annot_df$DonorID, levels = c("MTS002","MTS003","MTS005"))
# Sort Annot_df based on Treatment and Location
sorted_Annot_df <- Annot_df[order(Annot_df$Treatment, Annot_df$Location,Annot_df$DonorID), ]

# View the sorted annotated dataframe
head(sorted_Annot_df)
#head(sorted_Annot_df)
#                                           Timepoint Treatment DonorID Location          Histology
#MTS002_None_Tubulointerstitium_Edge_T0_031        T0     Naive  MTS002     Edge Tubulointerstitium
#MTS002_None_Tubulointerstitium_Edge_T0_032        T0     Naive  MTS002     Edge Tubulointerstitium
#MTS002_None_Tubulointerstitium_Edge_T0_040        T0     Naive  MTS002     Edge Tubulointerstitium


ordered.tub.data<-tub.ngois[,order(Annot_df$Treatment, Annot_df$Location)]

# Check if the columns match
identical(colnames(ordered.tub.data), rownames(sorted_Annot_df)) # Should return TRUE

#save it 
png("Tubules.T16vT0.heatmap.png", width = 1200, height = 1500, res = 150)

heatmap_plot <- pheatmap(ordered.tub.data,
                         color = my.colors,
                         cellwidth = 3,
                         cellheight = 0.02,
                         lwd = 0.5,
                         border_color = "darkgrey",
                         fontsize_col = 5,
                         fontsize_row = 5,
                         scale = "row",
                         breaks = my.breaks,
                         treeheight_row = 10,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         annotation_col = sorted_Annot_df,
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean")

# force the plot to render
print(heatmap_plot)

dev.off()






##########################################get up and down regulated genes ###################################

# Generate the heatmap object without drawing it yet
heatmap_obj <- pheatmap(ordered.tub.data, color = my.colors, cellwidth = 3,
                        cellheight = 0.02, lwd = 0.5,border_color = "darkgrey", fontsize_col = 5,
                        fontsize_row = 5, scale = "row", breaks = my.breaks, treeheight_row = 10,
                        show_rownames = FALSE, show_colnames = TRUE, annotation_col = sorted_Annot_df,
                        cluster_cols = FALSE, cluster_rows = TRUE, clustering_distance_rows = "euclidean",
                        silent = TRUE)  # <- do not plot yet

# Extract gene order from heatmap object
ordered_genes <- rownames(ordered.tub.data)[heatmap_obj$tree_row$order]

# Inspect or save
head(ordered_genes)
write.csv(ordered_genes, "Ordered_Genes_in_Tub.T0vT16.Heatmap.csv", row.names = FALSE)




# ---- START: Identify Upregulated Genes in T16 vs T0 ----

results.test<-results

# Add signed log2 fold change (preserve direction)
results.test$Signed_log2FC <- log2(abs(results$T16_vs_T0_linearFC)) * sign(results$T16_vs_T0_linearFC)

# Set FDR threshold
fdr_threshold <- 0.05

# Add Significance and Direction columns
results.test$Sig <- results.test$T16_vs_T0_FDR < fdr_threshold
results.test$Direction <- ifelse(results.test$Significance,
                                 ifelse(results.test$Signed_log2FC > 0, "Up in T16", "Down in T16"),
                                 "Not Significant")


# Subset: Upregulated genes in T16
upregulated_in_T16 <- subset(results.test, Direction == "Up in T16")
upregulated_in_T0<-subset(results.test, Direction == "Down in T16")


# View top results (optional)
head(upregulated_in_T16)

# Save to CSV
write.csv(upregulated_in_T16, "Upregulated_Tub_Genes-T16_T16_vs_T0.csv", row.names = FALSE)
write.csv(upregulated_in_T0, "Upregulated_Tub_Genes-T0_T16_vs_T0.csv", row.names = FALSE)
write.csv(results.test, "results.T16vT0.tubs.allgene+sig.csv", row.names = FALSE)    

##################################
#Volcanos 

##################################
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

res1<-results.test
res1$T16_vs_T0_pvals<-pvals_T0_vs_T16
res1$Signed_log2FC<-results.test$Signed_log2FC
 


# Set output file first
pdf(file = "MTS2,3,5-Quantile-T0-vs-T16-Volcano-Tubules.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_log2FC', # log2FC
                         y = 'T16_vs_T0_pvals', # p-value
                         title = 'T0 vs. T16 Tubules',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()

tub.res1<-res1[,c(1:6,8:9)]
tub.upndown.sig.genes<-ngoi
#write.csv(upregulated_in_T16, "Upregulated_GLom_Genes-T16_T16_vs_T0.csv", row.names = FALSE)


##################################################################################################################################

###################### Find shared up and down genes for Gloms and Tubes for pathway analysis 

############################################################################################################################


tub.res1
head(tub.res1)
#> head(tub.res1)
#Gene T16_vs_T0_FDR T16_vs_T0_linearFC T16_vs_T0_log2FC Significance Signed_log2FC       Direction T16_vs_T0_pvals
#    A2M  1.242387e-06        -73.4637844         6.198961         TRUE     -6.198961     Down in T16    5.918978e-08
#  ACADM  1.352359e-07        -37.0455546         5.211229         TRUE     -5.211229     Down in T16    3.579391e-09
#  ACAT1  3.146101e-08        -66.2191847         6.049177         TRUE     -6.049177     Down in T16    5.551350e-10

glom.res1
head(glom.res1)
#> head(glom.res1)
#Gene T16_vs_T0_FDR T16_vs_T0_linearFC T16_vs_T0_log2FC Significance Signed_log2FC       Direction T16_vs_T0_pvals
#    A2M  9.257025e-05        -358.213748         8.484677         TRUE     -8.484677     Down in T16    1.271096e-05
#  ACADM  4.605004e-02          -4.811989         2.266633         TRUE     -2.266633     Down in T16    2.307197e-02
#  ACAT1  2.395027e-05         -10.094643         3.335518         TRUE     -3.335518     Down in T16    2.528114e-06


# Merge both for easy comparison
merged.res1 <- merge(tub.res1, glom.res1, by = "Gene", suffixes = c(".tub", ".glom"))

# Shared upregulated: positive Signed_log2FC in both
Shared.up <- merged.res1[merged.res1$Direction.tub == "Up in T16" &
                           merged.res1$Direction.glom == "Up in T16", ]
dim(Shared.up) #should be 812 genes same as what I got in vendiagrams
#[1] 730  15 # this is new will have to make new vendiagram 

# Shared downregulated: negative Signed_log2FC in both
Shared.down <- merged.res1[merged.res1$Direction.tub =="Down in T16" &
                             merged.res1$Direction.glom == "Down in T16", ]

dim(Shared.down) #should be 938 genes same as what I got in vendiagrams
#[1] 800  15

Shared.sig.genes<-merged.res1[merged.res1$Significance.tub == "TRUE" &
                                merged.res1$Significance.glom == "TRUE", ]

Shared.up<-Shared.up[,"Gene"]
Shared.down<-Shared.down[,"Gene"]

Shared.sig.genes<-Shared.sig.genes[,"Gene"]
summary(Shared.sig.genes)
#Length     Class      Mode 
#1598 character character         # 1851 total shared genes may use this #heatmap? 

#total significant gene ( may be true for glom and tubs or individually)
glom.upndown.sig.genes<-glom.res1$Gene[glom.res1$Significance == TRUE]
total.sig.genes<- unique(c(glom.upndown.sig.genes,tub.upndown.sig.genes))
summary(total.sig.genes)
###################################### Pathway Analysis ##################################################
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("enrichplot", force = TRUE)
BiocManager::install("DOSE", force = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(magrittr)
library(dplyr)

Shared.up
Shared.down
bkgd.genes <- merged.res1[, "Gene"] 

######################### Shared gene list set#######################

up.gene_df <- bitr(Shared.up, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.gene_df <- bitr(Shared.down, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes_df<- bitr(bkgd.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


ego.up.Shared <- enrichGO( gene = up.gene_df$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 universe = bkgd.genes_df$ENTREZID,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 readable = TRUE)

ego.dn.Shared <- enrichGO(gene = dn.gene_df$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 universe = bkgd.genes_df$ENTREZID,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 readable = TRUE)



#Dot/box  plot for up  GO

png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 10, title = "GO Enrichment - Upregulated in T16")
dev.off()

#PDF of dot plot 
pdf("GO_Enrichment_Upregulated_T16.dotplot.pdf", width = 10, height = 7)
dotplot(ego.up.Shared, showCategory = 10, title = "GO Enrichment - Upregulated in T16")
dev.off()



#Dot/box  plot for down  GO

png("GO_Enrichment_Downregulated_T16.bar.png", width = 800, height = 600)
barplot(ego.dn.Shared, showCategory = 20, title = "GO Enrichment - Downregulated in T16")
dev.off()
png("GO_Enrichment_Downregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.dn.Shared, showCategory = 10, title = "GO Enrichment - Downregulated in T16")
dev.off()

#PDF of dot plot 
pdf("GO_Enrichment_Downregulated_T16.dotplot.pdf", width = 10, height = 7)
dotplot(ego.dn.Shared, showCategory = 10, title = "GO Enrichment - Downregulated in T16")
dev.off()












################################################################################################
#Heatmap combined 


################################################################################################

#ngois 
total.sig.genes

#data
mts.2.3.5.glom$Gene <- rownames(mts.2.3.5.glom)
mts.2.3.5.tub$Gene <- rownames(mts.2.3.5.tub)


# Merge both for easy comparison
all.T16vT0.Glom.Tub.compare <- merge(mts.2.3.5.glom, mts.2.3.5.tub, by = "Gene", all=TRUE)

rownames(all.T16vT0.Glom.Tub.compare) <- all.T16vT0.Glom.Tub.compare$Gene
all.T16vT0.Glom.Tub.compare$Gene <- NULL


All.ngois<-all.T16vT0.Glom.Tub.compare[total.sig.genes,]


Annot_df <- data.frame(
  Timepoint = character(length = ncol(All.ngois)),
  DonorID = character(length = ncol(All.ngois)),
  Location = character(length = ncol(All.ngois)),
  Histology = character(length = ncol(All.ngois))
)

Annot_df$Timepoint<-sub(".*_(T[0-9]+)_.*", "\\1",colnames(All.ngois))
Annot_df$DonorID<-sub("^([A-Za-z0-9]+)_.*", "\\1",  colnames(All.ngois))
Annot_df$Location <- sub("^[^_]+_[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(All.ngois))
Annot_df$Histology <- sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(All.ngois))



row.names(Annot_df) <- colnames(All.ngois)





# Ensure that 'Treatment' column is a factor with the desired order

Annot_df$Location <- factor(Annot_df$Location, levels = c("Edge","Central"))
Annot_df$DonorID <- factor(Annot_df$DonorID, levels = c("MTS002","MTS003","MTS005"))
Annot_df$Histology <- factor(Annot_df$Histology, levels = c("Tubulointerstitium","Glomerulus"))
Annot_df$Timepoint<- factor(Annot_df$Timepoint, levels = c("T0","T16"))
# Sort Annot_df based on Treatment and Location
sorted_Annot_df <- Annot_df[order(Annot_df$Timepoint,Annot_df$Histology, Annot_df$Location,Annot_df$DonorID), ]

# View the sorted annotated dataframe
head(sorted_Annot_df)
# head(sorted_Annot_df)
#Timepoint Treatment DonorID Location          Histology
#MTS002_None_Tubulointerstitium_Edge_T0_031        T0     Naive  MTS002     Edge Tubulointerstitium
#MTS002_None_Tubulointerstitium_Edge_T0_032        T0     Naive  MTS002     Edge Tubulointerstitium
#MTS002_None_Tubulointerstitium_Edge_T0_040        T0     Naive  MTS002     Edge Tubulointerstitium
#MTS002_None_Tubulointerstitium_Edge_T0_041        T0     Naive  MTS002     Edge Tubulointerstitium
#MTS002_None_Glomerulus_Edge_T0_001                T0     Naive  MTS002     Edge         Glomerulus


ordered.all.data<-All.ngois[,order(Annot_df$Timepoint,Annot_df$Histology, Annot_df$Location,Annot_df$DonorID)]

# Check if the columns match
identical(colnames(ordered.all.data), rownames(sorted_Annot_df)) # Should return TRUE

#save it 
png("Glom.n.Tub.T16vT0.heatmap.png", width = 1200, height = 1500, res = 150)

heatmap_plot <- pheatmap(ordered.all.data,
                         color = my.colors,
                         cellwidth = 3,
                         cellheight = 0.02,
                         lwd = 0.5,
                         border_color = "darkgrey",
                         fontsize_col = 5,
                         fontsize_row = 5,
                         scale = "row",
                         breaks = my.breaks,
                         treeheight_row = 10,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         annotation_col = sorted_Annot_df,
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean")

# force the plot to render
print(heatmap_plot)

dev.off()




# Save as PDF
pdf("Glom.n.Tub.T16vT0.heatmap.pdf", width = 10, height = 12)  # Adjust as needed

# Generate the heatmap
heatmap_plot <- pheatmap(ordered.all.data,
                         color = my.colors,
                         cellwidth = 3,
                         cellheight = 0.02,
                         lwd = 0.5,
                         border_color = "darkgrey",
                         fontsize_col = 5,
                         fontsize_row = 5,
                         scale = "row",
                         breaks = my.breaks,
                         treeheight_row = 10,
                         show_rownames = FALSE,
                         show_colnames = TRUE,
                         annotation_col = sorted_Annot_df,
                         cluster_cols = FALSE,
                         cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean")

# Force the plot to render
print(heatmap_plot)

# Close the device
dev.off()

################################################################################################
 # Up and down processes for tubules and gloms separately ( unique to tubules and gloms)

BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("enrichplot", force = TRUE)
BiocManager::install("DOSE", force = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(magrittr)
library(dplyr)

bkgd.genes <- merged.res1[, "Gene"] 

######################### generate unique gene list set#######################

#editing here 

# unique upregulated tubules: 
Unique.tub.up <- merged.res1[merged.res1$Direction.tub == "Up in T16" &
                           merged.res1$Direction.glom == "Down in T16"|
                           merged.res1$Direction.tub == "Up in T16" &
                           merged.res1$Direction.glom == "Not Significant", ]
dim(Unique.tub.up) #407  #correct checked ven 
Unique.tub.up<-Unique.tub.up[,"Gene"]

# unique Downregulated tubules: 
Unique.tub.dn <- merged.res1[merged.res1$Direction.tub == "Down in T16" &
                               merged.res1$Direction.glom == "Up in T16"|
                               merged.res1$Direction.tub == "Down in T16" &
                               merged.res1$Direction.glom == "Not Significant", ]
dim(Unique.tub.dn) # 550  #Correct!!
Unique.tub.dn<-Unique.tub.dn[,"Gene"]

# unique upregulated Glomeruli: 
Unique.glom.up <- merged.res1[merged.res1$Direction.tub == "Down in T16" &
                               merged.res1$Direction.glom == "Up in T16"|
                               merged.res1$Direction.tub == "Not Significant" &
                               merged.res1$Direction.glom == "Up in T16", ]
dim(Unique.glom.up) #744  #correct checked ven
Unique.glom.up<-Unique.glom.up[,"Gene"]

# unique Downregulated Glomeruli: 
Unique.glom.dn  <- merged.res1[merged.res1$Direction.tub == "Up in T16" &
                               merged.res1$Direction.glom == "Down in T16"|
                               merged.res1$Direction.tub == "Not Significant" &
                               merged.res1$Direction.glom == "Down in T16", ]
dim(Unique.glom.dn) # 905 #correct!
Unique.glom.dn<-Unique.glom.dn[,"Gene"]



Unique.tub.up_df <- bitr(Unique.tub.up, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Unique.tub.dn_df <- bitr(Unique.tub.dn, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Unique.glom.up_df <- bitr(Unique.glom.up, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
Unique.glom.dn_df <- bitr(Unique.glom.dn, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)



bkgd.genes_df<- bitr(bkgd.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)





ego.tub.up.genes<- enrichGO( gene = Unique.tub.up_df$ENTREZID,
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           universe = bkgd.genes_df$ENTREZID,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2,
                           readable = TRUE)

ego.tub.dn.genes <- enrichGO(gene = Unique.tub.dn_df$ENTREZID,
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          universe = bkgd.genes_df$ENTREZID,
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          readable = TRUE)

ego.glom.up.genes<- enrichGO( gene = Unique.glom.up_df$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             universe = bkgd.genes_df$ENTREZID,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2,
                             readable = TRUE)

ego.glom.dn.genes <- enrichGO(gene = Unique.glom.dn_df$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             universe = bkgd.genes_df$ENTREZID,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2,
                             readable = TRUE)


#Dot plot for up/down GO T16 glom and Tubs seperately  (unique gene list )



#PDF of dot plot 
pdf("GO_UP_T16_tubes.dotplot.pdf", width = 10, height = 7)
dotplot(ego.tub.up.genes, showCategory = 10, title = "GO Enrichment - Upregulated in T16 Tubules")
dev.off()

#PDF of dot plot 
pdf("GO_DOWN_T16_tubes.dotplot.pdf", width = 10, height = 7)
dotplot(ego.tub.dn.genes, showCategory = 10, title = "GO Enrichment - Downregulated in T16 Tubules")
dev.off()

#PDF of dot plot 
pdf("GO_UP_T16_gloms.dotplot.pdf", width = 10, height = 7)
dotplot(ego.glom.up.genes, showCategory = 10, title = "GO Enrichment - Upregulated in T16 Glomeruli")
dev.off()

#PDF of dot plot 
pdf("GO_DOWN_T16_gloms.dotplot.pdf", width = 10, height = 7)
dotplot(ego.glom.dn.genes, showCategory = 10, title = "GO Enrichment - Downregulated in T16 Glomeruli")
dev.off()



merged.res1

write.csv(merged.res1, "Timepoint.LMM.merged.TubnGlom.csv", row.names = FALSE)









