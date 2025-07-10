#LMM for Glom MTS 2,3,5 

# Load required package (stats is available by default in R)
library(stats)

# Read the data file

mts.2.3.5.glom <- read.csv("~/Desktop/summer25/mts.2.3.5.glom.quantilenorm.data.csv",
                           row.names = 1)
head(mts.2.3.5.glom)

timepoint_factor <- sub(".*_T(\\d+)_.*", "T\\1", colnames(mts.2.3.5.glom))

keep_columns <- timepoint_factor != "T0"

mts.2.3.5.glom_noT0 <- mts.2.3.5.glom[, keep_columns]

data<-mts.2.3.5.glom_noT0
########################## Comparisions on treatment #####################################


# Create treatment factor from column names
treatment_factor <- factor(sub(".*_(None|IFN|LPS|PIC)_.*", "\\1", colnames(mts.2.3.5.glom_noT0)),
                           levels = c("None", "IFN", "LPS", "PIC"))

# Set "None" as the reference level
treatment_factor <- relevel(treatment_factor, ref = "None")  # Reference level set to "None"

print(levels(treatment_factor))  # Verify levels are assigned correctly
# Expected output: [1] "None" "IFN" "LPS" "PIC"

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
contrast_ifn_vs_none <- c(0, 1, 0, 0)       # IFN vs. None
contrast_lps_vs_none <- c(0, 0, 1, 0)       # LPS vs. None
contrast_pic_vs_none <- c(0, 0, 0, 1)       # PIC vs. None

contrast_lps_vs_ifn <- c(0, -1, 1, 0)       # LPS vs. IFN
contrast_pic_vs_ifn <- c(0, -1, 0, 1)       # PIC vs. IFN
contrast_pic_vs_lps <- c(0, 0, -1, 1)       # PIC vs. LPS

# For each gene, calculate linear fold change and p-values for each contrast
results_ifn_vs_none <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, treatment_factor, contrast_ifn_vs_none))
results_lps_vs_none <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, treatment_factor, contrast_lps_vs_none))
results_pic_vs_none <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, treatment_factor, contrast_pic_vs_none))

results_lps_vs_ifn <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, treatment_factor, contrast_lps_vs_ifn))
results_pic_vs_ifn <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, treatment_factor, contrast_pic_vs_ifn))
results_pic_vs_lps <- apply(data, 1, function(x) perform_contrast_and_linear_fc(x, treatment_factor, contrast_pic_vs_lps))

# Extract p-values and linear fold changes for each contrast
pvals_ifn_vs_none <- sapply(results_ifn_vs_none, function(x) x$p_value)
linear_fc_ifn_vs_none <- sapply(results_ifn_vs_none, function(x) x$linear_fc)

pvals_lps_vs_none <- sapply(results_lps_vs_none, function(x) x$p_value)
linear_fc_lps_vs_none <- sapply(results_lps_vs_none, function(x) x$linear_fc)

pvals_pic_vs_none <- sapply(results_pic_vs_none, function(x) x$p_value)
linear_fc_pic_vs_none <- sapply(results_pic_vs_none, function(x) x$linear_fc)

pvals_lps_vs_ifn <- sapply(results_lps_vs_ifn, function(x) x$p_value)
linear_fc_lps_vs_ifn <- sapply(results_lps_vs_ifn, function(x) x$linear_fc)

pvals_pic_vs_ifn <- sapply(results_pic_vs_ifn, function(x) x$p_value)
linear_fc_pic_vs_ifn <- sapply(results_pic_vs_ifn, function(x) x$linear_fc)

pvals_pic_vs_lps <- sapply(results_pic_vs_lps, function(x) x$p_value)
linear_fc_pic_vs_lps <- sapply(results_pic_vs_lps, function(x) x$linear_fc)

# Apply FDR correction (Benjamini–Hochberg) for each set of p-values
fdr_ifn_vs_none <- p.adjust(pvals_ifn_vs_none, method = "BH")
fdr_lps_vs_none <- p.adjust(pvals_lps_vs_none, method = "BH")
fdr_pic_vs_none <- p.adjust(pvals_pic_vs_none, method = "BH")
fdr_lps_vs_ifn <- p.adjust(pvals_lps_vs_ifn, method = "BH")
fdr_pic_vs_ifn <- p.adjust(pvals_pic_vs_ifn, method = "BH")
fdr_pic_vs_lps <- p.adjust(pvals_pic_vs_lps, method = "BH")

# Combine results into a data frame, including gene names and log2FC calculations
results <- data.frame(
  Gene = rownames(data),
  IFN_vs_None_FDR = fdr_ifn_vs_none,
  IFN_vs_None_linearFC = linear_fc_ifn_vs_none,
  IFN_vs_None_log2FC = log2(abs(linear_fc_ifn_vs_none)),
  
  LPS_vs_None_FDR = fdr_lps_vs_none,
  LPS_vs_None_linearFC = linear_fc_lps_vs_none,
  LPS_vs_None_log2FC = log2(abs(linear_fc_lps_vs_none)),
  
  PIC_vs_None_FDR = fdr_pic_vs_none,
  PIC_vs_None_linearFC = linear_fc_pic_vs_none,
  PIC_vs_None_log2FC = log2(abs(linear_fc_pic_vs_none)),
  
  LPS_vs_IFN_FDR = fdr_lps_vs_ifn,
  LPS_vs_IFN_linearFC = linear_fc_lps_vs_ifn,
  LPS_vs_IFN_log2FC = log2(abs(linear_fc_lps_vs_ifn)),
  
  PIC_vs_IFN_FDR = fdr_pic_vs_ifn,
  PIC_vs_IFN_linearFC = linear_fc_pic_vs_ifn,
  PIC_vs_IFN_log2FC = log2(abs(linear_fc_pic_vs_ifn)),
  
  PIC_vs_LPS_FDR = fdr_pic_vs_lps,
  PIC_vs_LPS_linearFC = linear_fc_pic_vs_lps,
  PIC_vs_LPS_log2FC = log2(abs(linear_fc_pic_vs_lps)),
  
  row.names = NULL
)

head(results)

# Set significance thresholds
fdr_threshold <- 0.05


# Add significance columns
results$IFN_vs_None_Significance <- with(results, (IFN_vs_None_FDR < fdr_threshold))
results$LPS_vs_None_Significance <- with(results, (LPS_vs_None_FDR < fdr_threshold))
results$PIC_vs_None_Significance <- with(results, (PIC_vs_None_FDR < fdr_threshold))
results$LPS_vs_IFN_Significance <- with(results, (LPS_vs_IFN_FDR < fdr_threshold))
results$PIC_vs_IFN_Significance <- with(results, (PIC_vs_IFN_FDR < fdr_threshold))
results$PIC_vs_LPS_Significance <- with(results, (PIC_vs_LPS_FDR < fdr_threshold))


# Inspect first few rows of the results
head(results)
#> head(results)
#Gene IFN_vs_None_FDR IFN_vs_None_linearFC IFN_vs_None_log2FC LPS_vs_None_FDR LPS_vs_None_linearFC LPS_vs_None_log2FC PIC_vs_None_FDR
#    A2M      0.18508086           16.1785120          4.0160070      0.66732698           -5.8781824          2.5553701       0.9988350
#  ACADM      0.82598020           -0.7432661         -0.4280494      0.78406317           -0.8072909         -0.3088394       0.9994678
#  ACAT1      0.83133541           -0.7444678         -0.4257187      0.68738464           -1.1604467          0.2146803       0.9988350



# Save results to an Excel file
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)
write.xlsx(results, file = "Results_LMM_Treatment_Glom.xlsx")



IFN_vs_None.genes<-results[results$IFN_vs_None_Significance == TRUE, "Gene"]
#513
LPS_vs_None.genes<-results[results$LPS_vs_None_Significance == TRUE, "Gene"]
#598
PIC_vs_None.genes<-results[results$PIC_vs_None_Significance == TRUE, "Gene"]
#0 <- strange
LPS_vs_IFN.genes<-results[results$LPS_vs_IFN_Significance == TRUE, "Gene"]
#453
PIC_vs_IFN.genes<-results[results$PIC_vs_IFN_Significance == TRUE, "Gene"]
#429
PIC_vs_LPS.genes<-results[results$PIC_vs_LPS_Significance == TRUE, "Gene"]
#418


goi<- unique(c(IFN_vs_None.genes,LPS_vs_None.genes,PIC_vs_None.genes,
               LPS_vs_IFN.genes,PIC_vs_IFN.genes,PIC_vs_LPS.genes))

summary(goi)

ngoi <- results[results$IFN_vs_None_Significance == FALSE &
                  results$LPS_vs_None_Significance == FALSE &
                  results$PIC_vs_None_Significance == FALSE &
                  results$LPS_vs_IFN_Significance == FALSE &
                  results$PIC_vs_IFN_Significance == FALSE &
                  results$PIC_vs_LPS_Significance == FALSE,
                "Gene"]


summary(ngoi)

tdata <- t(data)
merged_df <- merge(tdata, Annot_df, by=0)
head(merged_df)
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL

  

#Length     Class      Mode 
#1266 character character              # 1266 (total sig treatment genes for gloms)

#####################  heatmap############################################

glom.gois<-mts.2.3.5.glom_noT0[goi,]

Annot_df <- data.frame(
  Treatment = character(length = ncol(glom.gois)),
  DonorID = character(length = ncol(glom.gois)),
  GlomLocation = character(length = ncol(glom.gois)),
  Histology = character(length = ncol(glom.gois))
)

Annot_df$Treatment<-sub(".*_(.*?)_.*", "\\1", colnames(glom.gois))
Annot_df$DonorID<-sub("^([A-Za-z0-9]+)_.*", "\\1",  colnames(glom.gois))
Annot_df$GlomLocation <- sub("^[^_]+_[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(glom.gois))
Annot_df$Histology <- sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(glom.gois))
Annot_df$Time <- sub("^[^_]+_[^_]+_([^_]+)_.*", "\\1", colnames(glom.gois))
Annot_df$Time <- sub("^(?:[^_]*_){4}([^_]+).*", "\\1", colnames(glom.gois))


row.names(Annot_df) <- colnames(glom.gois)



#apply to data 
desired_order <- c("None", "IFN", "LPS", "PIC")



# Ensure that 'Treatment' column is a factor with the desired order
Annot_df$Treatment <- factor(Annot_df$Treatment, levels = desired_order)
Annot_df$GlomLocation <- factor(Annot_df$GlomLocation, levels = c("Edge","Central"))
Annot_df$DonorID <- factor(Annot_df$DonorID, levels = c("MTS002","MTS003","MTS005"))
# Sort Annot_df based on Treatment and GlomLocation
sorted_Annot_df <- Annot_df[order(Annot_df$Treatment, Annot_df$GlomLocation,Annot_df$DonorID), ]

# View the sorted annotated dataframe
head(sorted_Annot_df)



ordered.glom.data<-glom.gois[,order(Annot_df$Treatment, Annot_df$GlomLocation,Annot_df$DonorID)]

# Check if the columns match
identical(colnames(ordered.glom.data), rownames(sorted_Annot_df)) # Should return TRUE

my.breaks <- c(seq(-2, 2, by=0.1)) 

my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("white", "orange", "red"))(length(my.breaks)/2))


library(pheatmap)
pheatmap(ordered.glom.data, color=my.colors, cellwidth=4, 
         cellheight=0.3, lwd=0.5, border_color="darkgrey", fontsize_col=5, fontsize_row=5, scale="row", 
         breaks= my.breaks, treeheight_row=10,  show_rownames = FALSE, show_colnames = TRUE, 
         annotation_col = sorted_Annot_df,cluster_cols=TRUE, cluster_rows=TRUE, clustering_distance_rows="euclidean")

#save it 
png("Glom.treatment.heatmap.png", width = 1200, height = 1500, res = 150)

heatmap_plot <- pheatmap(ordered.glom.data,
                         color = my.colors,
                         cellwidth = 4,
                         cellheight = 0.3,
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
pdf("Glom.treatment.heatmap.pdf", width = 10, height = 12)  # Adjust dimensions if needed

# Generate the heatmap
heatmap_plot <- pheatmap(ordered.glom.data,
                         color = my.colors,
                         cellwidth = 4,
                         cellheight = 0.3,
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

# Render the plot
print(heatmap_plot)

# Close the PDF device
dev.off()






# Extract gene order from heatmap object
ordered_genes <- rownames(ordered.glom.data)[heatmap_plot$tree_row$order]

# Inspect or save
head(ordered_genes)
write.csv(ordered_genes, "Ordered_Genes_in_Glom.treatment.Heatmap.csv", row.names = FALSE)

##########################################################################################################################

#                                ---- START: Identify Up/downregulated Genes in treatments  ----

###########################################################################################################################
results.test<-results


# Add signed log2 fold change (preserve direction)
results.test$Signed_IFN_vs_None_log2FC <- log2(abs(results$IFN_vs_None_linearFC)) * sign(results$IFN_vs_None_linearFC) #IFN_vs_None
results.test$Signed_LPS_vs_None_log2FC <- log2(abs(results$LPS_vs_None_linearFC)) * sign(results$LPS_vs_None_linearFC) #LPS_vs_None
results.test$Signed_PIC_vs_None_log2FC <- log2(abs(results$PIC_vs_None_linearFC)) * sign(results$PIC_vs_None_linearFC) #PIC_vs_None
results.test$Signed_LPS_vs_IFN_log2FC <- log2(abs(results$LPS_vs_IFN_linearFC)) * sign(results$LPS_vs_IFN_linearFC) #LPS_vs_IFN
results.test$Signed_PIC_vs_IFN_log2FC <- log2(abs(results$PIC_vs_IFN_linearFC)) * sign(results$PIC_vs_IFN_linearFC) #PIC_vs_IFN
results.test$Signed_PIC_vs_LPS_log2FC <- log2(abs(results$PIC_vs_LPS_linearFC)) * sign(results$PIC_vs_LPS_linearFC) #PIC_vs_LPS


# Add Direction columns

#IFN_vs_None
results.test$IFN_vs_None_Direction <- ifelse(results.test$IFN_vs_None_Significance,
                                             ifelse(results.test$Signed_IFN_vs_None_log2FC > 0, "Up in IFN","Down in IFN"),
                                             "Not Significant")

#LPS_vs_None
results.test$LPS_vs_None_Direction <- ifelse(results.test$LPS_vs_None_Significance,
                                             ifelse(results.test$Signed_LPS_vs_None_log2FC > 0, "Up in LPS","Down in LPS"),
                                             "Not Significant")
#PIC_vs_None
results.test$PIC_vs_None_Direction <- ifelse(results.test$PIC_vs_None_Significance,
                                             ifelse(results.test$Signed_PIC_vs_None_log2FC> 0, "Up in PIC","Down in PIC"),
                                             "Not Significant")

#LPS_vs_IFN
results.test$LPS_vs_IFN_Direction <- ifelse(results.test$LPS_vs_IFN_Significance,
                                            ifelse(results.test$Signed_LPS_vs_IFN_log2FC > 0, "Up in LPS","Up in IFN"),
                                            "Not Significant")

#PIC_vs_IFN
results.test$PIC_vs_IFN_Direction <- ifelse(results.test$PIC_vs_IFN_Significance,
                                            ifelse(results.test$Signed_PIC_vs_IFN_log2FC > 0, "Up in PIC", "Up in IFN"),
                                            "Not Significant")
#PIC_vs_LPS
results.test$PIC_vs_LPS_Direction <- ifelse(results.test$PIC_vs_LPS_Significance,
                                            ifelse(results.test$Signed_PIC_vs_LPS_log2FC > 0, "Up in PIC", "Up in LPS"),
                                            "Not Significant")

### Subset: Upregulated genes by treatment 

# treatment vs. control
up.IFN <- subset(results.test, IFN_vs_None_Direction == "Up in IFN")
down.IFN <- subset(results.test, IFN_vs_None_Direction == "Down in IFN")
up.LPS <- subset(results.test, LPS_vs_None_Direction == "Up in LPS")
down.LPS <- subset(results.test, LPS_vs_None_Direction == "Down in LPS")
up.PIC <- subset(results.test, PIC_vs_None_Direction == "Up in PIC")
down.PIC <- subset(results.test, PIC_vs_None_Direction == "Down in PIC")

#treatment vs treatment
up.LPS.LPSvIFN <- subset(results.test, LPS_vs_IFN_Direction == "Up in LPS")
up.IFN.LPSvIFN <- subset(results.test, LPS_vs_IFN_Direction == "Up in IFN")
up.PIC.PICvIFN <- subset(results.test, PIC_vs_IFN_Direction == "Up in PIC")
up.IFN.PICvIFN <- subset(results.test, PIC_vs_IFN_Direction == "Up in IFN")
up.PIC.PICvLPS <- subset(results.test, PIC_vs_LPS_Direction == "Up in PIC")
up.LPS.PICvLPS <- subset(results.test, PIC_vs_LPS_Direction == "Up in LPS")



# Save to CSV
write.csv(up.IFN, "Up.IFN.subset.IFNvNone.csv", row.names = FALSE)
write.csv(down.IFN, "down.IFN.subset.IFNvNone.csv", row.names = FALSE)
write.csv(up.LPS, "Up.LPS.subset.LPSvNone.csv", row.names = FALSE)
write.csv(down.LPS, "down.LPS.subset.LPSvNone.csv", row.names = FALSE)
#write.csv(up.PIC, "Up.IFN.subset.PICvNone.csv", row.names = FALSE)       # didnt save because no significant genes (PIC v None)
#write.csv(down.PIC, "down.IFN.subset.PICvNone.csv", row.names = FALSE)   # didnt save because no significant genes (PIC v None)
write.csv(up.LPS.LPSvIFN, "Up.LPS.subset.LPSvIFN.csv", row.names = FALSE)
write.csv(up.IFN.LPSvIFN, "Up.IFN.subset.LPSvIFN.csv", row.names = FALSE)
write.csv(up.PIC.PICvIFN, "Up.PIC.subset.PICvIFN.csv", row.names = FALSE)
write.csv(up.IFN.PICvIFN, "Up.IFN.subset.PICvIFN.csv", row.names = FALSE)
write.csv(up.PIC.PICvLPS, "Up.PIC.subset.PICvLPS.csv", row.names = FALSE)
write.csv(up.LPS.PICvLPS, "Up.LPS.subset.PICvLPS.csv", row.names = FALSE)
write.csv(results.test, "results.test.csv", row.names = FALSE)
##################################
#Volcanos 

##################################
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)


#IFN_vs_None

res1<-results.test

#add p-values to res1
res1$pvals_ifn_vs_none<-pvals_ifn_vs_none #IFN_vs_None
res1$pvals_lps_vs_none<-pvals_lps_vs_none #LPS_vs_None
res1$pvals_pic_vs_none<-pvals_pic_vs_none #PIC_vs_None
res1$pvals_lps_vs_ifn<-pvals_lps_vs_ifn #LPS_vs_IFN
res1$pvals_pic_vs_ifn<-pvals_pic_vs_ifn #PIC_vs_IFN
res1$pvals_pic_vs_lps<-pvals_pic_vs_lps #PIC_vs_LPS


# Set output file first
pdf(file = "MTS2,3,5-Quantile-IFN-vs-None-Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_IFN_vs_None_log2FC', # log2FC
                         y = 'pvals_ifn_vs_none', # p-value
                         title = 'Control vs. IFN Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()


#LPS_vs_None



# Set output file first
pdf(file = "MTS2,3,5-Quantile-LPS_vs_None-Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_LPS_vs_None_log2FC', # log2FC
                         y = 'pvals_lps_vs_none', # p-value
                         title = 'Control vs. LPS Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()


#PIC_vs_None



# Set output file first
pdf(file = "MTS2,3,5-Quantile-PIC_vs_None-Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_PIC_vs_None_log2FC', # log2FC
                         y = 'pvals_pic_vs_none', # p-value
                         title = 'Control vs. PIC Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()



#LPS_vs_IFN


# Set output file first
pdf(file = "MTS2,3,5-Quantile-LPS_vs_IFN-Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_LPS_vs_IFN_log2FC', # log2FC
                         y = 'pvals_lps_vs_ifn', # p-value
                         title = 'IFN vs. LPS Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()



#PIC_vs_IFN


# Set output file first
pdf(file = "MTS2,3,5-Quantile-PIC_vs_IFN-Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_PIC_vs_IFN_log2FC', # log2FC
                         y = 'pvals_pic_vs_ifn', # p-value
                         title = 'IFN vs. PIC Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()



#PIC_vs_LPS

# Set output file first
pdf(file = "MTS2,3,5-Quantile-PIC_vs_LPS-Volcano-Glom.pdf", width = 8, height = 7)

# Create and plot the volcano
plot1 <- EnhancedVolcano(res1,
                         lab = res1$Gene,
                         x = 'Signed_PIC_vs_LPS_log2FC', # log2FC
                         y = 'pvals_pic_vs_lps', # p-value
                         title = 'LPS vs. PIC Glomeruli',
                         pCutoff = 0.05,
                         FCcutoff = 0.584962501,
                         col = c('grey', 'grey', 'grey', 'red3'),
                         pointSize = 2,
                         labSize = 4)

print(plot1)  # Print to the PDF device

# Close the PDF device
dev.off()

##############################################################################################################

#######                      Pathway analysis 

#############################################################################################################

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

# treatment vs. control
up.IFN.only.genes <-  up.IFN[, "Gene"]
down.IFN.only.genes <-  down.IFN[, "Gene"]
up.LPS.only.genes <-  up.LPS[, "Gene"]
down.LPS.only.genes <-  down.LPS[, "Gene"]
#up.PIC.only.genes <-  up.PIC[, "Gene"]  # no pic up sig genes 
#down.PIC.only.genes <-  down.PIC[, "Gene"] # no pic up sig genes

#treatment vs treatment
up.LPS.LPSvIFN.genes <-  up.LPS.LPSvIFN[, "Gene"]
up.IFN.LPSvIFN.genes <-  up.IFN.LPSvIFN[, "Gene"]
up.PIC.PICvIFN.genes <-  up.PIC.PICvIFN[, "Gene"]
up.IFN.PICvIFN.genes <-  up.IFN.PICvIFN[, "Gene"]
up.PIC.PICvLPS.genes <-  up.PIC.PICvLPS[, "Gene"]
up.LPS.PICvLPS.genes <-  up.LPS.PICvLPS[, "Gene"]

bkgd.genes <- res1[, "Gene"]

#create gene ID dataframes

up.IFN.only.genes_df <- bitr(up.IFN.only.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
down.IFN.only.genes_df <- bitr(down.IFN.only.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.LPS.only.genes_df<- bitr(up.LPS.only.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
down.LPS.only.genes_df<- bitr(down.LPS.only.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

up.LPS.LPSvIFN.genes_df <- bitr(up.LPS.LPSvIFN.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.IFN.LPSvIFN.genes_df <- bitr(up.IFN.LPSvIFN.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.PIC.PICvIFN.genes_df <- bitr(up.PIC.PICvIFN.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.IFN.PICvIFN.genes_df <- bitr(up.IFN.PICvIFN.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.PIC.PICvLPS.genes_df <- bitr(up.PIC.PICvLPS.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
up.LPS.PICvLPS.genes_df <- bitr(up.LPS.PICvLPS.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

bkgd.genes_df<- bitr(bkgd.genes, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


#create enrichGo results 



#upregulated LPS/IFN vs NONE
ego.up.IFN.only <- enrichGO( gene = up.IFN.only.genes_df$ENTREZID,
                             OrgDb = org.Hs.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             universe = bkgd.genes_df$ENTREZID,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2,
                             readable = TRUE)


ego.up.LPS.only <- enrichGO(gene = up.LPS.only.genes_df$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            ont = "BP",
                            pAdjustMethod = "BH",
                            universe = bkgd.genes_df$ENTREZID,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2,
                            readable = TRUE)

#Downregulated LPS/IFN vs NONE
ego.down.IFN.only <- enrichGO( gene = down.IFN.only.genes_df$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               universe = bkgd.genes_df$ENTREZID,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               readable = TRUE)


ego.down.LPS.only <- enrichGO(gene = down.LPS.only.genes_df$ENTREZID,
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              pAdjustMethod = "BH",
                              universe = bkgd.genes_df$ENTREZID,
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2,
                              readable = TRUE)








#treatments 

ego.up.LPS.LPSvIFN <- enrichGO( gene = up.LPS.LPSvIFN.genes_df$ENTREZID,
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pAdjustMethod = "BH",
                                universe = bkgd.genes_df$ENTREZID,
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2,
                                readable = TRUE)

ego.up.IFN.LPSvIFN <- enrichGO(gene = up.IFN.LPSvIFN.genes_df$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               universe = bkgd.genes_df$ENTREZID,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               readable = TRUE)

ego.up.PIC.PICvIFN <- enrichGO(gene = up.PIC.PICvIFN.genes_df$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               universe = bkgd.genes_df$ENTREZID,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               readable = TRUE)


ego.up.IFN.PICvIFN <- enrichGO(gene = up.IFN.PICvIFN.genes_df$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               universe = bkgd.genes_df$ENTREZID,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               readable = TRUE)
ego.up.PIC.PICvLPS <- enrichGO(gene = up.PIC.PICvLPS.genes_df$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               universe = bkgd.genes_df$ENTREZID,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               readable = TRUE)
ego.up.LPS.PICvLPS <- enrichGO(gene = up.LPS.PICvLPS.genes_df$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               universe = bkgd.genes_df$ENTREZID,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2,
                               readable = TRUE)

# editing here plotting 


#Dot/box  plot for up  GO

#ego.up.IFN.only
png("GO_Enrichment_Up_IFN.only.png", width = 800, height = 600)
barplot(ego.up.IFN.only, showCategory = 20, title = "GO Enrichment - Upregulated in IFN IFN.v.None")
dev.off()
png("GO_Enrichment_Up_IFN.only.dot.png", width = 800, height = 600)
dotplot(ego.up.IFN.only, showCategory = 20, title = "GO Enrichment - Upregulated in IFN IFN.v.None")
dev.off()

# Save as PDF
pdf("GO_Enrichment_Up_IFN.only.dotplot.pdf", width = 10, height = 7)
dotplot(ego.up.IFN.only, showCategory = 10, title = "GO Enrichment - Upregulated in IFN IFN.v.None")
dev.off()


#ego.down.IFN.only
png("GO_Enrichment_Down_IFN.only.png", width = 800, height = 600)
barplot(ego.down.IFN.only, showCategory = 20, title = "GO Enrichment - Downregulated in IFN IFN.v.None")
dev.off()
png("GO_Enrichment_Down_IFN.only.dot.png", width = 800, height = 600)
dotplot(ego.down.IFN.only, showCategory = 20, title = "GO Enrichment - Downregulated in IFN IFN.v.None")
dev.off()

# Save as PDF
pdf("GO_Enrichment_Down_IFN.only.dotplot.pdf", width = 10, height = 7)
dotplot(ego.down.IFN.only, showCategory = 10, title = "GO Enrichment - Downregulated in IFN IFN.v.None")
dev.off()



#ego.up.LPS.only
png("GO_Enrichment_up.LPS.only.png", width = 800, height = 600)
barplot(ego.up.LPS.only, showCategory = 20, title = "GO Enrichment - Upregulated  LPS LPS.v.None")
dev.off()
png("GO_Enrichment_up.LPS.only.dot.png", width = 800, height = 600)
dotplot(ego.up.LPS.only, showCategory = 20, title = "GO Enrichment - Upregulated  LPS LPS.v.None")
dev.off()

# Save as PDF
pdf("GO_Enrichment_up.LPS.only.dotplot.pdf", width = 10, height = 7)
dotplot(ego.up.LPS.only, showCategory = 10, title = "GO Enrichment - Upregulated LPS LPS.v.None")
dev.off()



#ego.down.LPS.only
png("GO_Enrichment_down.LPS.only.png", width = 800, height = 600)
barplot(ego.down.LPS.only, showCategory = 20, title = "GO Enrichment - Downregulated  LPS LPS.v.None")
dev.off()
png("GO_Enrichment_down.LPS.only.dot.png", width = 800, height = 600)
dotplot(ego.down.LPS.only, showCategory = 20, title = "GO Enrichment - Downregulated  LPS LPS.v.None")
dev.off()


# Save as PDF
pdf("GO_Enrichment_down.LPS.only.dotplot.pdf", width = 10, height = 7)
dotplot(ego.down.LPS.only, showCategory = 10, title = "GO Enrichment - Downregulated LPS LPS.v.None")
dev.off()



#### STILL NEED TO EDIT!!!! 
# plots treatments  

#ego.up.LPS.LPSvIFN
png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()

#ego.up.IFN.LPSvIFN
png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()

#ego.up.PIC.PICvIFN
png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()

#ego.up.IFN.PICvIFN
png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()

#ego.up.PIC.PICvLPS
png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()

#ego.up.LPS.PICvLPS
png("GO_Enrichment_Upregulated_T16.png", width = 800, height = 600)
barplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()
png("GO_Enrichment_Upregulated_T16.dot.png", width = 800, height = 600)
dotplot(ego.up.Shared, showCategory = 20, title = "GO Enrichment - Upregulated in T16")
dev.off()


##################################################################################################

#PCA Plots

####################################################################################################

library(umap)
library(dplyr)


# Load data
PCA.data<-mts.2.3.5.glom_noT0

#shorten column names to treatment only 
colnames(PCA.data) <- sapply(colnames(PCA.data), function(x) {
  strsplit(x, "_")[[1]][2]})


# Extract numerical matrix
expr_matrix <- as.matrix(PCA.data)
expr_matrix_t <- t(expr_matrix)  # Transpose to have samples as rows and genes as columns


# Define conditions (map sample labels to one of the four levels)
raw_conditions <- colnames(PCA.data)

# Create a simplified conditions vector with 4 levels
conditions <- ifelse(grepl("^None", raw_conditions), "None",
                     ifelse(grepl("^IFN", raw_conditions), "IFN",
                            ifelse(grepl("^PIC", raw_conditions), "PIC",
                                   ifelse(grepl("^LPS", raw_conditions), "LPS", NA))))


# Check the condition labels
print(table(conditions))  # Ensure all samples are assigned to one of the four levels
#> print(table(conditions)) 
#conditions
#IFN  LPS None  PIC 
#16   18   18   15

# Perform PCA
pca_result <- prcomp(expr_matrix_t, scale. = TRUE)


# Create a data frame for PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Condition <- conditions  # Add conditions as a new column

###### best PCA plot for treatment - use principal 3 and 4 
# Create the PCA plot
pca.treatment <- ggplot(pca_scores, aes(x = PC3, y = PC4, color = Condition)) +
  geom_point(size = 3) +  # Add points
  labs(title = "PCA of Gene Expression Data",
       x = "Principal Component 3",
       y = "Principal Component 4") +
  theme_minimal() +
  scale_color_manual(values = c("None" = "#5B3794", 
                                "IFN" = "#9e002a", 
                                "LPS" = "#709437", 
                                "PIC" = "#009E73"))


# Set output file
pdf(file="Treatment.PCA.PC3.PC4.pdf", width=7, height=7)


# Plot the chart
print(pca.treatment)


#Save file
dev.off()


# PCA with principal 1 and principal 2 colored by treatment 
# Create the PCA plot
pca.treatment <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +  # Add points
  labs(title = "PCA of Gene Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_manual(values = c("None" = "#5B3794", 
                                "IFN" = "#9e002a", 
                                "LPS" = "#709437", 
                                "PIC" = "#009E73"))


# Set output file
pdf(file="Treatment.PCA.PC1.PC2.pdf", width=7, height=7)


# Plot the chart
print(pca.treatment)


#Save file
dev.off()






##################################################################################################
#maybe most of the varience is explain by patient? 


library(umap)
library(dplyr)


# Load data
PCA.data<-mts.2.3.5.glom_noT0

#shorten column names to treatment only 
colnames(PCA.data) <- sapply(colnames(PCA.data), function(x) {
  strsplit(x, "_")[[1]][1]
})


# Extract numerical matrix
expr_matrix <- as.matrix(PCA.data)
expr_matrix_t <- t(expr_matrix)  # Transpose to have samples as rows and genes as columns


# Define conditions (map sample labels to one of the four levels)
raw_conditions <- colnames(PCA.data)

# Create a simplified conditions vector with 4 levels
conditions <- ifelse(grepl("^MTS002", raw_conditions), "MTS002",
                     ifelse(grepl("^MTS003", raw_conditions), "MTS003",
                            ifelse(grepl("^MTS005", raw_conditions), "MTS005", NA)))


# Check the condition labels
print(table(conditions))  # Ensure all samples are assigned to one of the four levels
#> print(table(conditions)) 
#conditions
#MTS002 MTS003 MTS005 
#24     21     22 

# Perform PCA
pca_result <- prcomp(expr_matrix_t, scale. = TRUE)


# Create a data frame for PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Condition <- conditions  # Add conditions as a new column


# Create the PCA plot
pca.treatment <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +  # Add points
  labs(title = "PCA of Gene Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_manual(values = c("MTS002" = "#5B3794", 
                                "MTS003" = "#9e002a", 
                                "MTS005" = "#709437"))


# Set output file
pdf(file="PCA.by.Pt.PC1.PC2.pdf", width=7, height=7)


# Plot the chart
print(pca.treatment)


#Save file
dev.off()

