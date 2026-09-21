---
title: "DESeq D4"
author: "Esther Zhou"
date: "2026-04-03"
output: html_document
---

```{r}
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(stringr)
library(matrixStats)
library(knitr)
library(Seurat)
library(corrplot)
library(sparseMatrixStats)
library(SummarizedExperiment)
library(pheatmap)
library(corrplot)
library(ggrepel)
library(DESeq2)
library(EnhancedVolcano)
```


#Load the datasets
```{r}
dds_mouse <- readRDS("/work/users/x/c/xczhou/Chlamydia/processed/dds_mouse.rds")
```

#Excluding M52_D4 and M52_D6
```{r}
dds_mouse_filtered_filtered <- dds_mouse_filtered[, !(colnames(dds_mouse_filtered) %in% c("M52_D4", "M52_D6"))]
```

#Sample-level QC to assess overall similarity between samples
Using PCA and Hierarchical clustering
rlog: 
1. Preprocessing for plotting (PCA, Heatmaps) rather than differential expression testing (use DESeq for testing)
2. It handles the dependence of variance on the mean, improving upon simple transformations.

deseq2::plotPCA is a function in the DESeq2 R package used to visualize sample-to-sample distances using Principal Component Analysis (PCA). It is primarily used for quality control to identify batch effects, outliers, or expected grouping in RNA-seq data,
```{r}
#Run rlog transformation
rld <- rlog(dds_mouse_filtered, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "group")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "day")

PCA_data <- DESeq2::plotPCA(rld, ntop = 500,intgroup = "group", returnData = TRUE)
PCA_data$day <- colData(rld)$day

#Save plot
p_pca <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = group, shape = factor(day))) +
  geom_point(size = 4) +
  theme_bw()

ggsave(
  filename = "/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/PCA_plot_filtered.png",
  plot = p_pca,
  width = 6,
  height = 5,
  dpi = 300
)

```


#Correlation heatmap
sample vs sample correlation

```{r}

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Make sure sample_info rownames match rld_mat columns
rownames(sample_info) <- sample_info$sample

# Annotation data
annotation_col <- sample_info[colnames(rld_mat), c("group", "day"), drop = FALSE]

# Clearer annotation colors
ann_colors <- list(
  group = c(
    Naive = "#0072B2",
    Vaccinated = "#D55E00"
  ),
  day = c(
    "-2" = "#f2f0f7",
    "4"  = "#6a51a3",
    "6"  = "#9e9ac8",
    "11" = "#cbc9e2"
  )
)

# Plot heatmap
heat1_annotation <- pheatmap(
  rld_cor,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 18
)


ggsave(paste0("/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/heatmap_rlog_filtered.tiff"),
       heat1_annotation,
       device = "tiff",
       width =12,
       dpi = 300)
```


#PCA plot
check clustering by group/day
PCA summarizes global expression variation, where each point represents a sample, and clustering indicates similarity between biological conditions.
```{r}
#PCA plot
p1 <- ggplot(PCA_data,
             aes(x=PC1, y=PC2, color=group, shape=day)) +
        geom_point(aes(shape=factor(day), color=group), size=5)+
  scale_color_manual(values=c("blue", "orange"))+
  labs(color = "Group", x = "PC1", y = "PC2", shape = "dpi") + 
  scale_shape_manual(values = c(15, 16, 17, 18)) + 
  theme(axis.text.x = element_text(), #angle = 30
        panel.background = element_rect(fill = "white"),  # Black background
        plot.background = element_rect(fill = "white"),   # Black plot area
        axis.text = element_text(color = "black", size = 16), # Axis text white, bold, larger
        axis.title = element_text(color = "black", size = 18), # Axis titles white, bold, larger
        legend.position = "bottom",  # Bottom-right corner
        legend.justification = c(1, 0),
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.2, 'cm'),  # Reduce horizontal spacing between items
        legend.spacing.y = unit(0.2, 'cm'),  # Reduce vertical spacing between items
        legend.key.size = unit(0.4, 'cm'),   # Make the size of the legend keys smaller
        legend.margin = margin(2, 2, 2, 2),  
        legend.text = element_text(size = 18 ),  # Increase legend font size face = "bold"
        legend.title = element_text(size = 18, color = "black"),  # Bold legend title
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  #legend.box = "horizontal") +
  guides(shape = guide_legend(order = 2),
         color = guide_legend(order = 1))

p1

ggsave(
  "/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/PCA_plot_filtered.tiff",
  p1,
  device = "tiff",
  width = 12,
  height = 8,   # ← increase this (try 8–10)
  dpi = 300
)



```

#MANOVA
Multivariate Analysis of Variance (MANOVA) is a statistical technique used to compare multiple continuous dependent variables across two or more groups (independent variables) simultaneously. It determines if group means differ significantly by considering interrelationships between outcomes, offering greater power than running multiple separate ANOVAs
Qestion: Do group (Naive vs Vaccinated) and day explain variation in the multivariate outcome (PC1, PC2 together)?
H0: The mean vector (PC1, PC2) is the same across groups
H1: At least one group has a different (PC1, PC2) mean
```{r}
### Testing (PC1 & PC2) ~ Group + gpi
manova_model <- manova(cbind(PC1, PC2) ~ group + day, data = PCA_data)

# Get a summary of the MANOVA model
summary(manova_model)

# The summary of each individual ANOVA for the responses
summary.aov(manova_model)

```
## Interpretation:
Variation along PC1 is driven by time (day)
PC2 captures variation associated with group (strongly) and day (also contributes)


#Making a volcano plot for differential expression results
Keep only genes with pvalue < 0.05
label:the most strongly upregulated genes; the most strongly downregulated genes; the most statistically significant genes
Assign colors to gene: purple = significantly downregulated (log2FC < -1); orange = significantly upregulated (log2FC > 1); black = significant, but fold change is between -1 and 1

compute -log10(padj)

```{r}
#Function to visualize DE
plot_DE <- function(All.DE, name){
  
  All.DE <-  All.DE %>% filter(pvalue < 0.05)  %>%
    mutate(label.print = case_when(
      log2FoldChange > quantile(log2FoldChange, 0.975) ~ gene,
      log2FoldChange < quantile(log2FoldChange, 0.025) ~ gene,
      pvalue < quantile(pvalue, 0.05) ~ gene,
      TRUE ~ ""
    ),
    keys = case_when(log2FoldChange < -1 & pvalue < 0.05 ~ 'purple',
                     log2FoldChange > 1 & pvalue < 0.05 ~ "orange",
                     TRUE ~ "black"))
  
  keyvals <- All.DE$keys
  
  names(keyvals)[keyvals == 'orange'] <- 'Up-regulated'
  names(keyvals)[keyvals == 'purple'] <- 'Down-regulated'
  names(keyvals)[keyvals == 'black'] <- 'Significant'
  
  All.DE$log10_padj <- -log10(All.DE$padj)
  
  # Then, get the max value and add 0.5
  max_y <- max(All.DE$log10_padj, na.rm = TRUE) + 0.3
  
  DEplot <- EnhancedVolcano(All.DE,
                            lab = All.DE$gene,
                            selectLab = All.DE$label.print,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            colCustom = keyvals,
                            title = name, 
                            #legendLabels = c('', ''), 
                            drawConnectors = TRUE, 
                            widthConnectors = 0.5,
                            max.overlaps = 20,
                            subtitle = "",
                            colAlpha = 0.7) +
    theme(legend.position = "none",
          plot.title = element_text(
            hjust = 0.5, 
            size = 16,         # make the title larger
            face = "bold", 
            margin = margin(b = -10)  # reduce space below the title
          ))+
    scale_y_continuous(limits = c(1, max_y))
  
  ggsave(paste0("/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/",name,"_DE_filtered.png"), plot = DEplot, dpi = 300, width = 5.3, height = 5)
  # width = 4.3, height = 5
  
}

plot_DE_adj <- function(All.DE, name){
  
  All.DE <-  All.DE %>% filter(padj < 0.05)  %>%
    mutate(label.print = case_when(
      log2FoldChange > quantile(log2FoldChange, 0.975) ~ gene,
      log2FoldChange < quantile(log2FoldChange, 0.025) ~ gene,
      padj < quantile(padj, 0.05) ~ gene,
      TRUE ~ ""
    ),
    keys = case_when(log2FoldChange < -1 & padj < 0.05 ~ 'purple',
                     log2FoldChange > 1 & padj < 0.05 ~ "orange",
                     TRUE ~ "black"))
  
  keyvals <- All.DE$keys
  
  names(keyvals)[keyvals == 'orange'] <- 'Up-regulated'
  
  names(keyvals)[keyvals == 'purple'] <- 'Down-regulated'
  
  names(keyvals)[keyvals == 'black'] <- 'Significant'
  
  All.DE$log10_padj <- -log10(All.DE$padj)
  
  # Then, get the max value and add 0.5
  max_y <- max(All.DE$log10_padj, na.rm = TRUE) + 0.3
  
  DEplot <- EnhancedVolcano(All.DE,
                            lab = All.DE$gene,
                            selectLab = All.DE$label.print,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            colCustom = keyvals,
                            title = name, 
                            #legendLabels = c('', ''), 
                            drawConnectors = TRUE, 
                            widthConnectors = 0.5,
                            max.overlaps = 20,
                            subtitle = "",
                            colAlpha = 0.7) +
    theme(legend.position = "none",
          plot.title = element_text(
            hjust = 0.5, 
            size = 16,         # make the title larger
            face = "bold", 
            margin = margin(b = -10)  # reduce space below the title
          ))+
    scale_y_continuous(limits = c(1, max_y))
  
  ggsave(paste0("/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/",name,"_DE_padj_filtered.png"), plot = DEplot, dpi = 300, width = 5.3, height = 5)
  
}
```



```{r}
#DE Analysis Day4
dds4 <- dds_mouse_filtered[, dds_mouse_filtered$day==4]

# Remove genes with all zero counts in Day 4 samples
dds4 <- dds4[which(rowSums(dds4@assays@data@listData[["counts"]]) != 0),]
dds4

# Set group reference level
dds4$group <- factor(dds4$group, levels=c("Naive","Vaccinated"))

design(dds4) <- ~ group

dds4 <- DESeq(dds4)

plotDispEsts(dds4)

# Check the coefficients for the comparison
resultsNames(dds4)

# Generate results object
# res_vaccine_naive_day4 <- results(dds4, 
#                alpha = 0.05,
#                contrast=c("group","Vaccinated","Naive"))


res_vaccine_naive_day4 <- lfcShrink(dds=dds4 , coef=2, type="apeglm" )

head(res_vaccine_naive_day4 )

res_vaccine_naive_day4_table <- as.data.frame(res_vaccine_naive_day4)
res_vaccine_naive_day4_table$gene <- rownames(res_vaccine_naive_day4_table)

#Remove genes with baseMean = 0
res_vaccine_naive_day4_table <- res_vaccine_naive_day4_table %>% filter(baseMean != 0)
#Check for missing p-values
res_vaccine_naive_day4_table[which(is.na(res_vaccine_naive_day4_table$pvalue)), ]
#Remove genes with NA adjusted p-values
res_vaccine_naive_day4_table <- res_vaccine_naive_day4_table %>% filter(!is.na(padj))
#Count significant genes (raw p-value)
sum(res_vaccine_naive_day4_table$pvalue < 0.05)

plot_DE(All.DE = res_vaccine_naive_day4_table,name="Day 4")

write.csv(res_vaccine_naive_day4_table%>%
            arrange(-log2FoldChange)%>%
            filter(pvalue<0.05 & log2FoldChange > 0), 
          file = "/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/Day4_up_DE_p_filtered.csv")

write.csv(res_vaccine_naive_day4_table%>%
            arrange(log2FoldChange)%>%
            filter(pvalue<0.05 & log2FoldChange < 0), 
          file = "/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/Day4_down_DE_p_filtered.csv")


plot_DE_adj(All.DE = res_vaccine_naive_day4_table%>%dplyr::filter(!is.na(padj)),
            name="Day 4")

write.csv(res_vaccine_naive_day4_table%>%
            arrange(-log2FoldChange)%>%
            filter(padj<0.05 & log2FoldChange > 0), 
          file = "/work/users/x/c/xczhou/Chlamydia/Analysis Outputs/Day4_up_DE_padj_filtered.csv")

write.csv(res_vaccine_naive_day4_table%>%
            arrange(log2FoldChange)%>%
            filter(padj<0.05 & log2FoldChange < 0), 
          file = "/work/users/x/c/xczhou/Chlamydia/Analysis Outputs//Day4_down_DE_padj_filtered.csv")

```




