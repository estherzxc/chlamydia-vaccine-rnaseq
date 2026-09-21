#Reading the Cell Ranger files and building a DESeq2 dataset

library(Seurat)
library(Matrix)
library(DESeq2)

#Import data
# Path to Cell Ranger sample folders
folder_path <- "/work/users/x/c/xczhou/Chlamydia/CellRanger_QC_New"

# Get all sample directories
sample_dirs <- list.files(folder_path, pattern = "^M.*CellRanger$", full.names = TRUE)

# Read the first sample to get the full gene list
tmp <- Read10X_h5(
  file.path(sample_dirs[1], "outs", "raw_feature_bc_matrix.h5"),
  use.names = TRUE,
  unique.features = TRUE
)

# Create an empty count matrix: rows = genes, columns = samples
all_data <- matrix(0, nrow = nrow(tmp), ncol = length(sample_dirs))
rownames(all_data) <- rownames(tmp)
colnames(all_data) <- basename(sample_dirs)

# Loop through each sample and sum counts across all barcodes
for (i in seq_along(sample_dirs)) {
  
  mat <- Read10X_h5(
    file.path(sample_dirs[i], "outs", "raw_feature_bc_matrix.h5"),
    use.names = TRUE,
    unique.features = TRUE
  )
  
  # Collapse barcode-level counts into one count per gene per sample
  all_data[, i] <- as.matrix(rowSums(mat))
}

# Separate mouse and chlamydia genes
##Keep only mouse genes for host differential expression analysis
#Keep Chlamydia genes for QC
mouse_data <- all_data[!grepl("^Chlamydia_", rownames(all_data)), ]
chlam_data <- all_data[grepl("^Chlamydia_", rownames(all_data)), ]

dim(mouse_data)
dim(chlam_data)


#Clean the gene names and sample names
rownames(mouse_data) <- sub("^Mouse_+", "", rownames(mouse_data))
colnames(mouse_data) <- sub("^(M[0-9]+_D[0-9]+).*", "\\1", colnames(mouse_data))
rownames(chlam_data) <- sub("^Chlamydia_+", "", rownames(chlam_data))
colnames(chlam_data) <- sub("^(M[0-9]+_D[0-9]+).*", "\\1", colnames(chlam_data))



# Create metadata from sample names
sample_info <- data.frame(
  sample = colnames(mouse_data),
  stringsAsFactors = FALSE
)


# Extract mouse number (e.g., 80 from M80_D6)
sample_info$mouse <- as.numeric(sub("^M([0-9]+)_D.*", "\\1", sample_info$sample))


# Extract day
sample_info$day <- as.numeric(sub("^M[0-9]+_D([0-9]+)$", "\\1", sample_info$sample))
# Recode D0 to D-2
sample_info$day[sample_info$day == 0] <- -2
## Convert day to factor
sample_info$day <- factor(sample_info$day)
# Create a day label for plotting
sample_info$day_label <- paste0("D", sample_info$day)


# Create group variable from mouse number
sample_info$group <- ifelse(sample_info$mouse >= 71, "Naive", "Vaccinated")
# Convert group to factor and set reference level
sample_info$group <- factor(sample_info$group, levels = c("Naive", "Vaccinated"))


# Set row names to match count matrix columns
rownames(sample_info) <- sample_info$sample

# Check that metadata rows match count matrix columns exactly
all(rownames(sample_info) == colnames(mouse_data))

#Create Seurat Objects
# Create a Seurat object using the mouse gene counts
# counts = the mouse portion of the all_data matrix
# rows = mouse genes
# columns = swab samples
# project name is just a label for the Seurat object
# min.cells = 0 and min.features = 0 means we do NOT filter anything out
# assay = "MouseRNA" names the assay that stores mouse gene counts
data.seurat <- CreateSeuratObject(counts = mouse_data, 
                                  project = "Dual-seq", 
                                  min.cells = 0,
                                  min.features = 0, 
                                  assay = "MouseRNA")

# Create a second assay object that stores Chlamydia gene counts
chlamydia_assay <- CreateAssayObject(counts = chlam_data, 
                                     min.cells = 0, 
                                     min.features = 0)
# Add the Chlamydia assay into the Seurat object
# Now the Seurat object contains TWO assays:
#   1) MouseRNA
#   2) ChlamydiaRNA
data.seurat[["ChlamydiaRNA"]] <- chlamydia_assay

# Add metadata to Seurat
all(rownames(sample_info) == colnames(mouse_data))
data.seurat <- AddMetaData(
  object = data.seurat,
  metadata = sample_info
)

data.seurat
dim(data.seurat[["MouseRNA"]])
dim(data.seurat[["ChlamydiaRNA"]])



# Create DESeq2 dataset
dds_mouse <- DESeqDataSetFromMatrix(
  countData = mouse_data,
  colData = sample_info,
  design = ~ group + day
)

# Save Seurat object
saveRDS(data.seurat, 
        file = "/work/users/x/c/xczhou/Chlamydia/processed/data_seurat.rds")

# Save DESeq2 object
saveRDS(dds_mouse, 
        file = "/work/users/x/c/xczhou/Chlamydia/processed/dds_mouse.rds")

# Also save matrices for quick access
saveRDS(mouse_data, 
        file = "/work/users/x/c/xczhou/Chlamydia/processed/mouse_data.rds")

saveRDS(sample_info, 
        file = "/work/users/x/c/xczhou/Chlamydia/processed/sample_info.rds")
