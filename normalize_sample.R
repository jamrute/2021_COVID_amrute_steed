library(Seurat)

###### Preprocessing Steps

# Load the dataset
sample <- readRDS('/data/Junedh/COVID_PBMC/merged_annotated_raw.rds')

# Filtering
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# SCTransform
sample <- SCTransform(sample, vars.to.regress = "percent.mt")

saveRDS(sample, file = "/data/Junedh/COVID_PBMC/azimuth/sample_normalized.rds")

