library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

###### Preprocessing Steps

# Load the dataset
sample <- readRDS("/data/Junedh/COVID_PBMC/azimuth/sample_normalized.rds")

# Load the Reference
reference <- LoadH5Seurat("/data/Junedh/COVID_PBMC/azimuth/pbmc_multimodal.h5seurat")

# Mapping
anchors <- FindTransferAnchors(
  reference = reference,
  query = sample,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

sample <- MapQuery(
  anchorset = anchors,
  query = sample,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

#merge reference and query
reference$id <- 'reference'
sample$id <- 'query'
refquery <- merge(reference, sample)
refquery[["spca"]] <- merge(reference[["spca"]], sample[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)

saveRDS(sample, file = "/data/Junedh/COVID_PBMC/azimuth/mapped/merged_azimuth_mapped.rds")
saveRDS(refquery, file = "/data/Junedh/COVID_PBMC/azimuth/ref_query/refquery.rds")
