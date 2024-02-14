
library(SingleCellExperiment)
load("/nemo/lab/gandhis/data/STPs/babs/outputs/stephanie.strohbuecker/PD_snrnaseq_leiden_sce_anno.RData")


seurat <- as.Seurat(x = leid_sce)

samples_of_interest <- c("PD341", "PD732", "PDC05", "PD366", "PD413", "PDC22", "PD415", "PD416", "PDC34", "PD563", "PD523", "PDC87", "PD678", "PD747", "PD666", "PDC91", "PD531", "PD683")
gene_of_interest <- "ADAR"

cell_types <- seurat$annotation_level_1

FeaturePlot(object = seurat,
            features.plot = gene_of_interest)



AverageExpression(object = seurat, features = gene_of_interest, group.by = cell_types)

pt.size = 0.5,
group.by = cell_types)


# Trial manual conversion
seurat_2 <- CreateSeuratObject(counts = counts(leid_sce),
                               meta.data = as.data.frame(colData(leid_sce)))

seurat_2 <- SetAssayData(object = seurat_2, slot = "data", new.data = logcounts(leid_sce))

# Add feature-level metadata (for Seurat v3 or above)
seurat_2[["RNA"]] <- AddMetaData(seurat_2[["RNA"]], metadata = as.data.frame(rowData(leid_sce)))


# Create subset of seurat data, incorporating just the samples used in the editing analysis (note that DLB samples are not incorporated here).
seurat_subset <- subset(seurat_2, subset = sample_id %in% samples_of_interest)


VlnPlot(seurat_subset, features = "ADAR", group.by = "annotation_level_1", split.by = "disease")


# Create colour pallete for groups in seurat:
num_groups <- length(unique(seurat_2$annotation_level_1)) * length(unique(seurat_2$disease))
color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(3)


DotPlot(seurat_subset, features = "ADAR", group.by = "annotation_level_1",  cols = c("lightgrey", "blue", "darkgreen")) +
    RotatedAxis()
DotPlot(seurat_subset, features = "ADAR", group.by = "annotation_level_1", split.by = "disease", cols = c("lightgrey", "blue", "darkgreen")) +
    RotatedAxis()

RidgePlot()

