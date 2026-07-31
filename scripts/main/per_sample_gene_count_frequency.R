library(jazzPanda)
library(Seurat)
library(data.table)
library(dplyr)
library(glmnet)
library(caret)
library(limma)
library(edgeR)
library(speckle)
library(peakRAM)
library(dplyr) 
library(scCustomize)
library(SpatialExperimentIO)
library(SpatialExperiment)
library(here)
library(Matrix)
source(here("scripts/utils.R"))
args_all   = commandArgs(trailingOnly = FALSE)
script_arg = grep("^--file=", args_all, value = TRUE)
script_dir = if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())
}
project_root = normalizePath(file.path(script_dir, "..", ".."))


# Per-value frequency table for one count matrix (dense base matrix or sparse Matrix).
# Zero bin computed as prod(dim) - nnz (double) so it never overflows the integer range.
count_freq = function(cm, dataset_nm){
    if (is(cm, "sparseMatrix")){
        cm = drop0(cm)
        nz = cm@x
    } else {
        nz = cm[cm != 0]
    }
    total  = prod(dim(cm))
    counts = tabulate(nz + 1L)
    if (length(counts) == 0L) counts = 0
    counts[1] = total - length(nz)
    names(counts) = 0:(length(counts) - 1L)
    counts = counts[counts > 0]
    prop   = counts / sum(counts)
    data.frame(
        sample = dataset_nm,
        value  = names(counts),
        count  = as.numeric(counts),
        prop   = round(prop, 4),
        label  = paste0(format(counts, big.mark = ",", scientific = FALSE, trim = TRUE),
                        " (", round(100 * prop, 4), "%)"),
        row.names = NULL
    )
}

summary_tab = data.frame()

# ---- Xenium human breast ----
sp1_path = "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_breast_samples/Xenium_hbreast_sample1/"
sp1 = get_xenium_data(sp1_path, mtx_name = "cell_feature_matrix",
                      trans_name = "transcripts.csv.gz", cells_name = "cells.csv.gz")
sp1$trans_info = sp1$trans_info[sp1$trans_info$qv >= 20 &
                                    sp1$trans_info$cell_id != -1 &
                                    !(sp1$trans_info$cell_id %in% sp1$zero_cells), ]

sp2_path = "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_breast_samples/Xenium_hbreast_sample2/"
sp2 = get_xenium_data(sp2_path, mtx_name = "cell_feature_matrix",
                      trans_name = "transcripts.csv.gz", cells_name = "cells.csv.gz")
sp2$trans_info = sp2$trans_info[sp2$trans_info$qv >= 20 &
                                    sp2$trans_info$cell_id != -1 &
                                    !(sp2$trans_info$cell_id %in% sp2$zero_cells), ]

summary_tab = rbind(summary_tab, count_freq(sp1$cm, "xenium_hbreast_sp1"))
summary_tab = rbind(summary_tab, count_freq(sp2$cm, "xenium_hbreast_sp2"))

# ---- MERSCOPE human breast ----
se = readMerscopeSXE(dirName = "/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/HumanBreastCancerPatient1/",
                     countMatPattern = "cell_by_gene.csv", metaDataPattern = "cell_metadata.csv")
x_avg  = (se@colData$min_x + se@colData$max_x) / 2
y_avg  = (se@colData$min_y + se@colData$max_y) / 2
coords = cbind(x = x_avg, y = y_avg)
SpatialExperiment::spatialCoords(se) = coords

se = se[-grep("^Blank", rownames(se)), ]
lib_size = Matrix::colSums(assay(se, "counts"))
colData(se)$lib_size = lib_size
se = se[, lib_size > 30 & lib_size < 2500]

seu = as.Seurat(se, data = NULL)
cm  = GetAssayData(seu, assay = "originalexp", layer = "counts")
summary_tab = rbind(summary_tab, count_freq(cm, "merscope_sp1"))

# ---- CosMx liver (normal + cancer) ----
seu = readRDS("/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/CosMx_normal_and_diseased_liver_samples/LiverDataReleaseSeurat_newUMAP.RDS")
metadata = as.data.frame(seu@meta.data)
qc_cols = c("qcFlagsRNACounts", "qcFlagsCellCounts", "qcFlagsCellPropNeg",
            "qcFlagsCellComplex", "qcFlagsCellArea", "qcFlagsFOV")
metadata = metadata[!apply(metadata[, qc_cols], 1, function(x) any(x == "Fail")), ]
cell_info_cols = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm",
                   "nCount_negprobes", "nFeature_negprobes", "nCount_falsecode",
                   "nFeature_falsecode", "slide_ID_numeric", "Run_Tissue_name",
                   "fov", "cellType", "niche", "cell_id")
cellCoords = metadata[, cell_info_cols]
counts = seu[["RNA"]]@counts

normal_cells = row.names(cellCoords[cellCoords$slide_ID_numeric == 1, ])
cancer_cells = row.names(cellCoords[cellCoords$slide_ID_numeric == 2, ])

summary_tab = rbind(summary_tab, count_freq(counts[, normal_cells], "cosmx_healthy"))
summary_tab = rbind(summary_tab, count_freq(counts[, cancer_cells], "cosmx_cancer"))

out_dir = file.path(project_root, "data")
out_dir = normalizePath(out_dir)
write.csv(summary_tab,  file.path(out_dir,"gene_count_frequency.csv"),
          row.names = FALSE)

#################
