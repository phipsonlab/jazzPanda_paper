# rerun FindMarkers with default logFC=0.1 for every dataset
library(here)
library(Seurat)
library(peakRAM)
source(here("scripts/utils.R"))

data_nm_lst=c("xenium_hbreast","xenium_mbrain","xenium_hlc",
                  "merscope_hbreast", "cosmx_hhliver", "cosmx_hlc")
data_full_nm= c("xenium_human_breast_cancer", "xenium_mouse_brain",
                "xenium_human_lung_cancer", "merscope_human_breast_cancer",
                "cosmx_human_healthy_liver", "cosmx_human_liver_cancer")
for (data_nm in data_nm_lst){
    idx = which(data_nm_lst==data_nm)
    cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))
    complexity_df = read.csv(here(data_path,paste0(data_full_nm[idx], "_5core.csv")))
    seu = readRDS(here(data_path, paste0(data_nm, "_seu.Rds")))
    
    if (data_nm !="xenium_hlc") {
        rownames(cluster_info) <-cluster_info[[ intersect(c("cell_id", "cells"), 
                                                      colnames(cluster_info))[1] ]]}

    if (data_nm == "xenium_hbreast"){
        colnames(seu) <- sub("^_", "", colnames(seu))
        colnames(seu) <- sub("_2$", "-sp2", colnames(seu))
    }
    print(table(row.names(cluster_info)%in% colnames(seu)))
    # Reorder clusters$cluster to match Seurat cells
    Idents(seu) <- cluster_info$cluster[match(colnames(seu),
                                              rownames(cluster_info))]
    table(Idents(seu))
    set.seed(989)
    usage_fm= peakRAM({
        find_markers_result <- FindAllMarkers(seu, only.pos = TRUE,
                                              logfc.threshold = 0.1)
    })

    complexity_df$fm_Elapsed_Time_sec = usage_fm$Elapsed_Time_sec
    complexity_df$fm_Peak_RAM_Used_MiB =usage_fm$Peak_RAM_Used_MiB
    setwd(data_path)
    saveRDS(find_markers_result, paste0(data_nm, "_seu_markers.Rds"))
    write.csv(complexity_df,
              paste0(data_full_nm[idx], "_5core.csv"),
              row.names = FALSE)
    
}
