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
library(here)
source(here("scripts/utils.R"))
library(Banksy)
library(SpatialExperiment)
library(harmony)

###############################################################################
# load data

rep1=get_xenium_data(path="/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_mouse_brain/Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs/", 
                     mtx_name = "cell_feature_matrix/",
                     trans_name = "transcripts.csv.gz",
                     cells_name="cells.csv.gz")

rep2=get_xenium_data(path="/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_mouse_brain/Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs/", 
                     mtx_name = "cell_feature_matrix/",
                     trans_name = "transcripts.csv.gz",
                     cells_name="cells.csv.gz")

rep3=get_xenium_data(path="/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_mouse_brain/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/", 
                     mtx_name = "cell_feature_matrix/",
                     trans_name = "transcripts.csv.gz",
                     cells_name="cells.csv.gz")

rep1$trans_info = rep1$trans_info[rep1$trans_info$qv >=20 & rep1$trans_info$cell_id != -1 & !(rep1$trans_info$cell_id %in% rep1$zero_cells), ]
rep2$trans_info = rep2$trans_info[rep2$trans_info$qv >=20 &  rep2$trans_info$cell_id != -1 & !(rep2$trans_info$cell_id %in% rep2$zero_cells), ]
rep3$trans_info = rep3$trans_info[rep3$trans_info$qv >=20 & rep3$trans_info$cell_id != -1 & !(rep3$trans_info$cell_id %in% rep3$zero_cells), ]


cm1 = rep1$cm
colnames(cm1) = paste("r1-", colnames(cm1), sep="")

cm2 = rep2$cm
colnames(cm2) = paste("r2-", colnames(cm2), sep="")

cm3 = rep3$cm
colnames(cm3) = paste("r3-", colnames(cm3), sep="")

probe_coords <- as.data.frame(rbind(cbind(sample="replicate1",rep1$trans_info[rep1$trans_info$feature_name %in% rep1$probe, c("feature_name","x","y")]),
                                    cbind(sample="replicate2",rep2$trans_info[rep2$trans_info$feature_name %in% rep2$probe, c("feature_name","x","y")]),
                                    cbind(sample="replicate3",rep3$trans_info[rep3$trans_info$feature_name %in% rep3$probe, c("feature_name","x","y")])
))

codeword_coords <- as.data.frame(rbind(cbind(sample="replicate1",rep1$trans_info[rep1$trans_info$feature_name %in% rep1$codeword, c("feature_name","x","y")]),
                                       cbind(sample="replicate2",rep2$trans_info[rep2$trans_info$feature_name %in% rep2$codeword, c("feature_name","x","y")]),
                                       cbind(sample="replicate3",rep3$trans_info[rep3$trans_info$feature_name %in% rep3$codeword, c("feature_name","x","y")])
))


# List of source directories
sps <- c("1", "2", "3")

spe_list <- lapply(sps, function(id) {
    rp =  get(paste0("rep", id))
    cell_info <-rp$cell_info
    cell_info$cells = paste0("r", id, "-",cell_info$cell_id)
    cm <-get(paste0("cm", id))
    # Subset and rename cells
    sub_info <- cell_info[, c("x_centroid", "y_centroid", "cells")]
    colnames(sub_info) = c("x", "y", "cells")
    rownames(sub_info) <- sub_info$cells
    sub_info = sub_info[colnames(cm),]
    coords <- as.matrix(sub_info[, c("x", "y")])
    
    # Construct SpatialExperiment
    SpatialExperiment(
        assays = list(counts = cm),
        spatialCoords = coords,
        sample_id = paste0("replicate",as.character(id))
    )
})

se <- do.call(cbind, spe_list)

# Convert to Seurat and normalize
seu <- as.Seurat(se, data = NULL)
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu <- FindVariableFeatures(seu, nfeatures = nrow(seu))

# Back to SpatialExperiment
aname <- "logcounts"
# assay(se, aname) <- GetAssayData(seu)

logcounts_mat <- GetAssayData(seu, slot = "data")[, colnames(se)]
assay(se, "logcounts", withDimnames = FALSE) <- logcounts_mat


# Re-split by sample
spe_list <- split(seq_len(ncol(se)), se$sample_id) |>
    lapply(function(cols) se[, cols])

lambda <- c(0.2)
k_geom <- 15
use_agf <- TRUE
compute_agf <- TRUE

spe_list <- lapply(spe_list, computeBanksy, assay_name = aname, 
                   compute_agf = compute_agf, k_geom = k_geom)

spe_joint <- do.call(cbind, spe_list)

spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, 
                          lambda = lambda, group = "sample_id", seed = 1000)

spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, 
                           lambda = lambda, seed = 1000)

cat("Clustering starts\n")
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda,
                           resolution = c(0.1, 0.5), seed = 1000)

spe_joint <- connectClusters(spe_joint)
# 
# spe_joint$banksy_multisample = paste0("c",spe_joint$clust_M1_lam0.2_k50_res0.5)
# spe_joint$banksy_multisample= factor(spe_joint$banksy_multisample,
#                                      levels = paste0("c", 1:length(unique(spe_joint$clust_M1_lam0.2_k50_res0.5))))
# 

clusters <- data.frame(
    x = spatialCoords(spe_joint)[, 1],
    y = spatialCoords(spe_joint)[, 2],
    cell_id = colnames(spe_joint), 
    cluster = paste0("c", colData(spe_joint)$clust_M1_lam0.2_k50_res0.1),  
    sample = colData(spe_joint)$sample_id
)



###############################################################################
# FindMarkers
rownames(clusters) <- clusters$cell_id

# Reorder clusters$cluster to match Seurat cells
Idents(seu) <- clusters$cluster[match(colnames(seu), rownames(clusters))]

usage_fm= peakRAM({
find_markers_result <- FindAllMarkers(seu, only.pos = TRUE,
                                      logfc.threshold = 0.25)
})


###############################################################################\
# limma

all.bct <- factor(clusters$cluster)
sample <- spe_joint$sample_id
y <- DGEList(cbind(cm1, cm2, cm3))
y$genes <-row.names(rep1$cm)

logcounts.all <- normCounts(y,log=TRUE,prior.count=0.1)

design <- model.matrix(~0+all.bct+sample)
colnames(design)[1:(length(levels(all.bct)))] <- levels(all.bct)

mycont <- matrix(0,ncol=length(levels(all.bct)),nrow=length(levels(all.bct)))
colnames(mycont)<-levels(all.bct)
diag(mycont)<-1
mycont[upper.tri(mycont)]<- -1/(length(levels(all.bct))-1)
mycont[lower.tri(mycont)]<- -1/(length(levels(all.bct))-1)

# Fill out remaining rows with 0s
zero.rows <- matrix(0,ncol=length(levels(all.bct)),nrow=(ncol(design)-length(levels(all.bct))))
test <- rbind(mycont,zero.rows)
usage_limma= peakRAM({
fit <- lmFit(logcounts.all,design)
fit.cont <- contrasts.fit(fit,contrasts=test)
fit.cont <- eBayes(fit.cont,trend=TRUE,robust=TRUE)

fit.cont$genes <-row.names(rep1$cm)
limma_dt <- decideTests(fit.cont)
})
summary(limma_dt)

###############################################################################\
# jazzpanda
usage_sv= peakRAM({


all_genes =row.names(rep1$cm)
grid_length=50

# get spatial vectors
all_vectors = get_vectors(x= list("replicate1" = rep1$trans_info,
                                  "replicate2" = rep2$trans_info,
                                  "replicate3" = rep3$trans_info), 
                          sample_names=c("replicate1","replicate2","replicate3"),
                          cluster_info = clusters, bin_type="square",
                          bin_param=c(grid_length,grid_length), 
                          test_genes = all_genes , 
                          n_cores = 5)

nc_vectors = create_genesets(x=list("replicate1" = rep1$trans_info,
                                    "replicate2" = rep2$trans_info,
                                    "replicate3" = rep3$trans_info), 
                             sample_names=c("replicate1","replicate2","replicate3"),
                             name_lst=list(probe=rep1$probe, 
                                           codeword=rep1$codeword),
                             bin_type="square",
                             bin_param=c(grid_length,grid_length), 
                             cluster_info = NULL)
})

usage_glm = peakRAM({
set.seed(188)
jazzPanda_res_lst = lasso_markers(gene_mt=all_vectors$gene_mt,
                                  cluster_mt = all_vectors$cluster_mt,
                                  sample_names=c("replicate1","replicate2","replicate3"),
                                  keep_positive=TRUE, 
                                  background=nc_vectors)

})
#Create a data frame to store the results
results_df <- data.frame(
    technology = "xenium",
    dataset = "Xenium mouse brain",
    genes_n = nrow(cm1),
    transcript_n = nrow(rep1$trans_info)+ nrow(rep2$trans_info)+ nrow(rep3$trans_info),
    cells_n = nrow(clusters),
    clusters_n =length(unique(clusters$cluster)),
    grid_length = grid_length,
    fm_Elapsed_Time_sec = usage_fm$Elapsed_Time_sec,
    fm_Peak_RAM_Used_MiB =usage_fm$Peak_RAM_Used_MiB,
    limma_Elapsed_Time_sec = usage_limma$Elapsed_Time_sec,
    limma_Peak_RAM_Used_MiB =usage_limma$Peak_RAM_Used_MiB,
    sv_Elapsed_Time_sec = usage_sv$Elapsed_Time_sec,
    sv_Peak_RAM_Used_MiB =usage_sv$Peak_RAM_Used_MiB,
    glm_Elapsed_Time_sec = usage_glm$Elapsed_Time_sec,
    glm_Peak_RAM_Used_MiB =usage_glm$Peak_RAM_Used_MiB
)


output_file_name <- "xenium_mouse_brain_5core.csv"

args_all   <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())  # interactive fallback
}

out_dir <- file.path(script_dir, "..", "..", 
                     "data/dataset_computational_complexity")
out_dir <- normalizePath(out_dir)

setwd(out_dir)

saveRDS(clusters,"xenium_mbrain_clusters.Rds")
saveRDS(jazzPanda_res_lst,"xenium_mbrain_jazzPanda_res_lst.Rds")
saveRDS(fit.cont, "xenium_mbrain_fit_cont_obj.Rds")
saveRDS(find_markers_result, "xenium_mbrain_seu_markers.Rds")
saveRDS(seu, "xenium_mbrain_seu.Rds")
saveRDS(spe_joint, "xenium_mbrain_se.Rds")
saveRDS(all_vectors,"xenium_mbrain_sq50_vector_lst.Rds")

# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)

