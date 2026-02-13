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
library(SpatialExperimentIO)
library(SpatialExperiment)
library(Seurat)
library(Banksy)
source(here("scripts/utils.R"))

se =  readMerscopeSXE(dirName = "/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/HumanBreastCancerPatient1/",
                      countMatPattern = "cell_by_gene.csv", metaDataPattern = "cell_metadata.csv")
x_avg <- (se@colData$min_x + se@colData$max_x) / 2
y_avg <- (se@colData$min_y + se@colData$max_y) / 2

# create a matrix
coords <- cbind(x = x_avg, y = y_avg)

# assign to spatialCoordiantes
SpatialExperiment::spatialCoords(se) <- coords

blank_genes <- grep("^Blank", rownames(se), value = TRUE)


blank_idx <- grep("^Blank", rownames(se))
se <- se[-blank_idx, ]

lib_size <- Matrix::colSums(assay(se, "counts"))
colData(se)$lib_size <- lib_size

# Filter cells with library size between 30 and 2500
se <- se[, lib_size > 30 & lib_size < 2500]

seu <- as.Seurat(se, data = NULL)
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu <- FindVariableFeatures(seu, nfeatures = nrow(seu))

# copy log-normalized data back to se
aname <- "logcounts"
logcounts_mat <- GetAssayData(seu, slot = "data")[, colnames(se)]
assay(se, aname, withDimnames = FALSE) <- logcounts_mat

lambda <- 0.2
k_geom <- 15
use_agf <- TRUE

se <- computeBanksy(se, assay_name = aname,
                    compute_agf = TRUE, k_geom = k_geom)

se <- runBanksyPCA(se, use_agf = use_agf, lambda = lambda, seed = 1000)
se <- runBanksyUMAP(se, use_agf = use_agf, lambda = lambda, seed = 1000)

cat("Clustering starts\n")
se <- clusterBanksy(se, use_agf = use_agf, lambda = lambda,
                    resolution = c(0.5, 0.8), seed = 1000)

se <- connectClusters(se)

# saveRDS(se, "merscope_hbreast_se.Rds")

# se = readRDS("/vast/projects/xenium_5k/jazzPanda_paper/scripts/main/merscope_hbreast_se.Rds")
clusters_info <- data.frame(
    x = spatialCoords(se)[, 1],
    y = spatialCoords(se)[, 2],
    cell_id = colnames(se), 
    cluster = paste0("c", colData(se)$clust_M1_lam0.2_k50_res0.5),  
    sample ="sample01"
)

clusters_info$cluster = factor(clusters_info$cluster,
                               levels = paste0("c",sort(unique(colData(se)$clust_M1_lam0.2_k50_res0.5))))
cat("Loading transcripts\n")
transcript_df<-fread(paste0(MERSCOPE_RAW_DATA,"merscope_hbc_detected_transcripts.csv"))
transcript_df$x <- transcript_df$global_x
transcript_df$y <- transcript_df$global_y
transcript_df$gene = make.names(transcript_df$gene)
transcript_df$feature_name = transcript_df$gene
cat("Transcripts loaded \n")
real_genes <- row.names(se)
nc_coords <- transcript_df[!(transcript_df$gene %in% real_genes), ]

nc_coords$feature_name <- nc_coords$gene
nc_coords$feature_name <-factor(nc_coords$feature_name)

###############################################################################

grid_length =50
cat("Running get_vectors \n")
usage_sv= peakRAM({
hbreast_vector_lst<-get_vectors(x=list("sample01" = transcript_df), 
                                sample_names = "sample01",
                                cluster_info = clusters_info,
                                bin_type="square",
                                bin_param=c(grid_length,grid_length),
                                test_genes = row.names(se),
                                n_cores = 5)
})
cat("get_vectors completed \n")
###############################################################################
#saveRDS(hbreast_vector_lst, "hbreast_vector_lst.Rds")
### linear modelling approach 
nc_names <- unique(nc_coords$feature_name)
nc_coords$sample <- "sample01"
kpt_cols <- c("x","y","feature_name","sample","barcode_id")
nc_coords_mapped = as.data.frame(nc_coords)[,kpt_cols]


nc_vectors <- create_genesets(x=list("sample01" = nc_coords), 
                              sample_names = "sample01",
                              name_lst=list(blanks=nc_names),
                              bin_type="square",
                              bin_param=c(grid_length, grid_length),
                              cluster_info = NULL)
saveRDS(nc_vectors, "nc_vectors.Rds")
set.seed(589)
cat("Running lasso_markers \n")
usage_glm= peakRAM({
jazzPanda_res_lst <- lasso_markers(gene_mt=hbreast_vector_lst$gene_mt,
                                   cluster_mt = hbreast_vector_lst$cluster_mt,
                                   sample_names=c("sample01"),
                                   keep_positive=TRUE,
                                   background=nc_vectors,
                                   n_fold = 10)

})
cat("lasso_markers completed \n")
###############################################################################
## jazzPanda
### permutation approach
all_genes_names = row.names(se)
cat("Running compute_permp \n")
seed_number<-589
set.seed(seed_number)
usage_perm= peakRAM({
    perm_p <- compute_permp(x=list("sample01" = transcript_df),
                            cluster_info=clusters_info, 
                            perm.size=5000,
                            bin_type="square",
                            bin_param=c(grid_length,grid_length),
                            test_genes=all_genes_names,
                            correlation_method = "spearman", 
                            n_cores=5, 
                            correction_method="BH")
})
cat("compute_permp completed \n")
###############################################################################
# limma
cat("Running limma \n")
cm <- assay(se, "counts")

y <- DGEList(cm)
y$genes <-row.names(cm)

logcounts <- speckle::normCounts(y,log=TRUE,prior.count=0.1)
maxclust <- length(unique(clusters_info$cluster))

grp <- clusters_info$cluster

design <- model.matrix(~0+grp)
colnames(design) <- levels(grp)

mycont <- matrix(NA,ncol=length(levels(grp)),nrow=length(levels(grp)))
rownames(mycont)<-colnames(mycont)<-levels(grp)
diag(mycont)<-1
mycont[upper.tri(mycont)]<- -1/(length(levels(factor(grp)))-1)
mycont[lower.tri(mycont)]<- -1/(length(levels(factor(grp)))-1)
usage_limma= peakRAM({ 
fit <- lmFit(logcounts,design)
fit.cont <- contrasts.fit(fit,contrasts=mycont)
fit.cont <- eBayes(fit.cont,trend=TRUE,robust=TRUE)

limma_dt<-decideTests(fit.cont)

})
cat("limma completed \n")
###############################################################################

hbm_seu<-CreateSeuratObject(counts = cm,
                            project = "hbreast")
Idents(hbm_seu) <- clusters_info[match(colnames(hbm_seu), clusters_info$cell_id),"cluster"]
hbm_seu <- NormalizeData(hbm_seu, verbose = FALSE,
                         normalization.method = "LogNormalize")
hbm_seu <- FindVariableFeatures(hbm_seu, selection.method = "vst", 
                                nfeatures = 1000, verbose = FALSE)
hbm_seu <- ScaleData(hbm_seu, verbose = FALSE)

# set.seed(989)
# # print(ElbowPlot(seu, ndims = 50))
# hbm_seu <- RunPCA(hbm_seu, features = row.names(hbm_seu), 
#                   npcs = 50, verbose = FALSE)
# 
# hbm_seu <- RunUMAP(object = hbm_seu, dims = 1:20)

usage_fm= peakRAM({
seu_markers <- FindAllMarkers(hbm_seu, only.pos = TRUE,logfc.threshold = 0.1)
})

###############################################################################
results_df <- data.frame(
    technology = "MERSCOPE",
    dataset = "MERSCOPE human breast cancer",
    genes_n =  nrow(cm),
    transcript_n =  nrow(transcript_df),
    cells_n =  ncol(cm),
    clusters_n =length(unique(clusters_info$cluster)),
    grid_length = grid_length,
    fm_Elapsed_Time_sec = usage_fm$Elapsed_Time_sec,
    fm_Peak_RAM_Used_MiB =usage_fm$Peak_RAM_Used_MiB,
    limma_Elapsed_Time_sec = usage_limma$Elapsed_Time_sec,
    limma_Peak_RAM_Used_MiB =usage_limma$Peak_RAM_Used_MiB,
    perm_Elapsed_Time_sec = usage_perm$Elapsed_Time_sec,
    perm_Peak_RAM_Used_MiB =usage_perm$Peak_RAM_Used_MiB,
    sv_Elapsed_Time_sec = usage_sv$Elapsed_Time_sec,
    sv_Peak_RAM_Used_MiB =usage_sv$Peak_RAM_Used_MiB,
    glm_Elapsed_Time_sec = usage_glm$Elapsed_Time_sec,
    glm_Peak_RAM_Used_MiB =usage_glm$Peak_RAM_Used_MiB
)


output_file_name <- "merscope_human_breast_cancer_5core.csv"

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


saveRDS(clusters_info,"merscope_hbreast_clusters.Rds")
saveRDS(jazzPanda_res_lst,"merscope_hbreast_jazzPanda_res_lst.Rds")
saveRDS(fit.cont, "merscope_hbreast_fit_cont_obj.Rds")
saveRDS(seu_markers, "merscope_hbreast_seu_markers.Rds")
saveRDS(hbm_seu, "merscope_hbreast_seu.Rds")
saveRDS(se, "merscope_hbreast_se.Rds")
saveRDS(perm_p,"merscope_hbreast_perm_lst.Rds")
saveRDS(hbreast_vector_lst,"merscope_hbreast_sq50_vector_lst.Rds")


# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)
