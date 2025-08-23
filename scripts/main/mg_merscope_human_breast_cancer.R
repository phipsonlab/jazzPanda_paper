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

data_p = "/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/20231220/"
# load count matrix 

# count_raw =read.csv("/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/HumanBreastCancerPatient1/cell_by_gene.csv")
cm =  fread("/stornext/Bioinf/data/lab_phipson/melody/spaceMarker_analysis/analysis/merscope_processed_count_matrix.csv", header=TRUE)

all_gene_names <- cm[["V1"]]
cm = as.matrix(cm[, 2:ncol(cm)])
row.names(cm) <- all_gene_names

cells_meta =  fread("/stornext/Bioinf/data/lab_phipson/melody/spaceMarker_analysis/analysis/merscope_cell_meta.csv")

cells_meta$cell_id = cells_meta$V1


transcript_df<-fread("/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/20231220/detected_transcripts.csv")
transcript_df$x <- transcript_df$global_x
transcript_df$y <- transcript_df$global_y
transcript_df$feature_name = transcript_df$gene

clusters_info <- as.data.frame(cells_meta[, c("min_x","min_y","max_x", "max_y","leiden","cell_id")])
# clusters_info = cells_meta[, c("min_x","min_y","max_x", "max_y","cell_id")]
# seurat_cluster = as.data.frame(merscope_seu$seurat_clusters)
# seurat_cluster$cell_id = row.names(seurat_cluster)
# colnames(seurat_cluster) = c("cluster", "cell_id")
# clusters_info = merge(clusters_info, seurat_cluster, by = "cell_id", all.x = TRUE)
clusters_info$cluster <- paste("c", clusters_info$leiden,sep="")
# table(clusters_info$cluster)
#clusters_info[clusters_info$cluster %in% c("c2","c3","c4","c5"),"cluster"] = "c2"
#clusters_info$cluster = factor(clusters_info$cluster, levels=c("c0","c1","c2",paste("c", 6:15, sep="")))
clusters_info$cluster <- factor(clusters_info$cluster, levels=paste("c", 0:15, sep=""))
clusters_info$sample <- "sample1"


clusters_info$x <- (clusters_info$min_x + clusters_info$max_x)/2 + 40
clusters_info$y <- (clusters_info$min_y + clusters_info$max_y)/2 + 300
clusters_info <- clusters_info[!duplicated(clusters_info[,c("x","y","cluster")]),]

cluster_names <- paste("c",0:15, sep="")
clusters_info <- na.omit(clusters_info)
dim(clusters_info)

transcript_df$x <- transcript_df$x+ 40
transcript_df$y <- transcript_df$y + 300
# hb_s1 <-transcript_df[,c("x","y","feature_name")]
real_genes <- row.names(cm)
nc_coords <- transcript_df[!(transcript_df$gene %in% real_genes), ]


nc_coords$feature_name <- nc_coords$gene
nc_coords$feature_name <-factor(nc_coords$feature_name)

###############################################################################

grid_length =50

usage_sv= peakRAM({
hbreast_vector_lst<-get_vectors(x=list("sample1" = transcript_df), 
                                sample_names = "sample1",
                                cluster_info = clusters_info,
                                bin_type="square",
                                bin_param=c(grid_length,grid_length),
                                test_genes = row.names(cm),n_cores = 5)
})
###############################################################################
## jazzPanda
### permutation approach
all_genes_names = row.names(cm)

seed_number<-589
set.seed(seed_number)
usage_perm= peakRAM({
perm_p <- compute_permp(x=list("sample1" = transcript_df),
                        cluster_info=clusters_info, 
                        perm.size=5000,
                        bin_type="square",
                        bin_param=c(grid_length,grid_length),
                        test_genes=all_genes_names,
                        correlation_method = "spearman", 
                        n_cores=5, 
                        correction_method="BH")
})
###############################################################################
### linear modelling approach 
nc_names <- unique(nc_coords$feature_name)
nc_coords$sample <- "sample1"
kpt_cols <- c("x","y","feature_name","sample","barcode_id")
nc_coords_mapped = as.data.frame(nc_coords)[,kpt_cols]


nc_vectors <- create_genesets(x=list("sample1" = nc_coords), 
                              sample_names = "sample1",
                              name_lst=list(blanks=nc_names),
                              bin_type="square",
                              bin_param=c(grid_length, grid_length),
                              cluster_info = NULL)

set.seed(589)
usage_glm= peakRAM({
jazzPanda_res_lst <- lasso_markers(gene_mt=hbreast_vector_lst$gene_mt,
                                   cluster_mt = hbreast_vector_lst$cluster_mt,
                                   sample_names=c("sample1"),
                                   keep_positive=TRUE,
                                   background=nc_vectors,
                                   n_fold = 10)

})

###############################################################################
# limma
y <- DGEList(cm[,clusters_info$cell_id])
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
###############################################################################

hbm_seu<-CreateSeuratObject(counts = cm[,clusters_info$cell_id], 
                            project = "hbreast")
Idents(hbm_seu) <- clusters_info[match(colnames(hbm_seu), clusters_info$cell_id),"cluster"]
hbm_seu <- NormalizeData(hbm_seu, verbose = FALSE,
                         normalization.method = "LogNormalize")
hbm_seu <- FindVariableFeatures(hbm_seu, selection.method = "vst", 
                                nfeatures = 1000, verbose = FALSE)
hbm_seu <- ScaleData(hbm_seu, verbose = FALSE)
usage_fm= peakRAM({
seu_markers <- FindAllMarkers(hbm_seu, only.pos = TRUE,logfc.threshold = 0.25)
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


## Output folder name
output_dir_nm <- "dataset_computational_complexity"

out_dir <- file.path(script_dir, output_dir_nm)

setwd(out_dir)


saveRDS(clusters_info,"merscope_hbreast_clusters.Rds")
saveRDS(jazzPanda_res_lst,"merscope_hbreast_jazzPanda_res_lst.Rds")
saveRDS(fit.cont, "merscope_hbreast_fit_cont_obj.Rds")
saveRDS(seu_markers, "merscope_hbreast_seu_markers.Rds")
saveRDS(hbm_seu, "merscope_hbreast_seu.Rds")
saveRDS(perm_p,"merscope_hbreast_perm_lst.Rds")
saveRDS(hbreast_vector_lst,"merscope_hbreast_sq50_vector_lst.Rds")


# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)
