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

path="/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_lung_cancer/"
hl_cancer =get_xenium_data(path=path, 
                           mtx_name = "cell_feature_matrix/", trans_name="transcripts.csv.gz", cells_name="cells.csv.gz")


una_tran = nrow(hl_cancer$trans_info[hl_cancer$trans_info$cell_id == "UNASSIGNED", ])
lqc_tran = nrow(hl_cancer$trans_info[hl_cancer$trans_info$qv <20, ])
hl_cancer$trans_info = hl_cancer$trans_info[hl_cancer$trans_info$cell_id != "UNASSIGNED" & 
                                                !(hl_cancer$trans_info$cell_id %in% hl_cancer$zero_cells) & 
                                                hl_cancer$trans_info$qv>=20, ]

# create seurat object
set.seed(20394)

biorep=factor(rep(c("hlc"),c(ncol(hl_cancer$cm))))
names(biorep) = colnames( hl_cancer$cm)
hlc_seu=CreateSeuratObject(counts = hl_cancer$cm, project = "hlc")
hlc_seu=AddMetaData(object=hlc_seu, metadata = biorep, col.name="biorep")
hlc_seu = NormalizeData(hlc_seu, verbose = FALSE,
                        normalization.method = "LogNormalize")
hlc_seu = FindVariableFeatures(hlc_seu, selection.method = "vst", 
                               nfeatures = 392, verbose = FALSE)

hlc_seu=ScaleData(hlc_seu, verbose = FALSE)
hlc_seu=RunPCA(hlc_seu, npcs = 30, verbose = FALSE, 
               features = row.names(hlc_seu))
# print(ElbowPlot(seu, ndims = 50))
hlc_seu <- RunUMAP(object = hlc_seu, dims = 1:20)


# load provided cluster labels
graphclust_off=read.csv(paste(path, "graphclust.csv",sep=""))
graphclust_off$new_cluster = graphclust_off$Cluster
graphclust_off[graphclust_off$Cluster %in% c(7,10,11,17),"new_cluster"] = "6"
graphclust_off[graphclust_off$Cluster==14,"new_cluster"] = "4"
graphclust_off[graphclust_off$Cluster==21,"new_cluster"] = "18"
graphclust_off[graphclust_off$Cluster %in% c(20,22),"new_cluster"] = "15"
graphclust_off[graphclust_off$Cluster %in% c(19,23),"new_cluster"] = "5"

## refined cluster
rep_clusters = as.data.frame(paste("c",graphclust_off$new_cluster, sep=""))
row.names(rep_clusters) = graphclust_off$Barcode
colnames(rep_clusters) = "cluster"
rep_clusters$cluster=factor(rep_clusters$cluster)
cells= hl_cancer$cell_info
row.names(cells) = cells$cell_id
rp_names =  row.names(rep_clusters) 
rep_clusters[rp_names,"x"] = cells[match(rp_names,row.names(cells)),"x"]
rep_clusters[rp_names,"y"] = cells[match(rp_names,row.names(cells)),"y"]

clusters = rep_clusters
clusters$sample = "hl_cancer"
###############################################################################
## permutation approach

seed_number=589
grid_length=70

set.seed(seed_number)
usage_perm= peakRAM({
perm_p = compute_permp(x=list("hl_cancer" = hl_cancer$trans_info),
                       cluster_info=clusters, 
                       perm.size=5000,
                       bin_type="square",
                       bin_param=c(grid_length,grid_length),
                       test_genes=row.names(hl_cancer$cm),
                       correlation_method = "spearman", 
                       n_cores=5, 
                       correction_method="BH")
})
###############################################################################
## linear modelling appraoch 
clusters = rep_clusters
clusters$sample = "hl_cancer"


all_genes = row.names(hl_cancer$cm)
grid_length=70
usage_sv= peakRAM({
rep1_sq70_vectors = get_vectors(x=list("hl_cancer" = hl_cancer$trans_info),
                                sample_names="hl_cancer",
                                cluster_info = clusters,bin_type="square",
                                test_genes = all_genes, 
                                bin_param=c(grid_length,grid_length), 
                                n_cores = 5)

rep1_bc_vectors = create_genesets(x=list("hl_cancer" = hl_cancer$trans_info),
                                  sample_names="hl_cancer", 
                                  name_lst=list(probe = hl_cancer$probe,
                                                codeword = hl_cancer$codeword), 
                                  bin_type="square",
                                  bin_param=c(70, 70), 
                                  cluster_info=NULL)
})

set.seed(seed_number)
usage_glm= peakRAM({
rep1_lasso_lst= lasso_markers(gene_mt=rep1_sq70_vectors$gene_mt,
                              cluster_mt = rep1_sq70_vectors$cluster_mt,
                              sample_names=c("hl_cancer"),
                              keep_positive=TRUE,
                              background=rep1_bc_vectors,
                              n_fold = 10)
})

###############################################################################
## FindMarkers

Idents(hlc_seu) = rep_clusters[match(colnames(hlc_seu), 
                                     row.names(rep_clusters)),"cluster"]
usage_fm= peakRAM({
seu_markers <- FindAllMarkers(hlc_seu, only.pos = TRUE,logfc.threshold = 0.1)
table(seu_markers$cluster)

})

###############################################################################
## limma 

y <- DGEList(hl_cancer$cm[,setdiff(colnames(hl_cancer$cm),
                                   c(hl_cancer$zero_cells,"cate"))])
y$genes <-row.names(hl_cancer$cm)

logcounts <- speckle::normCounts(y,log=TRUE,prior.count=0.1)
maxclust <- length(unique(rep_clusters$cluster))

grp <- rep_clusters$cluster

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

#treat_res <- treat(fit.cont,lfc=0.5)
limma_dt<-decideTests(fit.cont)
})
###############################################################################

#Create a data frame to store the results
results_df <- data.frame(
    technology = "xenium",
    dataset = "Xenium human lung cancer",
    genes_n = nrow(hl_cancer$cm),
    transcript_n =  nrow(hl_cancer$trans_info),
    cells_n = ncol(hl_cancer$cm),
    clusters_n =length(unique(clusters$cluster)),
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


output_file_name <- "xenium_human_lung_cancer_5core.csv"


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


saveRDS(clusters,"xenium_hlc_clusters.Rds")
saveRDS(rep1_lasso_lst,"xenium_hlc_jazzPanda_res_lst.Rds")
saveRDS(fit.cont, "xenium_hlc_fit_cont_obj.Rds")
saveRDS(seu_markers, "xenium_hlc_seu_markers.Rds")
saveRDS(hlc_seu, "xenium_hlc_seu.Rds")
saveRDS(perm_p,"xenium_hlc_perm_lst.Rds")
saveRDS(rep1_sq70_vectors,"xenium_hlc_sq70_vector_lst.Rds")

write.csv(results_df, output_file_name, row.names = FALSE)


