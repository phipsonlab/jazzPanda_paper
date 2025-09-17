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
###############################################################################
# load data

rep1=get_xenium_data(path="/stornext/Bioinf/data/lab_phipson/data/xenium_Mm_brain/xenium_prerelease_jun20_mBrain_replicates/mBrain_ff_rep1/", 
                     mtx_name = "cell_feature_matrix_mtx/")

rep2=get_xenium_data(path="/stornext/Bioinf/data/lab_phipson/data/xenium_Mm_brain/xenium_prerelease_jun20_mBrain_replicates/mBrain_ff_rep2/", 
                     mtx_name = "cell_feature_matrix_mtx/")

rep3=get_xenium_data(path="/stornext/Bioinf/data/lab_phipson/data/xenium_Mm_brain/xenium_prerelease_jun20_mBrain_replicates/mBrain_ff_rep3/", 
                     mtx_name = "cell_feature_matrix_mtx/")

rep1$trans_info = rep1$trans_info[rep1$trans_info$qv >=20 & rep1$trans_info$cell_id != -1 & !(rep1$trans_info$cell_id %in% rep1$zero_cells), ]
rep2$trans_info = rep2$trans_info[rep2$trans_info$qv >=20 &  rep2$trans_info$cell_id != -1 & !(rep2$trans_info$cell_id %in% rep2$zero_cells), ]
rep3$trans_info = rep3$trans_info[rep3$trans_info$qv >=20 & rep3$trans_info$cell_id != -1 & !(rep3$trans_info$cell_id %in% rep3$zero_cells), ]


cm1 = rep1$cm
colnames(cm1) = paste("r1-", colnames(cm1), sep="")
cm1 =cm1[, colnames(cm1)!="r1-6128"]

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
ordered_feature = probe_coords %>% group_by(feature_name) %>% count() %>% arrange(desc(n))%>% pull(feature_name) 
probe_tb = as.data.frame(probe_coords %>% group_by(sample, feature_name) %>% count())
colnames(probe_tb) = c("sample","feature_name","value_count")
probe_tb$feature_name = factor(probe_tb$feature_name, levels= ordered_feature)



codeword_coords <- as.data.frame(rbind(cbind(sample="replicate1",rep1$trans_info[rep1$trans_info$feature_name %in% rep1$codeword, c("feature_name","x","y")]),
                                       cbind(sample="replicate2",rep2$trans_info[rep2$trans_info$feature_name %in% rep2$codeword, c("feature_name","x","y")]),
                                       cbind(sample="replicate3",rep3$trans_info[rep3$trans_info$feature_name %in% rep3$codeword, c("feature_name","x","y")])
))
ordered_feature = codeword_coords %>% group_by(feature_name) %>% count() %>% arrange(desc(n))%>% pull(feature_name) 
codeword_tb = as.data.frame(codeword_coords %>% group_by(sample, feature_name) %>% count())
colnames(codeword_tb) = c("sample","feature_name","value_count")
codeword_tb$feature_name = factor(codeword_tb$feature_name, levels= ordered_feature)

codeword_tb = codeword_tb[order(codeword_tb$feature_name), ]

rep1_seu <- CreateSeuratObject(counts = cm1, project = "replicate1")
rep2_seu <- CreateSeuratObject(counts = cm2, project = "replicate2")
rep3_seu <- CreateSeuratObject(counts = cm3, project = "replicate3")

set.seed(9858)
rep1_seu <- NormalizeData(rep1_seu,normalization.method ="LogNormalize")
rep1_seu = FindVariableFeatures(rep1_seu, selection.method = "vst", 
                                nfeatures = nrow(rep1$cm), verbose = FALSE)

rep2_seu <- NormalizeData(rep2_seu,normalization.method ="LogNormalize")
rep2_seu = FindVariableFeatures(rep2_seu, selection.method = "vst", 
                                nfeatures = nrow(rep2$cm), verbose = FALSE)

rep3_seu <- NormalizeData(rep3_seu,normalization.method ="LogNormalize")
rep3_seu = FindVariableFeatures(rep3_seu, selection.method = "vst", 
                                nfeatures = nrow(rep3$cm), verbose = FALSE)

all_seu<- list(rep1_seu, rep2_seu, rep3_seu)
features <- SelectIntegrationFeatures(object.list = all_seu)
all_seu <- lapply(X = all_seu, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = all_seu, 
                                  anchor.features = features, reduction = "rpca")

# Integrate data
mbrain_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
mbrain_integrated <- ScaleData(mbrain_integrated, verbose = FALSE)
mbrain_integrated <- RunPCA(mbrain_integrated, npcs = 50, verbose = FALSE)

mbrain_integrated <- FindNeighbors(mbrain_integrated, dims = 1:30)
mbrain_integrated <- RunUMAP(mbrain_integrated, dims = 1:30, 
                             verbose = FALSE)
mbrain_integrated <- FindClusters(mbrain_integrated, resolution = 0.05)
#saveRDS(mbrain_integrated, "mbrain_integrated.Rds")

mbrain_integrated$clusters = paste("c", mbrain_integrated$seurat_clusters, sep="")

clusters = as.data.frame(mbrain_integrated$seurat_clusters)
colnames(clusters) = "cluster"
clusters$cluster = paste("c",clusters$cluster, sep="")
clusters$cells = row.names(clusters)
clusters$sample =clusters$cells
clusters$sample =  sub("-.*", "", clusters$cells)
clusters[clusters$sample=="r1","sample"]="replicate1"
clusters[clusters$sample=="r2","sample"]="replicate2"
clusters[clusters$sample=="r3","sample"]="replicate3"
clusters$sample=factor(clusters$sample, 
                       levels=c("replicate1","replicate2","replicate3"))

clusters$x = 0
clusters$y = 0

all_data = list(rep1, rep2, rep3)

for (i in 1:3){
    rp=all_data[[i]]
    rp_nm = paste("replicate", i, sep="")
    curr_cells = rp$cell_info
    row.names(curr_cells)= paste("r",paste(as.character(i),curr_cells$cell_id, sep="-"),sep="")
    curr_cells = curr_cells[row.names(curr_cells) %in% clusters$cells,]
    clusters[clusters$sample==rp_nm, "x"] = curr_cells[match(row.names(curr_cells),clusters[clusters$sample==rp_nm, "cells"] ),"x_centroid"]
    clusters[clusters$sample==rp_nm, "y"] = curr_cells[match(row.names(curr_cells),clusters[clusters$sample==rp_nm, "cells"] ),"y_centroid"]
}


###############################################################################
# FindMarkers

Idents(mbrain_integrated)=clusters$cluster
usage_fm= peakRAM({
find_markers_result <- FindAllMarkers(mbrain_integrated, only.pos = TRUE,
                                      logfc.threshold = 0.25)
})


###############################################################################\
# limma

all.bct <- factor(clusters$cluster)
sample <- mbrain_integrated$orig.ident
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
#treat.all <- treat(fit.cont,lfc=0.25)
#limma_dt <- decideTests(treat.all)
})
summary(limma_dt)

###############################################################################\
# jazzpanda
usage_sv= peakRAM({


all_genes =row.names(rep1$cm)
grid_length=20

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
                             bin_param=c(20,20), 
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
saveRDS(mbrain_integrated, "xenium_mbrain_seu.Rds")
saveRDS(all_vectors,"xenium_mbrain_sq20_vector_lst.Rds")

# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)

