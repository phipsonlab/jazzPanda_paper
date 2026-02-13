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
library(here)
source(here("scripts/utils.R"))
args_all   <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())  # interactive fallback
}


project_root <- normalizePath(file.path(script_dir, "..", ".."))


sp1_path <- "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_breast_samples/Xenium_hbreast_sample1/"

sp1 = get_xenium_data(sp1_path, 
                      mtx_name="cell_feature_matrix",
                      trans_name="transcripts.csv.gz", 
                      cells_name="cells.csv.gz" )
sp1$trans_info = sp1$trans_info[sp1$trans_info$qv >=20 & 
                                    sp1$trans_info$cell_id != -1 & 
                                    !(sp1$trans_info$cell_id %in% sp1$zero_cells), ]

## sample2
sp2_path <- "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_breast_samples/Xenium_hbreast_sample2/"
sp2 = get_xenium_data(sp2_path, 
                      mtx_name="cell_feature_matrix",
                      trans_name="transcripts.csv.gz", 
                      cells_name="cells.csv.gz" )

sp2$trans_info = sp2$trans_info[sp2$trans_info$qv >=20 & 
                                    sp2$trans_info$cell_id != -1 & 
                                    !(sp2$trans_info$cell_id %in% sp1$zero_cells), ]



graphclust_sp1=read.csv(file.path(project_root, "data", "Xenium_rep1_supervised_celltype.csv"))
graphclust_sp1$anno =as.character(graphclust_sp1$Cluster)
target_clusters = c("Tumor", "Stromal","Macrophages","Myoepithelial",
                    "T_Cells", "B_Cells","Endothelial", "Dendritic", "Mast_Cells")

t_cells =  c("CD4+_T_Cells","CD8+_T_Cells","T_Cell_&_Tumor_Hybrid",
             "Stromal_&_T_Cell_Hybrid")
dc_cells = c("LAMP3+_DCs","IRF7+_DCs")
macro_cells = c("Macrophages_1","Macrophages_2")
myo_cells = c("Myoepi_KRT15+", "Myoepi_ACTA2+")
tumor_cells = c("Invasive_Tumor", "Prolif_Invasive_Tumor")
graphclust_sp1[graphclust_sp1$Cluster %in% t_cells, "anno"] = "T_Cells"
graphclust_sp1[graphclust_sp1$Cluster %in% macro_cells,"anno"] = "Macrophages"
graphclust_sp1[graphclust_sp1$Cluster %in% c("DCIS_1", "DCIS_2"),"anno"] = "DCIS"
graphclust_sp1[graphclust_sp1$Cluster %in% dc_cells,"anno"] = "Dendritic"
graphclust_sp1[graphclust_sp1$Cluster %in% myo_cells,"anno"] = "Myoepithelial"
graphclust_sp1[graphclust_sp1$Cluster %in% tumor_cells,"anno"] = "Tumor"

graphclust_sp1$anno=factor(graphclust_sp1$anno,
                           levels=c("Tumor", "DCIS", "Stromal",
                                    "Macrophages","Myoepithelial",
                                    "T_Cells", "B_Cells",
                                    "Endothelial", 
                                    "Dendritic", "Mast_Cells",
                                    "Perivascular-Like","Unlabeled"))


sp1_clusters = as.data.frame(cbind(as.character(graphclust_sp1$anno),
                                   paste("c",as.numeric(factor(graphclust_sp1$anno)), 
                                         sep="")))

row.names(sp1_clusters) = graphclust_sp1$Barcode
colnames(sp1_clusters) = c("anno","cluster")

cells= sp1$cell_info
row.names(cells) =cells$cell_id
rp_names =  row.names(sp1_clusters) 
sp1_clusters[rp_names,"x"] = cells[match(rp_names,row.names(cells)),
                                   "x"]
sp1_clusters[rp_names,"y"] = cells[match(rp_names,row.names(cells)),
                                   "y"]

sp1_clusters$anno=factor(sp1_clusters$anno,
                         levels=c("Tumor", "DCIS", "Stromal",
                                  "Macrophages","Myoepithelial",
                                  "T_Cells", "B_Cells",
                                  "Endothelial", 
                                  "Dendritic", "Mast_Cells",
                                  "Perivascular-Like","Unlabeled"))
sp1_clusters$cluster=factor(sp1_clusters$cluster,
                            levels=paste("c", 1:12, sep=""))

sp1_clusters$cells =paste(row.names(sp1_clusters),"_1",sep="")
sp1_clusters$sample="sample1"
selected_sp1 = sp1_clusters
selected_sp1$cluster_comb = as.character(selected_sp1$anno)
selected_sp1[selected_sp1$cluster_comb =="DCIS","cluster_comb"]  = "Tumor"
selected_sp1 = selected_sp1[selected_sp1$cluster_comb %in% target_clusters,]
selected_sp1$cluster_comb = factor(selected_sp1$cluster_comb, levels=target_clusters)
selected_sp1$cluster =paste("c",as.numeric(factor(selected_sp1$cluster_comb)),
                            sep="")
selected_sp1 = selected_sp1[,c("cluster","x","y","cells","sample","cluster_comb")]
colnames(selected_sp1)[6]="anno"

sample2_anno = read.csv(file.path(project_root, "data", "Sample2_Xenium_cell_type_manual.csv"),
                        row.names = 1)
sp2_clusters = as.data.frame(cbind(row.names(sample2_anno),
                                   as.character(sample2_anno$cell_type)))
colnames(sp2_clusters) = c("cells","anno")
row.names(sp2_clusters) = sp2_clusters$cells
cells= sp2$cell_info
row.names(cells) = cells$cell_id
rp_names = sp2_clusters$cells
sp2_clusters[rp_names,"x"] = cells[match(rp_names,row.names(cells)),
                                   "x"]
sp2_clusters[rp_names,"y"] = cells[match(rp_names,row.names(cells)),
                                   "y"]

sp2_clusters$sample="sample2"

sp2_clusters[sp2_clusters$anno %in% c("B", "Plasma"),"anno"]="B_Cells" 
sp2_clusters[sp2_clusters$anno %in% c("Macrophage"),"anno"]="Macrophages" 
sp2_clusters[sp2_clusters$anno %in% c("T","NK"),"anno"]="T_Cells" 
sp2_clusters[sp2_clusters$anno %in% c("Fibroblast"),"anno"]="Stromal" 
sp2_clusters[sp2_clusters$anno %in% c("Mast"),"anno"]="Mast_Cells" 
sp2_clusters$anno = factor(sp2_clusters$anno,
                           levels=c("Tumor", "Stromal","Macrophages",
                                    "Myoepithelial", "T_Cells", "B_Cells",
                                    "Endothelial", "Dendritic", "Mast_Cells",
                                    "ML","LP","Adipocyte"))

sp2_clusters$cells=paste(sp2_clusters$cells,"-sp2", sep="")


sp2_clusters$anno = as.character(sp2_clusters$anno)
selected_sp2 = sp2_clusters[sp2_clusters$anno %in% target_clusters,]

selected_sp2$anno = factor(selected_sp2$anno, levels=target_clusters)
selected_sp2$cluster =paste("c",as.numeric(factor(selected_sp2$anno)), 
                            sep="")

selected_sp2 = selected_sp2[,c("cluster","x","y","cells","sample","anno")]
clusters = rbind(selected_sp1, selected_sp2)
table(clusters$sample, clusters$anno)
clusters$anno = factor(clusters$anno, target_clusters)

###############################################################################


shared_genes = intersect(row.names(sp1$cm), row.names(sp2$cm))
seed_number= 589

grid_length=40
keep_targets = c(shared_genes,
                 intersect(sp1$probe,sp2$probe),
                 intersect(sp1$codeword,sp2$codeword))

usage_sv= peakRAM({
sp1_sp2_vectors = get_vectors(x= list("sample1" = sp1$trans_info, 
                                      "sample2" = sp2$trans_info),
                              sample_names = c("sample1", "sample2"),
                              cluster_info = clusters, bin_type="square",
                              bin_param=c(grid_length,grid_length), 
                              test_genes = shared_genes, 
                              n_cores = 5)
})
sp1_sp2_nc_vectors = create_genesets(x= list("sample1" = sp1$trans_info, 
                                             "sample2" = sp2$trans_info),
                                     sample_names = c("sample1", "sample2"),
                                     name_lst=list(probe=intersect(sp1$probe,sp2$probe), 
                                                   codeword=intersect(sp1$codeword,sp2$codeword)),
                                     bin_type="square",
                                     bin_param=c(40,40), 
                                     cluster_info = NULL)

set.seed(seed_number)

usage_glm= peakRAM({
sp1_sp2_lasso_with_nc = lasso_markers(gene_mt=sp1_sp2_vectors$gene_mt,
                                      cluster_mt = sp1_sp2_vectors$cluster_mt,
                                      sample_names=c("sample1","sample2"),
                                      keep_positive=TRUE, 
                                      background=sp1_sp2_nc_vectors)
})
###############################################################################
# limma

cm_new1=sp1$cm[,setdiff(colnames(sp1$cm),c(sp1$zero_cells,"cate"))]
colnames(cm_new1)=paste(colnames(cm_new1),"_1", sep="")
biorep1=factor(rep(c("hb1"),c(ncol(cm_new1))))
names(biorep1) = colnames(cm_new1)

cm_sp2=sp2$cm[,setdiff(colnames(sp2$cm),c(sp2$zero_cells,"cate"))]
colnames(cm_sp2)=paste(colnames(cm_sp2),"-sp2", sep="")
biosp2=factor(rep(c("hb2"),c(ncol(cm_sp2))))
names(biosp2) = colnames(cm_sp2)

rep_clusters = clusters
y=DGEList(cbind(cm_new1[shared_genes,intersect(colnames(cm_new1),rep_clusters$cells)], cm_sp2[shared_genes,intersect(colnames(cm_sp2),rep_clusters$cells)]))
sample=c(biorep1[intersect(colnames(cm_new1),rep_clusters$cells)],biosp2[intersect(colnames(cm_sp2),rep_clusters$cells)])
logcounts.all <- normCounts(y,log=TRUE,prior.count=0.1)
all.ct <- factor(rep_clusters$cluster)

design <- model.matrix(~0+all.ct+sample)
colnames(design)[1:(length(levels(all.ct)))] <- levels(all.ct)

mycont <- matrix(0,ncol=length(levels(all.ct)),nrow=length(levels(all.ct)))
colnames(mycont)<-levels(all.ct)
diag(mycont)<-1
mycont[upper.tri(mycont)]<- -1/(length(levels(all.ct))-1)
mycont[lower.tri(mycont)]<- -1/(length(levels(all.ct))-1)

# Fill out remaining rows with 0s
zero.rows <- matrix(0,ncol=length(levels(all.ct)),nrow=(ncol(design)-length(levels(all.ct))))
test <- rbind(mycont,zero.rows)
usage_limma= peakRAM({
fit <- lmFit(logcounts.all,design)
fit.cont <- contrasts.fit(fit,contrasts=test)
fit.cont <- eBayes(fit.cont,trend=TRUE,robust=TRUE)

fit.cont$genes <- shared_genes
limma_dt<-decideTests(fit.cont)
})

###############################################################################
set.seed(1939)


cm_new1=sp1$cm[,setdiff(colnames(sp1$cm),c(sp1$zero_cells,"cate"))]
colnames(cm_new1)=paste(colnames(cm_new1),"_1", sep="")
biorep1=factor(rep(c("hb1"),c(ncol(cm_new1))))
names(biorep1) = colnames(cm_new1)

sp1_seu=CreateSeuratObject(counts = cm_new1, project = "hb1")
sp1_seu=AddMetaData(object=sp1_seu, metadata = biorep1, col.name="biorep")
DefaultAssay(sp1_seu) <- 'RNA'
sp1_seu = NormalizeData(sp1_seu,normalization.method ="LogNormalize")
sp1_seu = FindVariableFeatures(sp1_seu, selection.method = "vst",
                               nfeatures = nrow(sp1$cm), verbose = FALSE)

sp1_seu=ScaleData(sp1_seu)

sp1_seu=RunPCA(sp1_seu, npcs = 50, verbose = FALSE)

set.seed(1939)
cm_sp2=sp2$cm[,setdiff(colnames(sp2$cm),c(sp2$zero_cells,"cate"))]
colnames(cm_sp2)=paste(colnames(cm_sp2),"_2", sep="")
biosp2=factor(rep(c("hb2"),c(ncol(cm_sp2))))
names(biosp2) = colnames(cm_sp2)
sp2_seu=CreateSeuratObject(counts = cm_sp2, project = "hb2")
sp2_seu=AddMetaData(object=sp2_seu, metadata = biosp2, col.name="biorep")
DefaultAssay(sp2_seu) <- 'RNA'
sp2_seu <- NormalizeData(sp2_seu,normalization.method ="LogNormalize")
sp2_seu = FindVariableFeatures(sp2_seu, selection.method = "vst",
                               nfeatures = nrow(sp2$cm), verbose = FALSE)

sp2_seu <- ScaleData(sp2_seu)
sp2_seu=RunPCA(sp2_seu, npcs = 50, verbose = FALSE)

merged_seu = Merge_Seurat_List(
    list(sp1_seu[shared_genes,], sp2_seu[shared_genes,]),
    # add.cell.ids = c(),
    merge.data = FALSE,
    project = "hbreast"
)
Idents(merged_seu) = clusters$cluster

# re-join layers after integration
merged_seu[["RNA"]] <- JoinLayers(merged_seu[["RNA"]])

VariableFeatures(merged_seu) <- row.names(merged_seu)
merged_seu <- ScaleData(merged_seu, features = row.names(merged_seu),
                        verbose = FALSE)
set.seed(1939)
# print(ElbowPlot(merged_seu, ndims = 50))
merged_seu <- RunPCA(merged_seu, features = row.names(merged_seu),
                     npcs = 50, verbose = FALSE)
merged_seu <- RunUMAP(merged_seu, dims = 1:20, verbose = FALSE)


usage_fm= peakRAM({
FM_res <- FindAllMarkers(merged_seu, only.pos = TRUE,logfc.threshold = 0.1)
})
###############################################################################
results_df <- data.frame(
    technology = "xenium",
    dataset = "Xenium human breast cancer",
    genes_n =  nrow(sp1$cm),
    transcript_n =  nrow(sp1$trans_info)+ nrow(sp2$trans_info),
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


output_file_name <- "xenium_human_breast_cancer_5core.csv"


out_dir <- file.path(script_dir, "..", "..", 
                     "data/dataset_computational_complexity")
out_dir <- normalizePath(out_dir)

setwd(out_dir)


saveRDS(clusters,"xenium_hbreast_clusters.Rds")
saveRDS(sp1_sp2_lasso_with_nc,"xenium_hbreast_jazzPanda_res_lst.Rds")
saveRDS(fit.cont, "xenium_hbreast_fit_cont_obj.Rds")
saveRDS(FM_res, "xenium_hbreast_seu_markers.Rds")
saveRDS(merged_seu, "xenium_hbreast_seu.Rds")
saveRDS(sp1_sp2_vectors,"xenium_hbreast_sq40_vector_lst.Rds")
# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)
