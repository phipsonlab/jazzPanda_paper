
library(jazzPanda)
library(SpatialExperiment)
library(Seurat)
library(ggplot2)
library(data.table)
library(spatstat)
library(dplyr)
library(glmnet)
library(caret)
library(tidyr)

PEAK_MEM_MB =6

# A help function to load Xenium data
get_xenium_data<-function(path,mtx_name, trans_name="transcript_info.csv.gz",
                          cells_name="cell_info.csv.gz"){

  transcript_info <- as.data.frame(fread(paste(path, trans_name,sep="")))
  cell_info <- as.data.frame(fread(paste(path,cells_name,sep="")))

  data <- Read10X(data.dir = paste(path,mtx_name, sep=""))

  cm <- as.matrix(data$`Gene Expression`)
  r_codeword <- as.matrix(data$`Negative Control Codeword`)

  r_probe <- as.matrix(data$`Negative Control Probe`)
  # merge negative control genes and real genes
  cm_neg <- as.data.frame(rbind(r_probe, r_codeword))
  zero_cells <- colnames(cm)[colSums(cm)==0]

  transcript_info$x <- as.numeric(transcript_info$x_location)
  transcript_info$y <- as.numeric(transcript_info$y_location)

  cell_info$x <- as.numeric(cell_info$x_centroid)
  cell_info$y <- as.numeric(cell_info$y_centroid)

  return (list(cm = cm, cm_neg=cm_neg, zero_cells = zero_cells,
               trans_info=transcript_info, cell_info=cell_info,
               probe = row.names(r_probe),
               codeword=row.names(r_codeword)))

}

args_all   <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())  # interactive fallback
}

project_root <- normalizePath(dirname(dirname(dirname(script_dir))), 
                              mustWork = TRUE)

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
                           levels=c("Tumor", "Stromal","Macrophages","Myoepithelial",
                                    "T_Cells", "B_Cells","Endothelial", "Dendritic", "Mast_Cells","ML","LP","Adipocyte"))

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

shared_genes = intersect(row.names(sp1$cm), row.names(sp2$cm))
seed_number= 589
# Time and memory complexity vs number of tiles


args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
repeat_times <- as.integer(args[2])
tile_lens <- as.numeric(strsplit(args[3], ",")[[1]])

# total number of tiles
num_tiles <- length(tile_lens)

# Calculate the tile_len index for the current task
tile_len_index <- ((task_id - 1) %% num_tiles) + 1


# Extract the current tile length
curr_size <- tile_lens[tile_len_index]

cat("curr_size =",curr_size,"\n")
grid_length = curr_size

gc(reset = TRUE)
sv_st = Sys.time()
sv_imem <- sum(gc()[, PEAK_MEM_MB])

# get spatial vectors
sp1_sp2_vectors = get_vectors(x= list("sample1" = sp1$trans_info, 
                                      "sample2" = sp2$trans_info),
                              sample_names = c("sample1", "sample2"),
                              cluster_info = clusters, bin_type="square",
                              bin_param=c(grid_length,grid_length), 
                              test_genes = shared_genes , 
                              n_cores = 2)
sv_ed= Sys.time()


sv_fmem <-sum(gc()[, PEAK_MEM_MB])
sv_dmem <- sv_fmem - sv_imem
sv_time = difftime(sv_ed,sv_st, units = "mins")
cat("creating vectors for",nrow(sp1$cm),"genes and",
    length(unique(clusters$cluster)),"clusters took", sv_time,"mins and",
    round(sv_dmem / (1024),5),"GB memory")

# create noise vectors
sp1_sp2_nc_vectors =  create_genesets(x= list("sample1" = sp1$trans_info, 
                                              "sample2" = sp2$trans_info),
                                      sample_names = c("sample1", "sample2"),
                                      name_lst=list(probe=intersect(sp1$probe,sp2$probe), 
                                                    codeword=intersect(sp1$codeword,sp2$codeword)),
                                      bin_type="square",
                                      bin_param=c(grid_length,grid_length), 
                                      cluster_info = NULL)

set.seed(seed_number)
gc(reset=TRUE)
glm_st = Sys.time()
glm_imem <-sum(gc()[, PEAK_MEM_MB])
jazzPanda_res_lst = lasso_markers(gene_mt=sp1_sp2_vectors$gene_mt,
                                    cluster_mt = sp1_sp2_vectors$cluster_mt,
                                    sample_names=c("sample1","sample2"),
                                    keep_positive=TRUE, 
                                    background=sp1_sp2_nc_vectors)
glm_fmem <- sum(gc()[, PEAK_MEM_MB])
glm_ed <- Sys.time()
glm_dmem <- glm_fmem - glm_imem
glm_time <- difftime(glm_ed, glm_st, units = "mins")
cat("jazzPanda-linear model took", glm_time,"mins and",
    round(glm_dmem / (1024),5),"GB memory")


curr_top_res =jazzPanda_res_lst$top_result


saveRDS(sp1_sp2_vectors, paste(paste("xenium_hbreast_vectors_gr", curr_size, sep=""),".Rds", sep=""))
saveRDS(jazzPanda_res_lst, paste(paste("xenium_hbreast_jazzPanda_res_lst_gr", curr_size, sep=""),".Rds", sep=""))

rm(sp1_sp2_nc_vectors,jazzPanda_res_lst,sp1_sp2_vectors)

#Create a data frame to store the results
results_df <- data.frame(
    genes = shared_genes,
    task_id = rep(task_id, length(shared_genes)),
    curr_size = rep(curr_size, length(shared_genes)),
    top_cluster = curr_top_res[match(shared_genes,curr_top_res$gene),"top_cluster"],
    glm_coef = curr_top_res[match(shared_genes,curr_top_res$gene),"glm_coef"],
    seed_number=rep(589, length(shared_genes))
)
colnames(results_df)[4] <- paste("top_cluster_gr", curr_size, sep="")
colnames(results_df)[5] <- paste("glm_coef_gr", curr_size, sep="")

output_file_name <- sprintf("xenium_hbreast_ntiles_mg_gr%d_id%d.csv",curr_size,task_id)

write.csv(results_df, output_file_name, row.names = FALSE)
