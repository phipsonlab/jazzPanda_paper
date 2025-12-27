
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(reshape2)
library(jazzPanda)
library(tidyr)
library(dplyr)
library(corrplot)
library(edgeR)
library(ComplexUpset)
library(dplyr)


# load data
seu <- readRDS("/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/CosMx_normal_and_diseased_liver_samples/LiverDataReleaseSeurat_newUMAP.RDS")
# local cell metadata
metadata <- as.data.frame(seu@meta.data)
cat("Number of cells in cancer and normal samples: ",nrow(metadata))
## remove low quality cells
qc_cols <- c("qcFlagsRNACounts", "qcFlagsCellCounts", "qcFlagsCellPropNeg",
             "qcFlagsCellComplex", "qcFlagsCellArea","qcFlagsFOV")
metadata <- metadata[!apply(metadata[, qc_cols], 1, function(x) any(x == "Fail")), ]
cat("Total number of cells remaining after removing low quality cells: ", nrow(metadata))


cell_info_cols = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm",
                   "nCount_negprobes","nFeature_negprobes","nCount_falsecode",
                   "nFeature_falsecode","slide_ID_numeric", "Run_Tissue_name",
                   "fov","cellType","niche", "cell_id")
cellCoords <- metadata[, cell_info_cols]
px_to_mm <- function(data){
    all_fv = unique(data$fov)
    parm_df = as.data.frame(matrix(0, ncol=5, nrow=length(all_fv)))
    colnames(parm_df) = c("fov","y_slope","y_intcp","x_slope","x_intcp")
    parm_df$fov = all_fv
    for (fv in all_fv){
        curr_fov = data[data$fov == fv, ]
        curr_fov = curr_fov[order(curr_fov$x_FOV_px), ]
        curr_fov = curr_fov[c(1, nrow(curr_fov)), ]
        curr_fov = curr_fov[c(1,2),c("y_slide_mm","x_slide_mm","y_FOV_px","x_FOV_px") ]
        # mm to px for y
        y_slope = (curr_fov[2,"y_slide_mm"] - curr_fov[1,"y_slide_mm"]) / (curr_fov[2,"y_FOV_px"] - curr_fov[1,"y_FOV_px"])
        y_intcp = curr_fov[2,"y_slide_mm"]- (y_slope*curr_fov[2,"y_FOV_px"])
        # mm to px for x
        x_slope = (curr_fov[2,"x_slide_mm"] - curr_fov[1,"x_slide_mm"]) / (curr_fov[2,"x_FOV_px"] - curr_fov[1,"x_FOV_px"])
        x_intcp = curr_fov[2,"x_slide_mm"]- (x_slope*curr_fov[2,"x_FOV_px"])
        parm_df[parm_df$fov==fv,"y_slope"] = y_slope
        parm_df[parm_df$fov==fv,"y_intcp"] = y_intcp
        parm_df[parm_df$fov==fv,"x_slope"] = x_slope
        parm_df[parm_df$fov==fv,"x_intcp"] = x_intcp
    }
    return (parm_df)
    
}

# number of cells per fov
# all fovs contain at least 2 cells
fov_summary = as.data.frame(table(cellCoords[cellCoords$slide_ID_numeric==2,"fov"]))

# keep cancer tissue only
liver_cancer = cellCoords[cellCoords$slide_ID_numeric==2 ,]

# caculate the slope and intercept parameters for each fov
parm_df = px_to_mm(liver_cancer)

# convert px to mm for each cell based on the calculated params
liver_cancer <- liver_cancer %>%
    left_join(parm_df, by = 'fov') %>%
    mutate(
        x_mm = x_FOV_px * x_slope + x_intcp,
        y_mm = y_FOV_px * y_slope + y_intcp
    ) %>%
    select(-x_slope, -y_slope, -x_intcp, -y_intcp)


## load transcript coordinates
gc(reset = TRUE)
load_tr_imem <- sum(gc()[, "(Mb)"])
transcriptCoords <-seu@misc$transcriptCoords

all_transcripts_cancer <- transcriptCoords[transcriptCoords$slideID == 2,]
# remove transctipt dfetections for low quality cells
all_transcripts_cancer <- all_transcripts_cancer[all_transcripts_cancer$cell_id %in% liver_cancer$cell_id, ]
rm(transcriptCoords)

load_tr_fmem <- sum(gc()[, "(Mb)"])
load_tr_dmem <- load_tr_fmem - load_tr_imem
ncells_tr = length(unique(all_transcripts_cancer$cell_id))
cat("loading transcript coordinates data for",1000,
    "genes and", ncells_tr,"cells took",
    round(load_tr_dmem / (1024),5), "GB memory")

all_transcripts_cancer <- all_transcripts_cancer %>%
    left_join(parm_df, by = 'fov') %>%
    mutate(
        x_mm = x_FOV_px * x_slope + x_intcp,
        y_mm = y_FOV_px * y_slope + y_intcp
    ) %>%
    select(-x_slope, -y_slope, -x_intcp, -y_intcp)

all_transcripts_cancer$x = all_transcripts_cancer$x_mm
all_transcripts_cancer$y = all_transcripts_cancer$y_mm
all_transcripts_cancer$feature_name = all_transcripts_cancer$target
# table(is.na(all_transcripts_cancer$x))

# refine clustering
selected_cols = c("x_FOV_px", "y_FOV_px","x_slide_mm", "y_slide_mm", "slide_ID_numeric", "Run_Tissue_name", "fov","cellType","niche" )
clusters_info = cellCoords[cellCoords$slide_ID_numeric=="2" & (row.names(cellCoords) %in% row.names(metadata)),selected_cols ]

colnames(clusters_info) =c("x_FOV_px","y_FOV_px" , "x", "y", "slide_ID_numeric", "Run_Tissue_name", "fov","cellTypes","niche")
clusters_info$cluster = as.character(clusters_info$cellTypes)
clusters_info[clusters_info$cellTypes %in% c("Antibody.secreting.B.cells", "Mature.B.cells"),"cluster"] = "B"
clusters_info[clusters_info$cellTypes %in% c("CD3+.alpha.beta.T.cells", "gamma.delta.T.cells.1"),"cluster"] = "T"
clusters_info[clusters_info$cellTypes %in% c("Non.inflammatory.macrophages", "Inflammatory.macrophages"),"cluster"] = "Macrophages"
# clusters_info[clusters_info$cellTypes %in% c("Hep.3","Hep.4", "Hep.5","Hep.6"),"cluster"] = "Hep3"
ig_clusters = c("NotDet")

clusters_info = clusters_info[clusters_info$cluster != "NotDet",]
colnames(clusters_info)[10] = "anno"
clusters_info$anno = factor(clusters_info$anno,
                            levels=c("tumor_1","tumor_2","Macrophages","T",
                                     "Periportal.LSECs","Stellate.cells","B",
                                     "Central.venous.LSECs","Cholangiocytes","Hep",
                                     "Portal.endothelial.cells",
                                     "NK.like.cells","Erthyroid.cells"))
clusters_info$cluster =  paste("c",as.numeric(factor(clusters_info$anno)), sep="")
clusters_info$cluster =factor(clusters_info$cluster, levels=paste("c",1:13, sep=""))
clusters_info$cell_id = row.names(clusters_info)

clusters_info$sample = "cancer"
# dim = 332873     11
clusters_info = clusters_info[!duplicated(clusters_info[,c("x","y")]),]
clusters_info$x = clusters_info$x * 1000
clusters_info$y = clusters_info$y * 1000


hl_cancer = all_transcripts_cancer[,c("x","y","feature_name")]
all_genes = row.names(seu[["RNA"]]@counts)
rm(all_transcripts_cancer, seu)
hl_cancer$x = hl_cancer$x * 1000
hl_cancer$y = hl_cancer$y * 1000

## load negative control coordinates
### negprobes

negprobes_coords <- hl_cancer[grepl("^NegPrb", hl_cancer$feature_name), ]

### falsecode
falsecode_coords <- hl_cancer[grepl("^FalseCode", hl_cancer$feature_name), ]



falsecode_coords$sample = "cancer"
negprobes_coords$sample = "cancer"
kpt_cols = c("x","y","feature_name","sample")
nc_dff = as.data.frame(rbind(falsecode_coords[,kpt_cols],
                                 negprobes_coords[,kpt_cols]))
falsecode_names = unique(falsecode_coords$feature_name)
negprobe_names = unique(negprobes_coords$feature_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
cluster_g1 <- as.integer(args[2])
tile_lens <- as.numeric(strsplit(args[3], ",")[[1]])

# total number of tiles
num_tiles <- length(tile_lens)

# Calculate the tile_len index for the current task
tile_len_index <- ((task_id - 1) %% num_tiles) + 1


# Extract the current tile length
curr_size <- tile_lens[tile_len_index]

cat("curr_size =",curr_size,"\n")
grid_length = curr_size


perm_p = compute_permp(x= list("cancer" = hl_cancer),
                       cluster_info=clusters_info, 
                       perm.size=5000,
                       bin_type="square",
                       bin_param=c(grid_length,grid_length),
                       test_genes= all_genes,
                       correlation_method = "pearson", 
                       n_cores = 5,
                       correction_method="BH")


saveRDS(perm_p, paste(paste("cosmx_hliver_cancer_perm_lst_gr", 
                            curr_size, sep=""),".Rds", sep=""))
rm(perm_p)
###############################################################################
gc(reset = TRUE)
sv_st = Sys.time()
sv_imem <- sum(gc()[, "(Mb)"])
hliver_vector_lst = get_vectors(x = list(cancer = hl_cancer),
                                sample_names = "cancer",
                                cluster_info = clusters_info,
                                bin_type = "square",
                                bin_param=c(grid_length,grid_length),
                                test_genes = all_genes,n_cores = 2)
sv_ed= Sys.time()
sv_fmem <-sum(gc()[, "(Mb)"])
sv_dmem <- sv_fmem - sv_imem
sv_time = difftime(sv_ed,sv_st, units = "mins")

cat("creating vectors for",length(all_genes),"genes and",
    length(unique(clusters_info$cluster)),"clusters with",
    grid_length,"x",grid_length, "tiles took", sv_time,"mins and",
    round(sv_dmem / (1024),5),"GB memory")

nc_vectors = create_genesets(x=list("cancer" = nc_dff),
                             name_lst=list(falsecode=falsecode_names,
                                           negprobe=negprobe_names),
                             sample_names = c("cancer"),
                             cluster_info=NULL,
                             bin_type="square",
                             bin_param=c(grid_length, grid_length))


set.seed(989)
gc(reset=TRUE)
glm_st = Sys.time()
glm_imem <-sum(gc()[, "(Mb)"])
jazzPanda_res_lst = lasso_markers(gene_mt=hliver_vector_lst$gene_mt,
                                  cluster_mt = hliver_vector_lst$cluster_mt,
                                  sample_names=c("cancer"),
                                  keep_positive=TRUE,
                                  background=nc_vectors,
                                  n_fold = 10)

glm_fmem <- sum(gc()[, "(Mb)"])
glm_ed <- Sys.time()
glm_dmem <- glm_fmem - glm_imem
glm_time <- difftime(glm_ed, glm_st, units = "mins")
cat("jazzPanda-linear model took", glm_time,"mins and",
    round(glm_dmem / (1024),5),"GB memory \n")

curr_top_res =jazzPanda_res_lst$top_result


saveRDS(hliver_vector_lst, paste(paste("cosmx_hliver_cancer_vectors_gr", curr_size, sep=""),".Rds", sep=""))
saveRDS(jazzPanda_res_lst, paste(paste("cosmx_hliver_cancer_jazzPanda_res_lst_gr", curr_size, sep=""),".Rds", sep=""))

rm(nc_vectors,jazzPanda_res_lst,hliver_vector_lst)

#Create a data frame to store the results
glm_results_df <- data.frame(
    genes = all_genes,
    task_id = rep(task_id, length(all_genes)),
    curr_size = rep(curr_size, length(all_genes)),
    top_cluster = curr_top_res[match(all_genes,curr_top_res$gene),"top_cluster"],
    glm_coef = curr_top_res[match(all_genes,curr_top_res$gene),"glm_coef"],
    seed_number=rep(989, length(all_genes))
)
colnames(glm_results_df)[4] <- paste("top_cluster_gr", curr_size, sep="")
colnames(glm_results_df)[5] <- paste("glm_coef_gr", curr_size, sep="")

output_file_name <- sprintf("cosmx_hlc_glm_ntiles_mg_gr%d_id%d.csv",curr_size,task_id)
write.csv(glm_results_df, output_file_name, row.names = FALSE)
