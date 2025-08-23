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

seu <- readRDS("/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/CosMx_normal_and_diseased_liver_samples/LiverDataReleaseSeurat_newUMAP.RDS")


# local cell metadata   
metadata <- as.data.frame(seu@meta.data)

## remove low quality cells
qc_cols <- c("qcFlagsRNACounts", "qcFlagsCellCounts", "qcFlagsCellPropNeg",
             "qcFlagsCellComplex", "qcFlagsCellArea","qcFlagsFOV")
metadata <- metadata[!apply(metadata[, qc_cols], 1, function(x) any(x == "Fail")), ]

cell_info_cols = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm",
                   "nCount_negprobes","nFeature_negprobes","nCount_falsecode",
                   "nFeature_falsecode","slide_ID_numeric", "Run_Tissue_name",
                   "fov","cellType","niche","cell_id")
cellCoords <- metadata[, cell_info_cols]

# keep cancer tissue only 
liver_normal = cellCoords[cellCoords$slide_ID_numeric==1 ,]

# load count matrix 
counts <-seu[["RNA"]]@counts
normal_cells = row.names(liver_normal)

counts_normal_sample = counts[, normal_cells]

dim(counts_normal_sample)



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
# fov 21 contains 1 cell only 
fov_summary = as.data.frame(table(cellCoords[cellCoords$slide_ID_numeric==1,"fov"]))


# caculate the slope and intercept parameters for each fov 
parm_df = px_to_mm(liver_normal)

# convert px to mm for eeach cell based on the calculated params 
liver_normal <- liver_normal %>%
    left_join(parm_df, by = 'fov') %>%
    mutate(
        x_mm = x_FOV_px * x_slope + x_intcp,
        y_mm = y_FOV_px * y_slope + y_intcp
    ) %>%
    select(-x_slope, -y_slope, -x_intcp, -y_intcp) 


transcriptCoords <-seu@misc$transcriptCoords


all_transcripts_normal <- transcriptCoords[transcriptCoords$slideID == 1,]
# remove 76 transcripts from fov21
all_transcripts_normal <- all_transcripts_normal[all_transcripts_normal$cell_id %in% liver_normal$cell_id, ]

#all_transcripts_normal = all_transcripts_normal[all_transcripts_normal$fov!=21,]

rm(transcriptCoords)

all_transcripts_normal <- all_transcripts_normal %>%
    left_join(parm_df, by = 'fov') %>%
    mutate(
        x_mm = x_FOV_px * x_slope + x_intcp,
        y_mm = y_FOV_px * y_slope + y_intcp
    ) %>%
    select(-x_slope, -y_slope, -x_intcp, -y_intcp) 

all_transcripts_normal$x = all_transcripts_normal$x_mm
all_transcripts_normal$y = all_transcripts_normal$y_mm
all_transcripts_normal$feature_name = all_transcripts_normal$target

hl_normal = all_transcripts_normal[,c("x","y","feature_name")]
all_genes = row.names(seu[["RNA"]]@counts)

rm(all_transcripts_normal, seu)
hl_normal$x = hl_normal$x * 1000
hl_normal$y = hl_normal$y * 1000

negprobes_coords <- hl_normal[grepl("^NegPrb", hl_normal$feature_name), ]

falsecode_coords <- hl_normal[grepl("^FalseCode", hl_normal$feature_name), ]

selected_cols = c("x_FOV_px", "y_FOV_px","x_slide_mm", "y_slide_mm", "slide_ID_numeric", "Run_Tissue_name", "fov","cellType","niche","cell_id" )
clusters_info = liver_normal[,selected_cols ]
# [1] "x_FOV_px"         "y_FOV_px"         "x_slide_mm"       "y_slide_mm"       "slide_ID_numeric"
# [6] "Run_Tissue_name"  "fov"   
colnames(clusters_info) =c("x_FOV_px","y_FOV_px" , "x", "y", "slide_ID_numeric", "Run_Tissue_name", "fov","cellTypes","niche","cell_id")
clusters_info$cluster = as.character(clusters_info$cellTypes)
clusters_info[clusters_info$cellTypes %in% c("Antibody.secreting.B.cells", "Mature.B.cells"),"cluster"] = "B"
clusters_info[clusters_info$cellTypes %in% c("CD3+.alpha.beta.T.cells", "gamma.delta.T.cells.1"),"cluster"] = "T"
clusters_info[clusters_info$cellTypes %in% c("Non.inflammatory.macrophages", "Inflammatory.macrophages"),"cluster"] = "Macrophages"
# clusters_info[clusters_info$cellTypes %in% c("Hep.3","Hep.4", "Hep.5","Hep.6"),"cluster"] = "Hep3"
ig_clusters = c("NotDet")

clusters_info = clusters_info[clusters_info$cluster != "NotDet",]
clusters_info$cluster = factor(clusters_info$cluster)
clusters_info$sample = "normal"
# dim = 332873     11
clusters_info = clusters_info[!duplicated(clusters_info[,c("x","y")]),]
clusters_info$x = clusters_info$x * 1000
clusters_info$y = clusters_info$y * 1000

###############################################################################

seed_number=589
set.seed(seed_number)
usage_sv = peakRAM({
hliver_vector_lst = get_vectors(x=list("normal"=hl_normal),
                                sample_names = "normal",
                                cluster_info = clusters_info,
                                bin_type="square",
                                bin_param=c(40,40), 
                                test_genes = row.names(counts_normal_sample),
                                n_cores = 5)

})

###############################################################################

grid_length = 40
seed_number=589
set.seed(seed_number)


usage_perm= peakRAM({
perm_p = compute_permp(x=list("normal"=hl_normal),
                       cluster_info=clusters_info, 
                       perm.size=5000,
                       bin_type="square",
                       bin_param=c(grid_length,grid_length),
                       test_genes=row.names(counts_normal_sample), 
                       correlation_method = "spearman", 
                       n_cores=5, 
                       correction_method="BH")

})


falsecode_coords$sample = "normal"
negprobes_coords$sample = "normal"
kpt_cols = c("x","y","feature_name","sample")
nc_coords = rbind(falsecode_coords[,kpt_cols],negprobes_coords[,kpt_cols])


falsecode_names = unique(falsecode_coords$feature_name)
negprobe_names = unique(negprobes_coords$feature_name)
nc_vectors = create_genesets(x=list("normal"=nc_coords),
                             sample_names = "normal",
                             name_lst=list(falsecode=falsecode_names,
                                           negprobe=negprobe_names),
                             bin_type="square",
                             bin_param=c(grid_length, grid_length),
                             cluster_info=NULL)
set.seed(seed_number)
###############################################################################

usage_glm = peakRAM({
jazzPanda_res_lst = lasso_markers(gene_mt=hliver_vector_lst$gene_mt,
                                  cluster_mt = hliver_vector_lst$cluster_mt,
                                  sample_names=c("normal"),
                                  keep_positive=TRUE, 
                                  background=nc_vectors,
                                  n_fold = 10)

})

###############################################################################

y <- DGEList(counts_normal_sample[, clusters_info$cell_id])
y$genes <-row.names(counts_normal_sample)

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
usage_limma = peakRAM({
fit <- lmFit(logcounts,design)
fit.cont <- contrasts.fit(fit,contrasts=mycont)
fit.cont <- eBayes(fit.cont,trend=TRUE,robust=TRUE)

#treat_res <- treat(fit.cont,lfc=0.5)
limma_dt<-decideTests(fit.cont)
})

###############################################################################
hln_seu=CreateSeuratObject(counts = counts_normal_sample[, clusters_info$cell_id], project = "hl_normal")
Idents(hln_seu) = clusters_info[match(colnames(hln_seu), 
                                      clusters_info$cell_id),"cluster"]
hln_seu = NormalizeData(hln_seu, verbose = FALSE,
                        normalization.method = "LogNormalize")
hln_seu = FindVariableFeatures(hln_seu, selection.method = "vst", 
                               nfeatures = 1000, verbose = FALSE)

hln_seu=ScaleData(hln_seu, verbose = FALSE)
usage_fm= peakRAM({
seu_markers <- FindAllMarkers(hln_seu, only.pos = TRUE,logfc.threshold = 0.25)

})

#Create a data frame to store the results
results_df <- data.frame(
    technology = "cosmx",
    dataset = "CosMx healthy human liver",
    genes_n = nrow(counts_normal_sample),
    transcript_n =  nrow(hl_normal),
    cells_n = ncol(counts_normal_sample),
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


output_file_name <- "cosmx_human_healthy_liver_5core.csv"

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


saveRDS(clusters_info,"cosmx_hhliver_clusters.Rds")
saveRDS(jazzPanda_res_lst,"cosmx_hhliver_jazzPanda_res_lst.Rds")
saveRDS(fit.cont, "cosmx_hhliver_fit_cont_obj.Rds")
saveRDS(seu_markers, "cosmx_hhliver_seu_markers.Rds")
saveRDS(hln_seu, "cosmx_hhliver_seu.Rds")
saveRDS(perm_p,"cosmx_hhliver_perm_lst.Rds")
saveRDS(hliver_vector_lst,"cosmx_hhliver_sq40_vector_lst.Rds")

write.csv(results_df, output_file_name, row.names = FALSE)
