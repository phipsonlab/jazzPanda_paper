library(matrixStats)
library(Seurat)
library(reshape2)
library(jazzPanda)
library(tidyr)
library(dplyr)

# calculate memory
PEAK_MEM_ID=which(colnames(gc())=="max used")
PEAK_MEM_MB = PEAK_MEM_ID+1

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



# refine clustering
selected_cols = c("x_FOV_px", "y_FOV_px","x_slide_mm", "y_slide_mm", "slide_ID_numeric", "Run_Tissue_name", "fov","cellType","niche" )
clusters_info = cellCoords[cellCoords$slide_ID_numeric=="2" & (row.names(cellCoords) %in% row.names(metadata)),selected_cols ]

colnames(clusters_info) =c("x_FOV_px","y_FOV_px" , "x", "y", "slide_ID_numeric", "Run_Tissue_name", "fov","cellTypes","niche")
clusters_info$cluster = as.character(clusters_info$cellTypes)
clusters_info[clusters_info$cellTypes %in% c("Antibody.secreting.B.cells", "Mature.B.cells"),"cluster"] = "B"
clusters_info[clusters_info$cellTypes %in% c("CD3+.alpha.beta.T.cells", "gamma.delta.T.cells.1"),"cluster"] = "T"
clusters_info[clusters_info$cellTypes %in% c("Non.inflammatory.macrophages", "Inflammatory.macrophages"),"cluster"] = "Macrophages"
ig_clusters = c("NotDet")

clusters_info = clusters_info[clusters_info$cluster != "NotDet",]
clusters_info$cluster = factor(clusters_info$cluster)
clusters_info$cell_id = row.names(clusters_info)
clusters_info$sample = "cancer"
# dim = 332873     11
clusters_info = clusters_info[!duplicated(clusters_info[,c("x","y")]),]
clusters_info$x = clusters_info$x * 1000
clusters_info$y = clusters_info$y * 1000


hl_cancer = all_transcripts_cancer[,c("x","y","feature_name")]
all_genes = row.names(seu[["RNA"]]@counts)
n_transcripts = nrow(all_transcripts_cancer[all_transcripts_cancer$feature_name %in% all_genes, ])
rm(all_transcripts_cancer, seu)

hl_cancer$x = hl_cancer$x * 1000
hl_cancer$y = hl_cancer$y * 1000



################ Time and memory complexity vs number of cores ################
cat("\n Complexity measure starts \n")


args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
repeat_times <- as.integer(args[2])
n_core_lst <- as.numeric(strsplit(args[3], ",")[[1]])

# total number of tested cores
num_test_cores <- length(n_core_lst)

# Calculate the n_core index for the current task
n_core_id <- ((task_id - 1) %% num_test_cores) + 1

# Calculate the repetition number for the current task
repetition_number <- floor((task_id - 1) / num_test_cores) + 1

n_core = n_core_lst[n_core_id]
n_gene = 1000
cat("\n task_id =", task_id, ", repeat_times =", repeat_times,", n_core_lst =", n_core_lst, "\n")
cat("\n Repeat =", repetition_number, ", n_core =", n_core,", n_gene =", n_gene, "\n")

# Initialize memory and time tracking
cat("\n")
gc(reset = TRUE)
sv_st = Sys.time()
sv_imem <- sum(gc()[, PEAK_MEM_MB])
grid_length = 40
# creat sv based on n_genes only
hliver_vector_lst = get_vectors(x = list(cancer = hl_cancer),
                                sample_names = "cancer",
                                cluster_info = NULL,
                                bin_type = "square",
                                bin_param=c(grid_length,grid_length),
                                test_genes = all_genes,
                                n_cores=n_core)
sv_ed= Sys.time()
sv_fmem <-sum(gc()[, PEAK_MEM_MB])
sv_dmem <- sv_fmem-sv_imem
sv_time = difftime(sv_ed,sv_st, units = "mins")

#Create a data frame to store the results
results_df <- data.frame(
    task_id = task_id,
    rd = repetition_number,
    ncore = n_core,
    gene_n = n_gene,
    grid_length = grid_length,
    time_min_sv = sv_time,
    memory_MB_sv = sv_dmem,
    transcript_n = n_transcripts
)



args_all   <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())  # interactive fallback
}


## Output folder name
output_dir_nm <- "ncores_result"

## Path to /scripts/main/cosmx_hlc_simulation_result
out_dir <- file.path(script_dir, output_dir_nm)

setwd(out_dir)

# Output file name uniquely identified by SLURM_ARRAY_TASK_ID
output_file_name <- sprintf("complexity_ncores%d_id%d.csv", n_core,task_id)

# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)
