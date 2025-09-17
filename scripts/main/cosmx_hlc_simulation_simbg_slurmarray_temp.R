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
library(ComplexUpset)
library(spatstat)
library(data.table)
library(purrr)

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

colnames(clusters_info) =c("x_FOV_px","y_FOV_px" , "x", "y", 
                           "slide_ID_numeric", "Run_Tissue_name", "fov",
                           "cellTypes","niche")
clusters_info$cluster = as.character(clusters_info$cellTypes)
clusters_info[clusters_info$cellTypes %in% c("Antibody.secreting.B.cells", "Mature.B.cells"),"cluster"] = "B"
clusters_info[clusters_info$cellTypes %in% c("CD3+.alpha.beta.T.cells", "gamma.delta.T.cells.1"),"cluster"] = "T"
clusters_info[clusters_info$cellTypes %in% c("Non.inflammatory.macrophages", "Inflammatory.macrophages"),"cluster"] = "Macrophages"
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
# dim = 460164 12
# duplicated 3169 rows

#counts_cancer_sample = counts_cancer_sample[, clusters_info$cell_id]

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

# template starts here

grid_length= 40
sim_length = 100

################ negative control background
falsecode_coords$sample = "cancer"
negprobes_coords$sample = "cancer"
kpt_cols = c("x","y","feature_name","sample")
nc_dff = as.data.frame(rbind(falsecode_coords[,kpt_cols],
                                 negprobes_coords[,kpt_cols]))
falsecode_names = unique(falsecode_coords$feature_name)
negprobe_names = unique(negprobes_coords$feature_name)

# prepare background #detection for simulationn
n_tr_bg = as.data.frame(nc_dff %>% dplyr::group_by(feature_name) %>% dplyr::count())

total_rows_bg <- sum(n_tr_bg$n)
falsecode_names_sim =paste("sim", falsecode_names, sep="_")
negprobe_names_sim =paste("sim", negprobe_names, sep="_")

sim_length_bg =length(falsecode_names) + length(negprobe_names)

################
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
x_rng <- range(hl_cancer$x, na.rm = TRUE)
y_rng <- range(hl_cancer$y, na.rm = TRUE)
w_x <- c(floor(x_rng[1]), ceiling(x_rng[2]))
w_y <- c(floor(y_rng[1]), ceiling(y_rng[2]))
W <- owin(w_x, w_y)


ntr_summary = as.data.frame(hl_cancer %>% dplyr::group_by(feature_name) %>% dplyr::count())

all_genes <- intersect(all_genes,ntr_summary$feature_name)
set.seed(989)
global_seed_lst = sample(1:999999, size=100,replace = FALSE)
global_seed_sub = global_seed_lst[((task_id-1)*10+1):(task_id*10)]
args <- commandArgs(trailingOnly = FALSE)
script_dir <- if (any(grepl("^--file=", args))) {
    dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
} else {
    normalizePath(getwd())   # fallback if run interactively
}

## Output folder name
output_dir_nm <- "cosmx_hlc_simulation_result"

## Path to /scripts/main/cosmx_hlc_simulation_result
out_dir <- file.path(script_dir, output_dir_nm)

setwd(out_dir)

generate_uniform_points <- function(n_points, W, abs_jitter = 0.1) {
    stopifnot(is.numeric(W$xrange), length(W$xrange) == 2,
              is.numeric(W$yrange), length(W$yrange) == 2)
    side_n <- ceiling(sqrt(n_points))
    
    # Create base grid
    x_grid <- seq(W$xrange[1], W$xrange[2], length.out = side_n)
    y_grid <- seq(W$yrange[1], W$yrange[2], length.out = side_n)
    
    grid <- expand.grid(
        x = x_grid,
        y = y_grid
    )[1:n_points, ]  # Truncate to exact number of points
    
    
    grid$x <- grid$x + runif(n_points, -abs_jitter, abs_jitter)
    grid$y <- grid$y + runif(n_points, -abs_jitter, abs_jitter)
    
    # Clip to window
    grid$x <- pmin(pmax(grid$x, W$xrange[1]), W$xrange[2])
    grid$y <- pmin(pmax(grid$y, W$yrange[1]), W$yrange[2])
    
    # Assert output size
    stopifnot(nrow(grid) == n_points)
    
    return(grid)
}
for (gs in 1:length(global_seed_sub)){
  #  999488
  global_seed = global_seed_sub[gs]
  # keep coordinates of simulation
  simulated_tr = as.data.frame(matrix(0, ncol=4))
  colnames(simulated_tr) = c("x","y","feature_name","seed")

  sim_file_nm = paste("cosmx_hlc_simulation_s",global_seed, sep="")
  set.seed(global_seed)
  seed_lst = sample(1:999999, size=sim_length,replace = FALSE)
  # selected 100 real genes
  selected_genes = sample(all_genes,size = sim_length, replace = FALSE)

  # the number of transcripts of the selected 100 real genes
  n_tr <- ntr_summary %>%
      filter(feature_name %in% selected_genes) %>%
      slice(match(selected_genes, feature_name))
  
  simulated_list <- map2(seq_along(selected_genes), seed_lst, function(i, seed) {
      set.seed(seed)
      n <- n_tr$n[i]
      
      # Generate uniform grid points instead of random
      pts <- generate_uniform_points(n, W)
      
      data.frame(
          x = pts$x,
          y = pts$y,
          feature_name = paste0("sim_", n_tr$feature_name[i]),
          sim_genetr_seed = seed,
          stringsAsFactors = FALSE
      )
  })
  
  # Combine all simulated points into one data frame
  simulated_tr <- rbindlist(simulated_list)
  
  # Get unique simulated gene names
  sim_gene_names <- unique(simulated_tr$feature_name)

  all_vectors = get_vectors(x= list("cancer"=simulated_tr),
                            sample_names = "cancer",
                            cluster_info = clusters_info,
                            bin_type="square",
                            bin_param=c(grid_length,grid_length),
                            test_genes =sim_gene_names)
  saveRDS(all_vectors, 
          paste(sim_file_nm,"spatial_vectors.Rds", sep="_"))
  
  set.seed(global_seed)
  seed_lst_bg <- sample(1:999999, size = sim_length_bg, replace = FALSE)
    
  simulated_list_bg <- map2(
      seed_lst_bg,
      seq_len(sim_length_bg),
      function(seed, i) {
          set.seed(seed)
          n_points <- n_tr_bg$n[i]
          
          pts <- generate_uniform_points(n_points, W)
          
          data.frame(
              x = pts$x,
              y = pts$y,
              feature_name = paste0("sim_", n_tr_bg$feature_name[i]),
              sim_bctr_seed = seed,
              stringsAsFactors = FALSE
          )
      }
  )
  # Combine into one data frame
  simulated_tr_bg <- rbindlist(simulated_list_bg)
  
 
  nc_vectors_sim = create_genesets(x=list("cancer" = simulated_tr_bg),
                                   sample_names = c("cancer"),
                               name_lst=list(falsecode=falsecode_names_sim,
                                             negprobe=negprobe_names_sim),
                               bin_type="square",
                               cluster_info=NULL,
                               bin_param=c(grid_length, grid_length))

  saveRDS(nc_vectors_sim, 
          paste(sim_file_nm,"nc_spatial_vectors.Rds", sep="_"))
  
  gc(reset=TRUE)
  glm_st = Sys.time()
  jazzPanda_res_lst = lasso_markers(gene_mt=all_vectors$gene_mt,
                                      cluster_mt = all_vectors$cluster_mt,
                                      sample_names=c("cancer"),
                                      keep_positive=TRUE,
                                      background=nc_vectors_sim)

  saveRDS(jazzPanda_res_lst, 
          paste(sim_file_nm,"lasso_res_lst.Rds", sep="_"))
  

  set.seed(global_seed)
  perm_p = compute_permp(x=list("cancer"=simulated_tr),
                         cluster_info=clusters_info,
                         perm.size=1000,
                         bin_type="square",
                         bin_param=c(grid_length,grid_length),
                         test_genes=sim_gene_names,
                         correlation_method = "pearson",
                         correction_method="BH")

  saveRDS(perm_p, 
          paste(sim_file_nm,"perm_res_lst.Rds", sep="_"))

}

