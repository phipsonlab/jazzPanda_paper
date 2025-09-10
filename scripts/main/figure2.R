################################################################################
# Figure 2: Cosmx human liver normal sample"

# w_x[2] - w_x[1] --> 10726
# w_y[2] - w_y[1] --> 8142

library(Seurat)
library(ggplot2)
library(matrixStats)
library(patchwork)
library(reshape2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/utils.R"))


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
                   "fov","cellType","niche","cell_id")
cellCoords <- metadata[, cell_info_cols]
gc(reset = TRUE)
load_cell_imem <- sum(gc()[, PEAK_MEM_MB])
# keep cancer tissue only 
liver_normal = cellCoords[cellCoords$slide_ID_numeric==1 ,]
load_cell_fmem <- sum(gc()[, PEAK_MEM_MB])
load_cell_dmem <- load_cell_fmem - load_cell_imem

cat("loading cell coordinates data for",
    nrow(cellCoords[cellCoords$Run_Tissue_name=="NormalLiver",]),"cells took",
    round(load_cell_dmem / (1024),5),"GB memory")


# load count matrix 
counts <-seu[["RNA"]]@counts
normal_cells = row.names(liver_normal)

gc(reset = TRUE)
load_cm_imem <- sum(gc()[, PEAK_MEM_MB])

counts_normal_sample = counts[, normal_cells]

load_cm_fmem <- sum(gc()[, PEAK_MEM_MB])
load_cm_dmem <- load_cm_fmem - load_cm_imem
rm(counts)
cat("loading count matrix data for",nrow(counts_normal_sample),"genes and",
    ncol(counts_normal_sample),"cells took",round(load_cm_dmem / (1024),5),"GB memory")

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

gc(reset = TRUE)
load_tr_imem <- sum(gc()[, PEAK_MEM_MB])

all_transcripts_normal <- transcriptCoords[transcriptCoords$slideID == 1,]
# remove 76 transcripts from fov21
all_transcripts_normal <- all_transcripts_normal[all_transcripts_normal$cell_id %in% liver_normal$cell_id, ]

load_tr_fmem <- sum(gc()[, PEAK_MEM_MB])
load_tr_dmem <- load_tr_fmem - load_tr_imem
rm(transcriptCoords)
cat("loading transcript coordinates data for",nrow(counts_normal_sample),"genes and",
    ncol(counts_normal_sample),"cells took",round(load_tr_dmem / (1024),5),"GB memory")

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

rm(all_transcripts_normal)
hl_normal$x = hl_normal$x * 1000
hl_normal$y = hl_normal$y * 1000


################################################################################
data_nm  <- "cosmx_hhliver"
# load generated data
cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))
cluster_info$anno_name = cluster_info$cluster

cluster_names = c( "B", "Central.venous.LSECs", "Cholangiocytes", "Erthyroid.cells", 
                   "Hep.1", "Hep.3", "Hep.4", "Hep.5", "Hep.6", "Macrophages",
                   "NK.like.cells", "Periportal.LSECs", "Portal.endothelial.cells", 
                   "Stellate.cells", "T")  

cluster_info$cluster = factor(cluster_info$cluster,
                              levels=cluster_names)

cluster_info$anno_name = factor(cluster_info$anno_name,
                                levels=cluster_names)
cluster_info = cluster_info[order(cluster_info$anno_name), ]


fit.cont = readRDS(here(data_path,paste0(data_nm, "_fit_cont_obj.Rds")))

perm_lst = readRDS(here(data_path,paste0(data_nm, "_perm_lst.Rds")))


FM_result= readRDS(here(data_path,paste0(data_nm, "_seu_markers.Rds")))
sv_lst = readRDS(here(data_path,paste0(data_nm, "_sq40_vector_lst.Rds")))
seu = readRDS(here(data_path,paste0(data_nm, "_seu.Rds")))
seu <- subset(seu, cells = cluster_info$cell_id)
nbins = 1600
Idents(seu)=cluster_info$anno_name[match(colnames(seu), cluster_info$cell_id)]
seu$sample = cluster_info$sample[match(colnames(seu), cluster_info$cell_id)]

################################################################################
# Figure 2(c) UMAP
plot_umap_seu(cluster_info=cluster_info,seu=seu, file_prefix = data_nm, 
              out_dir = fig2,
              my_colors = my_colors, ct_nm = "anno_name", fig_w = 1400)


################################################################################
# Figure 2(a)
plot_cluster_props(cluster_info=cluster_info,file_prefix = data_nm, 
                   out_dir = fig2,
                   my_colors = my_colors, ct_nm = "anno_name")


################################################################################
# Figure 2(b)
plot_data_sp(cluster_info=cluster_info,file_prefix = data_nm, out_dir = fig2,
             my_colors = my_colors, ct_nm = "anno_name")


################################################################################
# Figure 2(d) 
## spatial coordinates for one cluster
cl = "Stellate.cells"
p_cl<- ggplot(data = cluster_info[cluster_info$cluster==cl, ],
              aes(x = x, y = y))+ 
    facet_wrap(~cluster)+
    #geom_hex(bins = 120)+
    geom_point(size=0.01, alpha=0.1,color="#f58231")+
    guides(fill = guide_colorbar(height= unit(0.6, "cm")))+
    scale_fill_gradient(low="white", high="black") + 
    #scale_fill_gradient(low="white", high="maroon4") + 
    defined_theme+ 
    theme(legend.position = "bottom",
          strip.background = element_rect(colour = "white", 
                                          fill = "white", 
                                          linetype = "solid"), 
          panel.border = element_rect(colour = "NA", fill="NA"),
          strip.text = element_text(size = 15),
          panel.spacing = unit(0, "lines"))


jpeg(file.path(fig2, "figure2d_cosmx_hhliver_stellate_xy.jpg"),
     height = 1000, width = 1100, res=200)
p_cl 
dev.off()


################################################################################
# Figure 2(e) 
## spatial coordinates for top marker gene

perm_res = get_perm_adjp(perm_lst)
obs_corr = get_cor(perm_lst)
head(perm_res)

obs_cutoff = quantile(obs_corr[, cl], 0.75)
perm_cl=intersect(row.names(perm_res[perm_res[,cl]<0.05,]),
                  row.names(obs_corr[obs_corr[, cl]>obs_cutoff,]))
inters=perm_cl
rounded_val=signif(as.numeric(obs_corr[inters,cl]), digits = 3)
inters_df = as.data.frame(cbind(gene=inters, value=rounded_val))
inters_df$value = as.numeric(inters_df$value)
inters_df=inters_df[order(inters_df$value, decreasing = TRUE),]
inters = inters_df$gene[1]
iters_sp1=hl_normal[hl_normal$feature_name %in%  inters,]
vis_r1 =iters_sp1[,c("x","y","feature_name")]
p1<- ggplot(data = vis_r1,
            aes(x = x, y = y))+ 
    #geom_hex(bins = tr_hex[cl_ids])+
    geom_point(size=0.01, alpha=0.1)+
    facet_wrap(~feature_name, scales="free", ncol=3)+
    scale_fill_gradient(low="white", high="maroon4") + 
    guides(fill = guide_colorbar(height= unit(0.6, "cm")))+
    defined_theme+  theme(legend.position = "bottom",
                          strip.background = element_rect(colour = "white", 
                                                          fill = "white", 
                                                          linetype = "solid"), 
                          panel.border = element_rect(colour = "NA", fill="NA"),
                          strip.text = element_text(size = 15),
                          panel.spacing = unit(0, "lines"))


jpeg(file.path(fig2, "figure2e_cosmx_hhliver_top1_xy_perm.jpg"),
     height = 1000, width = 1100, res=200)
p1

dev.off()

################################################################################
# Figure 2(f) 
mk_gene ="IGFBP7"
perm_res[mk_gene,cl]
obs_corr[mk_gene,cl]
dff = as.data.frame(cbind(sv_lst$cluster_mt[,cl],
                          sv_lst$gene_mt[,mk_gene]))
colnames(dff) = c("cluster", mk_gene)
dff$vector_id = c(1:nbins)
long_df <- dff %>% 
    pivot_longer(cols = -c(cluster, vector_id), names_to = "gene", 
                 values_to = "vector_count")
long_df$gene = factor(long_df$gene, levels=mk_gene)

p<-ggplot(long_df, aes(x = cluster, y = vector_count, color =gene )) +
    geom_point(color="black", size=0.01) +
    facet_wrap(~gene, scales = "free_y", nrow=1) +
    labs(x = paste(cl, " cluster vector ", sep=""), y = "IGFBP7 gene vector") +
    theme_minimal()+
    scale_y_continuous(expand = c(0.01,0.01))+ 
    scale_x_continuous(expand =  c(0.01,0.01))+ 
    theme(panel.grid = element_blank(),legend.position = "none",
          strip.text = element_blank(),#element_text(size = rel(1)),
          axis.line=element_blank(),
          axis.text=element_text(size = 12),
          axis.ticks=element_line(color="black"),
          axis.title=element_text(size = 15),
          panel.border  =element_rect(colour = "black", fill=NA, linewidth=0.5)
    )

pdf(file.path(fig2, "figure2f_cosmx_hhliver_top1_perm_vvplot.pdf"),
    height = 5, width = 5)
p
dev.off()


