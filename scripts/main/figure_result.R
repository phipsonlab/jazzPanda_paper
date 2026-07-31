################################################################################
# Figure 3: Cosmx human liver normal sample"

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

# keep cancer tissue only 
liver_normal = cellCoords[cellCoords$slide_ID_numeric==1 ,]




# load count matrix 
counts <-seu[["RNA"]]@counts
normal_cells = row.names(liver_normal)


counts_normal_sample = counts[, normal_cells]

rm(counts)
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

plot_data_sp(cluster_info=cluster_info,file_prefix = data_nm, out_dir = figure_result,
             my_colors = my_colors, ct_nm = "anno_name")


################################################################################
cl = "Stellate.cells"
p_cl<- ggplot(data = cluster_info[cluster_info$cluster==cl, ],
              aes(x = x, y = y))+ 
    facet_wrap(~cluster)+
    #geom_hex(bins = 120)+
    geom_point(size=0.01, alpha=0.5,color="steelblue")+
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


jpeg(file.path(figure_result, "figure_cosmx_hhliver_stellate_xy.jpg"),
     height = 1000, width = 1100, res=200)
p_cl 
dev.off()

plot_data_sp(cluster_info=cluster_info,file_prefix = data_nm, out_dir = figure_result,
             my_colors = my_colors, ct_nm = "anno_name")

plot_cluster_props(cluster_info=cluster_info,file_prefix = data_nm, 
                   out_dir = figure_result,
                   my_colors = my_colors, ct_nm = "anno_name")


plot_umap_seu(cluster_info=cluster_info,seu=seu, file_prefix = data_nm, 
              out_dir = figure_result,
              my_colors = my_colors, ct_nm = "anno_name", fig_w = 1400)


mg = "IGFBP7"
iters_sp1=hl_normal[hl_normal$feature_name %in%  mg,]
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


jpeg(file.path(figure_result, "figure_cosmx_hhliver_top1_xy_glm.jpg"),
     height = 1000, width = 1100, res=200)
p1

dev.off()

################################################################################

library(data.table)
library(Seurat)
library(ggplot2)
library(matrixStats)
library(patchwork)
library(reshape2)
library(tidyr)
library(dplyr)
library(here)
source(here("scripts/utils.R"))

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


data_nm  <- "xenium_hbreast"
# load geenrated data
cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))
colnames(cluster_info)[6] = "anno_name"

cluster_names = paste0("c", 1:9)
cluster_info$cluster = factor(cluster_info$cluster,
                              levels=cluster_names)
ct_names =c("Tumor", "Stromal","Macrophages","Myoepithelial", "T_Cells", 
            "B_Cells","Endothelial", "Dendritic", "Mast_Cells")
cluster_info$anno_name = factor(cluster_info$anno_name,
                                levels=ct_names)

cluster_info = cluster_info[order(cluster_info$anno_name), ]
cluster_info$cells = paste0("_", cluster_info$cells)
anno_df = unique(cluster_info[c("cluster", "anno_name")])
anno_df$anno_name = factor(anno_df$anno_name, levels = ct_names)
jazzPanda_res_lst = readRDS(here(data_path,paste0(data_nm, "_jazzPanda_res_lst.Rds")))

fit.cont = readRDS(here(data_path,paste0(data_nm, "_fit_cont_obj.Rds")))

FM_result= readRDS(here(data_path,paste0(data_nm, "_seu_markers.Rds")))
sv_lst = readRDS(here(data_path,paste0(data_nm, "_sq40_vector_lst.Rds")))
seu = readRDS(here(data_path,paste0(data_nm, "_seu.Rds")))
seu <- subset(seu, cells = cluster_info$cells)
nbins = 1600
Idents(seu)=cluster_info$anno_name[match(colnames(seu), cluster_info$cells)]
seu$sample = cluster_info$sample[match(colnames(seu), cluster_info$cells)]

xhb_color<- c("#FC8D62","#66C2A5" ,"#8DA0CB","#E78AC3",
              "#A6D854","skyblue","purple3","#E5C498","blue")
################################################################################


plot_data_sp(cluster_info=cluster_info,file_prefix = data_nm, 
             out_dir = figure_result,
             my_colors = xhb_color, ct_nm = "anno_name",
             reverse_y = TRUE)


plot_cluster_props(cluster_info=cluster_info,file_prefix = data_nm, 
                   out_dir = figure_result,
                   my_colors = xhb_color, ct_nm = "anno_name" )

################################################################################

## spatial coordinates for one cluster

ct_nm = "T_Cells"
p_cl<- ggplot(data = cluster_info[cluster_info$anno_name==ct_nm, ],
              aes(x = x, y = y))+ 
    facet_wrap(~sample, nrow=2)+
    #geom_hex(bins = 120)+
    geom_point(size=0.01, alpha=0.7,color="#A6D854")+
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


jpeg(file.path(figure_result, "figure4c_xenium_hbreast_tcell_xy.jpg"),
     width = 1000, height = 1800, res=200)
p_cl 
dev.off()

################################################################################
## spatial coordinates for top marker gene
mk_gene = "IL7R"

iters_sp1= sp1$trans_info$feature_name %in% mk_gene
vis_r1 =sp1$trans_info[iters_sp1,
                       c("x","y","feature_name")]
vis_r1$sample="sample1"
iters_rep2= sp2$trans_info$feature_name %in% mk_gene
vis_r2 =sp2$trans_info[iters_rep2,
                       c("x","y","feature_name")]
vis_r2$sample="sample2"
p1<- ggplot(data = vis_r1,
            aes(x = x, y = y))+ 
    #geom_hex(bins = tr_hex[cl_ids])+
    geom_point(size=0.01, alpha=0.7)+
    facet_wrap(~sample, scales="free", ncol=1)+
    #scale_fill_gradient(low="white", high="maroon4") + 
    guides(fill = guide_colorbar(height= unit(0.6, "cm")))+
    defined_theme+  theme(legend.position = "bottom",
                          strip.background = element_rect(colour = "white", 
                                                          fill = "white", 
                                                          linetype = "solid"), 
                          panel.border = element_rect(colour = "NA", fill="NA"),
                          strip.text = element_text(size = 15),
                          panel.spacing = unit(0, "lines"))

p2<- ggplot(data = vis_r2,
            aes(x = x, y = y))+ 
    #geom_hex(bins = tr_hex[cl_ids])+
    geom_point(size=0.01, alpha=0.7)+
    facet_wrap(~sample, scales="free", ncol=1)+
    #scale_fill_gradient(low="white", high="maroon4") + 
    guides(fill = guide_colorbar(height= unit(0.6, "cm")))+
    defined_theme+  theme(legend.position = "bottom",
                          strip.background = element_rect(colour = "white", 
                                                          fill = "white", 
                                                          linetype = "solid"), 
                          panel.border = element_rect(colour = "NA", fill="NA"),
                          strip.text = element_text(size = 15),
                          panel.spacing = unit(0, "lines"))

lyt = (p1 / p2) 
layout_design <- lyt + patchwork::plot_layout(heights = c(1,1),widths = c(1, 1)) 
jpeg(file.path(figure_result, "figure4d_xenium_hbreast_top1_xy_glm.jpg"),
     width = 1000, height = 1800, res=200)
layout_design
dev.off()

################################################################################
cluster_nm = anno_df[anno_df$anno_name==ct_nm, "cluster"]
dff = as.data.frame(cbind(sv_lst$cluster_mt[,cluster_nm],sv_lst$gene_mt[,mk_gene]))
colnames(dff) = c("cluster", mk_gene)
dff$sample= "sample1"
dff[nbins:(nbins*2),"sample"] = "sample2"
dff$vector_id = c(1:nbins, 1:nbins)
long_df <- dff %>% 
    pivot_longer(cols = -c(cluster, sample, vector_id), names_to = "gene", 
                 values_to = "vector_count")
long_df$gene = factor(long_df$gene, levels=mk_gene)
dff$sample = factor(dff$sample , levels=c("sample1","sample2"))
p<-ggplot(long_df, aes(x = cluster, y = vector_count, color =gene )) +
    geom_point(color="black", size=0.01) +
    facet_wrap(~sample, nrow=2) +
    labs(x = paste("T cell", " cluster vector ", sep=""), y = "IL7R gene vector") +
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

pdf(file.path(figure_result,"figure4e_xenium_hbreast_vvplot.pdf"),
    height = 9, width = 5)
p
dev.off()



