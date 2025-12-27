################################################################################
# Figure 4: Xenium human breast cancer samples
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
# Figure 4(a) 

plot_data_sp(cluster_info=cluster_info,file_prefix = data_nm, 
             out_dir = fig4,
             my_colors = xhb_color, ct_nm = "anno_name",
             reverse_y = TRUE)

# Figure 4(b) 
plot_cluster_props(cluster_info=cluster_info,file_prefix = data_nm, 
                   out_dir = fig4,
                   my_colors = xhb_color, ct_nm = "anno_name" )

################################################################################
# Figure 4(d) 
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


jpeg(file.path(fig4, "figure4c_xenium_hbreast_tcell_xy.jpg"),
     width = 1000, height = 1800, res=200)
p_cl 
dev.off()

################################################################################
# Figure 4(e) 
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
jpeg(file.path(fig4, "figure4d_xenium_hbreast_top1_xy_glm.jpg"),
     width = 1000, height = 1800, res=200)
layout_design
dev.off()

################################################################################
# Figure 2(f) 
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

pdf(file.path(fig4,"figure4e_xenium_hbreast_vvplot.pdf"),
    height = 9, width = 5)
p
dev.off()

