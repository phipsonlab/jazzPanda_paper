library(ggplot2)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(jazzPanda)
library(here)
source(here("scripts/utils.R"))
library(scales)
library(ComplexUpset)
library(cowplot)

selected_samples = c("xenium_hbreast_sp1","cosmx_healthy", "merscope_sp1")
df = fread(here::here("data","gene_count_frequency.csv"))
df = as.data.frame(df)
df = df[df$sample %in% c(selected_samples),]
df$value = as.integer(df$value)


cap = 10
lab_fun = function(v){
    v = as.integer(as.character(v))
    ifelse(v >= cap, paste0(cap, "+"), as.character(v))
}

dist_df = df %>%
    mutate(value_grp = ifelse(value >= cap, cap, value)) %>%
    group_by(sample, value_grp) %>%
    summarise(prop = sum(prop), .groups = "drop")

dist_df$sample = factor(dist_df$sample, levels = selected_samples)
p_dist = ggplot(dist_df, aes(factor(value_grp), prop)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = percent(prop, accuracy = 0.1)), vjust = -0.3, size = 2) +
    scale_x_discrete(labels = lab_fun) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ sample) +
    labs(x = "Gene count", y = "Proportion of gene count", fill="gene count") +
    theme_bw()+
    theme(strip.background = element_rect(fill=NA, color="black"))


ggsave(file.path(figure_intro,"fig_count_distribution.png"), p_dist, width = 10, height = 3,  dpi = 200)


data_nm  <- "xenium_hbreast"
cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))
colnames(cluster_info)[6] = "anno_name"
cluster_names = paste0("c", 1:9)
cluster_info$cluster = factor(cluster_info$cluster,
                              levels=cluster_names)
ct_names =c("Tumor", "Stromal","Macrophages","Myoepithelial", "T_Cells", 
            "B_Cells","Endothelial", "Dendritic", "Mast_Cells")
cluster_info$anno_name = factor(cluster_info$anno_name,
                                levels=ct_names)
rownames(cluster_info) <-cluster_info$cells
cluster_info = cluster_info[order(cluster_info$anno_name), ]
anno_df = unique(cluster_info[c("cluster", "anno_name")])
anno_df$anno_name = factor(anno_df$anno_name, levels = ct_names)
jazzPanda_res_lst = readRDS(here(data_path,paste0(data_nm, "_jazzPanda_res_lst.Rds")))

fit.cont = readRDS(here(data_path,paste0(data_nm, "_fit_cont_obj.Rds")))
seu = readRDS(here(data_path,paste0(data_nm, "_seu.Rds")))
FM_result= readRDS(here(data_path,paste0(data_nm, "_seu_markers.Rds")))

total_genes = length(rownames(seu))     
dt = decideTests(fit.cont)
total_genes = length(rownames(seu))

sig_df = data.frame(
    cluster = colnames(dt),
    n_sig   = colSums(dt == 1)
) %>%
    mutate(prop = n_sig / total_genes)

sig_df$ct = anno_df$anno_name[match(sig_df$cluster, anno_df$cluster)]
ct_size = sig_df$ct[order(sig_df$n_sig, decreasing = TRUE)]
# reorder the x-axis factor by cluster size
sig_df$ct = factor(sig_df$ct, levels = ct_size)
library(patchwork)


fill_col = "grey30"

p_bar = ggplot(sig_df, aes(ct, prop)) +
    geom_col(fill = fill_col) +
    geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
              vjust = -0.3, size = 2.5) +
    scale_y_continuous(limits = c(0, 0.5), labels = scales::percent,
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Cluster", y = "Proportion of significant genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
mat = as.data.frame((dt == 1) * 1)
colnames(mat) = anno_df$anno_name[match(colnames(mat), anno_df$cluster)]

set_queries = lapply(ct_names, function(s) {
    ComplexUpset::upset_query(set = s, fill = fill_col)
})

p_upset = ComplexUpset::upset(mat, intersect = ct_names,
                              min_size = 2,          # fewer intersections -> narrower
                              queries = set_queries,
                              wrap=TRUE, keep_empty_groups= FALSE, name="",
                              stripes='white',
                              sort_intersections_by ="cardinality",  sort_sets= "ascending",min_degree=1,
                              set_sizes =( 
                                  upset_set_size()
                                  + theme(axis.title= element_blank(),
                                          axis.ticks.y = element_blank(),
                                          axis.text.y = element_blank())),
                              sort_intersections= "descending", warn_when_converting=FALSE,
                              warn_when_dropping_groups=TRUE,encode_sets=TRUE,
                              
                              width_ratio=0.3, height_ratio=1/4)+
                                theme(plot.title = element_text( size=15))
jpeg(file.path(figure_intro, paste0(data_nm, "_limma_mg_prop_intersection.jpg")),
     width=4000, height=2000, res=300)
print(p_bar + p_upset + plot_layout(widths = c(0.8, 3)))
dev.off()

