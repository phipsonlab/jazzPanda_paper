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
library(data.table)
library(limma)
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

dist_df$sample = factor(dist_df$sample, levels = selected_samples,
                        labels = c("Xenium breast cancer", 
                                   "CosMx human liver", 
                                   "MERSCOPE breast cancer"))
p_dist = ggplot(dist_df, aes(factor(value_grp), prop)) +
    geom_col(fill = "grey35", colour = NA) +
    geom_text(aes(label = percent(prop, accuracy = 0.1)), vjust = -0.3, size = 2) +
    scale_x_discrete(labels = lab_fun) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ sample) +
    labs(x = "Gene count", y = "Proportion of gene count", fill="gene count") +
    theme_bw()+
    theme(strip.background = element_rect(fill=NA, color="black"))


ggsave(file.path(figure_intro,"fig_count_distribution.png"), p_dist, width = 10, height = 3,  dpi = 200)


data_nm="cosmx_hhliver"
cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))


if (data_nm=="xenium_hbreast"){
    ct_names =c("Tumor", "Stromal","Macrophages","Myoepithelial", "T_Cells", 
                "B_Cells","Endothelial", "Dendritic", "Mast_Cells")
    colnames(cluster_info)[6] = "anno_name"
    cluster_names = paste0("c", 1:9)
    cluster_info$cluster = factor(cluster_info$cluster,
                                  levels=cluster_names)
}else if(data_nm == "cosmx_hhliver"){
    ct_names = c( "B", "Central.venous.LSECs", "Cholangiocytes", "Erthyroid.cells", 
                       "Hep.1", "Hep.3", "Hep.4", "Hep.5", "Hep.6", "Macrophages",
                       "NK.like.cells", "Periportal.LSECs", "Portal.endothelial.cells", 
                       "Stellate.cells", "T")  
    cluster_info$anno_name = cluster_info$cluster
    
}

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
FM_result = FM_result[FM_result$p_val_adj<0.05,]
total_genes = length(rownames(seu))     
dt = decideTests(fit.cont)
total_genes = length(rownames(seu))

limma_sig = data.frame(
    cluster = colnames(dt),
    n_sig   = colSums(dt == 1)
) %>%
    mutate(prop = n_sig / total_genes)
limma_sig$method = "t-test"

wrst_sig = data.frame(
    table(FM_result$cluster)
) 
colnames(wrst_sig) = c("cluster", "n_sig")
wrst_sig$method = "WRST"
wrst_sig$prop = wrst_sig$n_sig / total_genes

sig_df = rbind(limma_sig, wrst_sig)
sig_df$ct = anno_df$anno_name[match(sig_df$cluster, anno_df$cluster)]
sig_df$ct = factor(sig_df$ct, levels = ct_names)

p_bar_1 = ggplot(sig_df, aes(ct, prop, fill = method)) +
    geom_col(position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent,
                       expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = c("t-test" = "#0072B2", "WRST" = "#E69F00"))+   # blue / orange
    labs(x = "Cluster", y = "Proportion of significant genes", fill = NULL,title = "CosMx human liver") +
    theme_bw() +
    theme( plot.title    = element_text(size = 15, hjust=0.5,margin = margin(t = 2, b = 8)),
        legend.position          = "inside",
        legend.position.inside   = c(0.99, 0.99),
        legend.justification     = c(1, 1),
        legend.title = element_blank(),
        axis.title.x  = element_blank(),
        legend.background        = element_rect(fill = alpha("white", 0.85), colour = "black"),
        legend.key.size          = unit(0.5, "cm"),
        legend.margin            = margin(2, 4, 2, 4),
        axis.text.x              = element_text(angle = 45, hjust = 1)
    )


library(limma)
res <- decideTests(fit.cont, adjust.method = "BH", p.value = 0.05, lfc = 0)

gene_tab <- data.frame(
    gene          = rownames(res),
    cluster_count = rowSums(res==1)
)
gene_tab <- gene_tab[gene_tab$cluster_count > 0, ]   # drop never-significant genes

gene_limma <- data.frame(table(gene_tab$cluster_count))
colnames(gene_limma) = c("count", "freq")
gene_limma$method = "t-test"

gene_tab=data.frame(table(FM_result$gene))
colnames(gene_tab) = c("gene", "cluster_count")
gene_fm=data.frame(table(gene_tab$cluster_count))
colnames(gene_fm) = c("count", "freq")
gene_fm$method = "WRST"

gene_df = rbind(gene_limma, gene_fm)

p_bar_2 = ggplot(gene_df, aes(x=count, freq, fill=method)) +
     geom_col(position = position_dodge(width = 0.9)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    scale_fill_manual(values = c("t-test" = "#0072B2", "WRST" = "#E69F00"))+  
    labs(x = "Number of clusters a gene is a marker for", y = "Number of genes",title = "CosMx human liver") +
    theme_bw() +
    theme( plot.title    = element_text(size = 15, hjust=0.5,
                                        margin = margin(t = 2, b = 8)),
        legend.position          = "inside",
        legend.title = element_blank(),
        legend.position.inside   = c(0.99, 0.99),
        legend.justification     = c(1, 1),
        legend.background        = element_rect(fill = alpha("white", 0.85), colour = "black"),
        legend.key.size          = unit(0.4, "cm"),
        legend.margin            = margin(2, 4, 2, 4),
        axis.text.x              = element_text( hjust =0.5)
    )

p_combined <- (p_bar_1 | free(p_bar_2, type = "label", side = "b")) +
    plot_layout(widths = c(0.6, 0.4))
pdf(file.path(figure_intro, paste0(data_nm, "_mg_prop_ncmg.pdf")),
    width=12, height=5)
print(p_combined)
dev.off()
