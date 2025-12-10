################################################################################
# Figure 5: Compare marker genes with FindMarkers and limma 
library(jazzPanda)
library(SpatialExperiment)
library(Seurat)
library(ggplot2)
library(data.table)
library(spatstat)
library(dplyr)
library(glmnet)
library(ComplexUpset)
library(patchwork)
library(RColorBrewer)
library(ggvenn)
library(tidyr)
library(xtable)
library(limma)
library(here)
source(here("scripts/utils.R"))
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir  <- dirname(normalizePath(script_path))


################################################################################
# One sample: CosMx healthy liver
chln_jazzPanda_res_lst = readRDS(file.path(here("data//dataset_computational_complexity/cosmx_hhliver_jazzPanda_res_lst.Rds")))
chln_jazzPanda_res = get_top_mg(chln_jazzPanda_res_lst,coef_cutoff=0.2)  
chln_perm_lst = readRDS(file.path(here("data/dataset_computational_complexity/cosmx_hhliver_perm_lst.Rds")))
chln_perm_res = chln_perm_lst$perm.pval.adj
chln_obs_corr = chln_perm_lst$obs.stat
hln_seu=readRDS(file.path(here("data/dataset_computational_complexity/cosmx_hhliver_seu.Rds")))
chln_fit = readRDS(file.path(here("data/dataset_computational_complexity/cosmx_hhliver_fit_cont_obj.Rds")))
chln_vecs= readRDS(file.path(here("data/dataset_computational_complexity/cosmx_hhliver_sq40_vector_lst.Rds")))
hln_seu_markers=readRDS(file.path(here("data/dataset_computational_complexity/cosmx_hhliver_seu_markers.Rds")))

limma_dt<-decideTests(chln_fit) 
plot_lst=list()
for (cl in c("Hep.4", "Stellate.cells", "Cholangiocytes")){
    obs_cutoff = quantile(chln_obs_corr[, cl], 0.75)
    findM_sig =hln_seu_markers[hln_seu_markers$cluster==cl & hln_seu_markers$p_val_adj<0.05,"gene"]
    limma_sig=row.names(limma_dt[limma_dt[,cl]==1,])
    perm_cl=intersect(row.names(chln_perm_res[chln_perm_res[,cl]<0.05,]),
                      row.names(chln_obs_corr[chln_obs_corr[, cl]>obs_cutoff,]))
    lasso_cl=chln_jazzPanda_res[chln_jazzPanda_res$top_cluster==cl, "gene"]
    df_mt =as.data.frame(matrix(FALSE,nrow=nrow(chln_perm_res),ncol=4))
    row.names(df_mt) =row.names(chln_perm_res)
    colnames(df_mt)=c("jazzPanda-glm",
                      "jazzPanda-correlation",
                      "Wilcoxon Rank Sum Test","limma")
    df_mt[findM_sig,"Wilcoxon Rank Sum Test"] = TRUE
    df_mt[limma_sig,"limma"] = TRUE
    df_mt[perm_cl,"jazzPanda-correlation"] = TRUE
    
    df_mt[lasso_cl,"jazzPanda-glm"] = TRUE
    
    cl = sub(cl, pattern = ".cells", replacement="")
    p = plot(upset(df_mt,
               intersect=c("Wilcoxon Rank Sum Test", "limma",
                           "jazzPanda-correlation","jazzPanda-glm"),
               wrap=TRUE, keep_empty_groups= FALSE, name="",
               #themes=theme_grey(),
               stripes='white',
               #stripes=upset_stripes(geom=geom_segment(size=8),colors=c('grey95', 'grey95', 'grey95',"grey95")),
               sort_intersections_by ="cardinality", sort_sets= FALSE,min_degree=1,
               set_sizes =( 
                   upset_set_size()
                   + theme(axis.title= element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.text.y = element_blank())),
               sort_intersections= "descending", warn_when_converting=FALSE,
               warn_when_dropping_groups=TRUE,encode_sets=TRUE,
               width_ratio=0.3, height_ratio=1/4)+
             ggtitle(paste(cl,"cells"))+
                 theme(plot.title = element_text( size=20))

    )
    plot_lst[[cl]] = p
    dev.off()
}
combined_plot <- wrap_plots(plot_lst, ncol = 3)
pdf(file.path(fig5, "cosmx_hhliver_upset.pdf"),height=6, width=18)
print(combined_plot)
dev.off()
################################################################################
cor_M = cor(chln_vecs$gene_mt,
            chln_vecs$cluster_mt,method = "spearman")
# cor_M[cor_M <= 0.7] = 0
# cor_M[cor_M > 0.7]=1
plot_lst=list()

for (cl in c("Hep.4", "Stellate.cells", "Cholangiocytes")){
    obs_cutoff = quantile(chln_obs_corr[, cl], 0.75)
    fm_cl=FindMarkers(hln_seu, ident.1 = cl, only.pos = TRUE,
                      logfc.threshold = 0.1)
    fm_cl = fm_cl[fm_cl$p_val_adj<0.05, ]
    fm_cl = fm_cl[order(fm_cl$avg_log2FC, decreasing = TRUE),]
    to_plot_fm =row.names(fm_cl)
    # FM_pt = data.frame("name"=to_plot,"rk"= 1:length(to_plot),
    #                    "y"= cumsum(cor_M[to_plot, cl]),
    #                    "type"="Wilcoxon Rank Sum Test")
    limma_cl=limma_dt[limma_dt[,cl]==1,]
    
    limma_cl<-topTable(chln_fit,coef=cl,p.value = 0.05, n=Inf, sort.by = "p")
    limma_cl = limma_cl[limma_cl$logFC>0, ]
    to_plot_lm = row.names(limma_cl)
    
    perm_cl=intersect(row.names(chln_perm_res[chln_perm_res[,cl]<0.05,]),
                      row.names(chln_obs_corr[chln_obs_corr[, cl]>obs_cutoff,]))
    rounded_val=signif(as.numeric(chln_obs_corr[perm_cl,cl]), digits = 3)
    roudned_pval=signif(as.numeric(chln_perm_res[perm_cl,cl]), digits = 3)
    perm_sorted = as.data.frame(cbind(gene=perm_cl, value=rounded_val, pval=roudned_pval))
    perm_sorted$value = as.numeric(perm_sorted$value)
    perm_sorted$pval = as.numeric(perm_sorted$pval)
    perm_sorted=perm_sorted[order(-perm_sorted$pval,
                                  perm_sorted$value, 
                                  decreasing = TRUE),]
    
    
    lasso_sig = chln_jazzPanda_res[chln_jazzPanda_res$top_cluster==cl,]
    lasso_sig = lasso_sig[order(lasso_sig$glm_coef, decreasing = TRUE),]
    # lasso_pt = data.frame("name"=lasso_sig$gene,"rk"= 1:nrow(lasso_sig),
    #                       "y"= cumsum(cor_M[lasso_sig$gene, cl]),
    #                       "type"="jazzPanda-glm")
    # 
    # limma_pt = data.frame("name"=to_plot_lm,"rk"= 1:length(to_plot_lm),
    #                       "y"= cumsum(cor_M[to_plot_lm, cl]),
    #                       "type"="limma")
    # 
    # corr_pt = data.frame("name"=perm_sorted$gene,"rk"= 1:length(perm_sorted$gene),
    #                      "y"= cumsum(cor_M[perm_sorted$gene, cl]),
    #                      "type"="jazzPanda-correlation")

    FM_pt = data.frame("name"=to_plot_fm,"rk"= 1:length(to_plot_fm),
                       "y"= get_cmr_ma(to_plot_fm,cor_M = cor_M, cl = cl),
                       "type"="Wilcoxon Rank Sum Test")
    
    
    lasso_pt = data.frame("name"=lasso_sig$gene,"rk"= 1:nrow(lasso_sig),
                          "y"= get_cmr_ma(lasso_sig$gene,cor_M = cor_M, cl = cl),
                          "type"="jazzPanda-glm")

    limma_pt = data.frame("name"=to_plot_lm,"rk"= 1:length(to_plot_lm),
                          "y"= get_cmr_ma(to_plot_lm,cor_M = cor_M, cl = cl),
                          "type"="limma")

    corr_pt = data.frame("name"=perm_sorted$gene,"rk"= 1:length(perm_sorted$gene),
                         "y"= get_cmr_ma(perm_sorted$gene,cor_M = cor_M, cl = cl),
                         "type"="jazzPanda-correlation")

    cl = sub(cl, pattern = ".cells", replacement="")
    #data_lst = rbind(limma_pt, lasso_pt,FM_pt)
    data_lst = rbind(limma_pt, FM_pt,corr_pt,lasso_pt)
    data_lst$type = factor(data_lst$type)
    data_lst$rk = as.numeric(data_lst$rk)
    p <-ggplot(data_lst, aes(x = rk, y = y, color = type)) +
        geom_step(size = 0.8) +  # type = "s"
        scale_color_manual(values = c("jazzPanda-correlation" = "orange",
                                      "jazzPanda-glm" = "red",
                                      "limma" = "black",
                                      "Wilcoxon Rank Sum Test" = "blue"

        )) +
        scale_x_continuous(limits = c(0, 50))+
        labs(title = paste(cl, "cells"), x = "Rank of marker genes",
             y = "Cumulative average correlation",
             color = NULL) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, size=16),
              axis.line = element_blank(),  
              legend.text  = element_text(size=13),
              panel.border = element_rect(color = "black", 
                                          fill = NA, linewidth = 1),
              legend.position = "inside",
              legend.position.inside = c(0.98, 0.02),
              legend.justification = c("right", "bottom"),
              legend.background = element_rect(color = "black", 
                                               fill = "white", linewidth = 0.5),
              legend.box.background = element_rect(color = "black",
                                                   linewidth = 0.5)
        )
    plot_lst[[cl]] <- p
}
combined_plot <- wrap_plots(plot_lst, ncol = 3)
pdf(file.path(fig5, "cosmx_hhliver_rank_ma.pdf"), height=5, width=15)
print(combined_plot)
dev.off()

# Two samples: Xenium human breast 
################################################################################
xhb_glm_lst = readRDS(file.path(here("data/dataset_computational_complexity/xenium_hbreast_jazzPanda_res_lst.Rds")))
xhb_jazzPanda_res = get_top_mg(xhb_glm_lst,coef_cutoff=0.2)  
xhb_fm=readRDS(file.path(here("data/dataset_computational_complexity/xenium_hbreast_seu_markers.Rds")))
xhb_seu=readRDS(file.path(here("data/dataset_computational_complexity/xenium_hbreast_seu.Rds")))
xhb_fit = readRDS(file.path(here("data/dataset_computational_complexity/xenium_hbreast_fit_cont_obj.Rds")))
xhb_vecs= readRDS(file.path(here("data/dataset_computational_complexity/xenium_hbreast_sq40_vector_lst.Rds")))
xhb_clusters= readRDS(file.path(here("data/dataset_computational_complexity/xenium_hbreast_clusters.Rds")))

limma_dt<-decideTests(xhb_fit)
panel_genes=row.names(xhb_seu)
plot_lst=list()
for (cl in c("c1", "c5", "c8")){
    anno_name = unique(xhb_clusters[xhb_clusters$cluster==cl,"anno"])
    anno_name = sub(anno_name,pattern = "_Cells", replacement="")
    limma_sig=row.names(limma_dt[limma_dt[,cl]==1,])
    
    findM_sig =xhb_fm[xhb_fm$cluster==cl & xhb_fm$p_val_adj<0.05,"gene"]
    lasso_sig = xhb_jazzPanda_res[xhb_jazzPanda_res$top_cluster==cl,"gene"]
    
    data_lst = list("limma"=limma_sig,
                    "Wilcoxon Rank Sum Test" =findM_sig,
                    "jazzPanda-glm"=lasso_sig)
    
    df_mt =as.data.frame(matrix(FALSE,nrow=length(panel_genes),ncol=3))
    row.names(df_mt) =panel_genes
    colnames(df_mt)=c("jazzPanda-glm", "Wilcoxon Rank Sum Test",
                      "limma")
    df_mt[findM_sig,"Wilcoxon Rank Sum Test"] = TRUE
    df_mt[limma_sig,"limma"] = TRUE
    df_mt[lasso_sig,"jazzPanda-glm"] = TRUE
    df_mt$gene_name =row.names(df_mt)
    p = plot(upset(df_mt,
                   intersect=c("Wilcoxon Rank Sum Test", "limma",
                               "jazzPanda-glm"),
                   wrap=TRUE, keep_empty_groups= FALSE, name="",
                   #themes=theme_grey(),
                   stripes='white',
                   #stripes=upset_stripes(geom=geom_segment(size=8),colors=c('grey95', 'grey95', 'grey95',"grey95")),
                   sort_intersections_by ="cardinality", sort_sets= FALSE,min_degree=1,
                   set_sizes =( 
                       upset_set_size()
                           # geom_text(aes(label=..count..), size = 4, hjust=1.1, stat='count')
                       + theme(axis.title= element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.text.y = element_blank())),
                   sort_intersections= "descending", warn_when_converting=FALSE,
                   warn_when_dropping_groups=TRUE,encode_sets=TRUE,
                   width_ratio=0.3, height_ratio=1/4)+
                 ggtitle(paste(anno_name,"cells"))+
                 theme(plot.title = element_text( size=20))
             
    )
   plot_lst[[cl]] = p
}
combined_plot <- wrap_plots(plot_lst, ncol = 3)
pdf(file.path(fig5, "xenium_hbreast_upset.pdf"), height=6, width=18)
print(combined_plot)
dev.off()
################################################################################
# cor_M = cor(xhb_vecs$gene_mt,
#             xhb_vecs$cluster_mt,method = "spearman")
# cor_M[cor_M <= 0.7] = 0
# cor_M[cor_M > 0.7]=1
cor1 <- cor(xhb_vecs$gene_mt[1:1600, ],
            xhb_vecs$cluster_mt[1:1600, paste0("c", 1:9)], method = "spearman")

cor2 <- cor(xhb_vecs$gene_mt[1600:3200, ],
            xhb_vecs$cluster_mt[1600:3200, paste0("c", 1:9)], method = "spearman")

cor_M <-(cor1 + cor2) / 2
cluster_names= colnames(limma_dt)
plot_lst=list()
colnames(xhb_seu) <- sub("^_", "", colnames(xhb_seu))
colnames(xhb_seu) <- sub("_2$", "-sp2", colnames(xhb_seu))
rownames(xhb_clusters) <-xhb_clusters$cells
Idents(xhb_seu) <- xhb_clusters$cluster[match(colnames(xhb_seu),
                                          rownames(xhb_clusters))]

for (cl in c("c1", "c5","c8")){
    anno_name = unique(xhb_clusters[xhb_clusters$cluster==cl,"anno"])
    anno_name = sub(anno_name,pattern = "_Cells", replacement="")
    fm_cl=FindMarkers(xhb_seu, ident.1 = cl, only.pos = TRUE,
                      logfc.threshold = 0.1)
    fm_cl = fm_cl[order(fm_cl$avg_log2FC, decreasing = TRUE),]
    to_plot_fm =row.names(fm_cl)
    # FM_pt = data.frame("name"=to_plot,"rk"= 1:length(to_plot),
    #                    "y"= cumsum(cor_M[to_plot, cl]),
    #                    "type"="Wilcoxon Rank Sum Test")
    
    limma_cl<-topTable(xhb_fit,coef=cl,p.value = 0.05, n=Inf, sort.by = "p")
    limma_cl = limma_cl[limma_cl$logFC>0, ]
    to_plot_lm = row.names(limma_cl)
    
    lasso_sig = xhb_jazzPanda_res[xhb_jazzPanda_res$top_cluster==cl,]
    lasso_sig = lasso_sig[order(lasso_sig$glm_coef, decreasing = TRUE),]

    # limma_pt = data.frame("name"=to_plot_lm,"rk"= 1:length(to_plot_lm),
    #                       "y"= cumsum(cor_M[to_plot_lm, cl]),
    #                       "type"="limma")
    # 
    # lasso_pt = data.frame("name"=lasso_sig$gene,"rk"= 1:nrow(lasso_sig),
    #                       "y"= cumsum(cor_M[lasso_sig$gene, cl]),
    #                       "type"="jazzPanda-glm")
    
    
    FM_pt = data.frame("name"=to_plot_fm,"rk"= 1:length(to_plot_fm),
                       "y"= get_cmr_ma(to_plot_fm,cor_M = cor_M, cl = cl),
                       "type"="Wilcoxon Rank Sum Test")
    
    
    lasso_pt = data.frame("name"=lasso_sig$gene,"rk"= 1:nrow(lasso_sig),
                          "y"= get_cmr_ma(lasso_sig$gene,cor_M = cor_M, cl = cl),
                          "type"="jazzPanda-glm")
    
    limma_pt = data.frame("name"=to_plot_lm,"rk"= 1:length(to_plot_lm),
                          "y"= get_cmr_ma(to_plot_lm,cor_M = cor_M, cl = cl),
                          "type"="limma")
    data_lst = rbind(limma_pt, lasso_pt,FM_pt)
    data_lst$type = factor(data_lst$type)
    data_lst$rk = as.numeric(data_lst$rk)
    p <-ggplot(data_lst, aes(x = rk, y = y, color = type)) +
        geom_step(size = 0.8) +  # type = "s"
        scale_color_manual(values = c("jazzPanda-glm" = "red",
                                      "limma" = "black",
                                      "Wilcoxon Rank Sum Test" = "blue")) +

        labs(title = paste(anno_name, "cells"), x = "Rank of marker genes",
             #y = "Cumulative count of genes with correlation >0.7",
             y = "Cumulative average correlation",
             color = NULL) +
        scale_x_continuous(limits = c(0, 50))+
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, size=16),
              axis.line = element_blank(),  
              legend.text  = element_text(size=13),
              panel.border = element_rect(color = "black", 
                                          fill = NA, linewidth = 1),
              legend.position = "inside",
              legend.position.inside = c(0.98, 0.02),
              legend.justification = c("right", "bottom"),
              legend.background = element_rect(color = "black", 
                                               fill = "white", linewidth = 0.5),
              legend.box.background = element_rect(color = "black",
                                                   linewidth = 0.5)
        )
    plot_lst[[cl]] = p
}
combined_plot <- wrap_plots(plot_lst, ncol = 3)
pdf(file.path(fig5, "xenium_hbreast_rank_ma.pdf"), height=5, width=15)
print(combined_plot)
dev.off()