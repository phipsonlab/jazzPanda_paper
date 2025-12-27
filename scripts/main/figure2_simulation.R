################################################################################
# CosMx human liver cancer - simulation

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

################################################################################
# Negative controls are extracted from the full set of transcripts.  
# For details on how different types of negative controls are visualized 
# for each dataset, please refer to the supplementary R Markdown files.
# figure 2(a)
# scripts/supp/supplementary_application_xenium_human_breast_cancer.Rmd
# figure 2(b)
# scripts/supp/supplementary_application_cosmx_healthy_human_liver.Rmd
# figure 2(c)
# scripts/supp/supplementary_application_merscope_human_breast_cancer.Rmd

################################################################################

wk_path = "scripts/main/cosmx_hlc_simulation_result/"

perm_res_files <- list.files(path =wk_path,
                             pattern = "^cosmx_hlc_simulation_s[0-9]+_perm_res_lst\\.Rds$",
                             full.names = TRUE)

p_value_lst <- list()
adjp_value_lst <- list()
obs_lst <- list()


# Loop over the files
for (f in perm_res_files) {
    res_obj <- readRDS(f)
    
    # Optionally extract the seed from the filename
    seed <- sub("^.*_s([0-9]+)_perm_res_lst\\.Rds$", "\\1", basename(f))
    
    # Add a column for the seed to keep track
    res_obj$perm.pval = as.data.frame(res_obj$perm.pval)
    res_obj$obs.stat = as.data.frame(res_obj$obs.stat)
    res_obj$obs.stat$seed <- seed
    res_obj$perm.pval.adj$seed <- seed
    res_obj$perm.pval$seed <- seed
    
    obs_lst[[length(obs_lst) + 1]] <- get_cor(res_obj)
    p_value_lst[[length(p_value_lst) + 1]] <- res_obj$perm.pval
    adjp_value_lst[[length(adjp_value_lst) + 1]] <- get_perm_adjp(res_obj)
}

# Combine into single data frames
obs_df <- do.call(rbind, obs_lst)
adjp_df <- do.call(rbind, adjp_value_lst)
p_value_df <- do.call(rbind, p_value_lst)

colSums(p_value_df[, 1:13] <0.05)
colSums(adjp_df[, 1:13] <0.05)

zz = rowSums(adjp_df[, 1:13] <0.05)

table(zz)
nsig_ps = data.frame(cluster = character(),
                     seed = character(),
                     nsig = numeric(),
                     stringsAsFactors = FALSE)

for (cl in paste("c", 1:13, sep="")){
    curr_df = adjp_df[,c(cl, "seed") ]
    #curr_df = p_value_df[,c(cl, "seed") ]
    dff = curr_df %>% group_by(seed) %>%
        summarise(nsig = sum(.data[[cl]] < 0.05), .groups = "drop")
    nsig_ps = rbind(nsig_ps,cbind(cluster=cl, dff))
    
}

nsig_ps$p = nsig_ps$nsig/100
nsig_df = nsig_ps %>% group_by(cluster) %>% summarise(mean_sig = mean(nsig),
                                            median_sig = median(nsig),
                                            mean_prob = mean(p),
                                            median_p = median(p),
                                            sd_p = sd(p))

data_nm  <- "cosmx_hlc"
# load geenrated data
cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))
colnames(cluster_info)[colnames(cluster_info) == "anno"] <- "anno_name"

cluster_names = c("tumor_1","tumor_2","Macrophages",
                  "T", "Periportal.LSECs","Stellate.cells",
                  "B","Central.venous.LSECs","Cholangiocytes",
                  "Hep","Portal.endothelial.cells","NK.like.cells",
                  "Erthyroid.cells")
cluster_info$cluster_ref =  paste("c",as.numeric(factor(cluster_info$anno_name)), sep="")
cluster_info$cluster_ref =factor(cluster_info$cluster_ref , levels=paste("c",1:13, sep=""))


plt_dff = as.data.frame(table(cluster_info$cluster_ref))
colnames(plt_dff) = c("cluster", "ncells")
dff = merge(plt_dff, nsig_df, by="cluster")
dff$log_ncells = log10(dff$ncells)
cl_anno = as.data.frame(unique(cluster_info[,c("anno_name","cluster_ref")]))
colnames(cl_anno) = c("anno_name","cluster")
dff = merge(dff, cl_anno, by="cluster")
dff$cluster = as.character(dff$cluster)
dff$anno_name = as.character(dff$anno_name)
# dff = rbind(dff, c("falsecode", 0, 0, 0, 0,0,0, 0, "falsecode"), 
#             c("negprobe", 0, 0, 0, 0,0,0, 0,  "negprobe"))
# all_clusters <- c(paste0("c", 1:13), "falsecode", "negprobe")
# 
# dff$cluster = factor(dff$cluster, levels =all_clusters )
dff$cluster = factor(dff$cluster, levels =paste0("c", 1:13) )
dff$log_ncells = as.numeric(dff$log_ncells)
dff$ncells = as.numeric(dff$ncells)
dff$mean_prob = as.numeric(dff$mean_prob)
dff = dff[order(dff$ncells, decreasing = TRUE),]
p1<- ggplot(dff, aes(x = log_ncells, y = mean_prob, color = cluster)) +
    geom_point(size = 5) +
    theme_classic() +
    scale_color_manual(values = my_colors)+
    geom_hline(yintercept = 0.05, colour = "grey80",linewidth = 1,
               linetype = "dashed")+
    scale_y_continuous(breaks = c(0, 0.05, 0.2, 0.4, 0.6, 0.8),
                       labels= c(0, 0.05, 0.2, 0.4, 0.6, 0.8))+
    labs(title = " ", x = " log(Cluster size)",
         y = "Mean false discovery rate", color="") +
         #y = "Mean p-value", color="") +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=1)))+
    theme(
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        axis.text.y = element_text(size=12, vjust=0.5),
        axis.text.x = element_text(size=12, hjust=0.5),
        axis.title.y = element_text(size=12),
        
        legend.text  = element_text(size=12),
        legend.title = element_blank(),          # optional
        legend.key.height = unit(0.5, "cm"),    # shrink vertical space
        legend.key.width  = unit(0.25, "cm"),    # optional: shrink horizontally
        legend.margin     = margin(0, 0, 0, 0)   # remove padding
    )

pdf(file.path(fig2, "simulation_ncells_nsig_adj_pval_scatter.pdf"), 
    width=6, height=4)
p1
dev.off()
################################################################################

jazzPanda_perm = as.data.frame(matrix(0, nrow=ncol(obs_df), ncol=2))
colnames(jazzPanda_perm) = c("cluster","jazzPanda_correlation")
jazzPanda_perm$cluster = colnames(obs_df)

for (cl in paste("c", 1:13, sep="")){
    obs_cutoff=0.05
    #obs_cutoff = quantile(obs_df[, cl], 0.95)
    perm_cl=intersect(row.names(adjp_df[adjp_df[,cl]<0.05,]),
                      row.names(obs_df[obs_df[, cl]>obs_cutoff,]))
    jazzPanda_perm[jazzPanda_perm$cluster==cl,"jazzPanda_correlation"] = length(perm_cl)
}

cat("False positive rate =", (sum(jazzPanda_perm[,2])/10000*100),"%")


all_files = list.files(path = wk_path, pattern = "\\.Rds$", full.names = TRUE)
res_list <- list()  # empty list to store results
seeds <- str_extract(all_files, "s\\d+")
seeds <- unique(seeds)

kep_cols= c("seed","gene","falsecode","negprobe",paste("c", 1:13, sep=""))
for (seed in seeds) {
    
    # Find all files related to this seed
    seed_files <- all_files[grepl(seed, all_files)]
    
    # Load each file
    nc_spatial_file <- seed_files[grepl("_nc_spatial_vectors\\.Rds$", seed_files)]  
    spatial_file <- seed_files[grepl("_spatial_vectors\\.Rds$", seed_files) & !grepl("_nc_spatial_vectors\\.Rds$", seed_files)]  
    nc_spatial = readRDS(nc_spatial_file)
    sv_spatial = readRDS(spatial_file)
    curr_corr= as.data.frame(cor(sv_spatial$gene_mt,
                                 cbind(nc_spatial,sv_spatial$cluster_mt[,paste("c", 1:13, sep="")]), 
                                 method = "pearson"))
    curr_corr$seed = seed
    curr_corr$gene = row.names(curr_corr)
    curr_corr=curr_corr[, kep_cols]
    res_list[[seed]] <- curr_corr
    #corr_df = rbind(curr_corr, corr_df)
}
corr_df <- do.call(rbind, res_list)


dff = reshape2::melt(corr_df, id.vars =c("gene", "seed"), variable.name = "cluster",
                     value.name = "observed_corr")
dff$cluster = factor(dff$cluster, levels = c(paste("c", 1:13, sep=""), "falsecode","negprobe"))

summary(dff[!(dff$cluster %in% c("falsecode","negprobe")),"observed_corr"])
summary(dff[dff$cluster == "falsecode","observed_corr"])
summary(dff[dff$cluster == "negprobe","observed_corr"])

avg_dff = dff %>% group_by(seed, cluster)%>% summarise(mean_corr = mean(observed_corr),
                                          median_corr= median(observed_corr),
                                          sd_corr= sd(observed_corr),
                                          .groups = "drop")
p1<- ggplot(avg_dff, aes(x = cluster, y = mean_corr)) +
    #geom_point(size = 1) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6)+
    #geom_violin(scale = "width")+
    theme_classic() +
    labs(title = " ", x = " ", y = "Mean pearson correlation") +
    theme( axis.line = element_blank(),
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
           axis.text.y = element_text(size=10, vjust=0.5),
           axis.text.x = element_text(size = 10, angle =45, hjust = 1),
           axis.title.y= element_text(size=12))

pdf(file.path(fig2, "figure2_cosmx_hlc_correlation.pdf"), width=5, height=4)
p1
dev.off()


################################################################################
lasso_files <- list.files(path =wk_path,
                          pattern = "^cosmx_hlc_simulation_s[0-9]+_lasso_res_lst\\.Rds$",
                          full.names = TRUE)

top_df_lst <- list()
full_df_lst<- list()

# Loop over the files
for (f in lasso_files) {
    res_obj <- readRDS(f)
    
    # Optionally extract the seed from the filename
    seed <- sub("^.*_s([0-9]+)_lasso_res_lst\\.Rds$", "\\1", basename(f))
    
    # Add a column for the seed to keep track
    res_obj$top_result$seed <- seed
    res_obj$full_result$seed <- seed
    
    top_df_lst[[length(top_df_lst) + 1]] <- get_top_mg(res_obj, coef_cutoff = 0)
    full_df_lst[[length(full_df_lst) + 1]] <- get_full_mg(res_obj, coef_cutoff = 0)
}

# Combine into single data frames
top_df <- do.call(rbind, top_df_lst)
full_df <- do.call(rbind, full_df_lst)

table(top_df$top_cluster)

# check the top cluster from full table
top_by_gene <- full_df %>%
    arrange(seed, gene, p_value, desc(glm_coef)) %>%  # break ties in p by larger coef
    group_by(seed, gene) %>%
    slice_head(n = 1) %>%
    ungroup()
table(top_by_gene$cluster)


## no genes with real cluster as only fetures 
genes_only_cluster <- full_df %>%
    group_by(seed,gene) %>%
    summarise(all_c = all(grepl("^c", cluster))) %>%
    filter(all_c) %>%
    nrow()

n_summary = full_df %>%
    mutate(
        type = case_when(
            grepl("^c", cluster) ~ "real",
            cluster %in% c("falsecode", "negprobe") ~ "background",
            TRUE ~ "other"
        )
    ) %>%
    group_by(seed, gene) %>%
    summarise(
        n_real_clusters = sum(type == "real"),
        n_background_clusters = sum(type == "background")
    )

nrow(n_summary[n_summary$n_real_clusters==0,])
nrow(n_summary[n_summary$n_real_clusters==1,])
nrow(n_summary[n_summary$n_real_clusters==2,])
nrow(n_summary[n_summary$n_real_clusters==3,])
nrow(n_summary[n_summary$n_real_clusters==4,])

ranked_df <- full_df %>%
    group_by(gene, seed) %>%
    #arrange(desc(glm_coef), .by_group = TRUE) %>%
    arrange(p_value, .by_group = TRUE) %>%
    mutate(p_rank = row_number()) %>%
    ungroup()

################################################################################
# # Ensure rank is treated as a factor (for full axis)
# top_df <- ranked_df %>%
#     mutate(rank = as.integer(p_rank),
#            cluster = as.factor(cluster))
# 
# # Get all combinations of cluster and rank (e.g., ranks 1:6)
# full_combos <- expand_grid(
#     cluster = unique(top_df$cluster),
#     rank = 1:6
# )
# 
# # Count frequencies and fill in missing ones with 0
# plot_df <- top_df %>%
#     dplyr::count(cluster, rank) %>%
#     right_join(full_combos, by = c("cluster", "rank")) %>%
#     replace_na(list(n = 0))
# 
# plot_df$cluster <- factor(plot_df$cluster, levels = c(
#     "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", 
#     "c10", "c11", "c12", "c13", "falsecode", "negprobe"
# ))
# 
# pdf(file.path(fig2, "rank_p_value_heatmap.pdf"), width=8,height=5)
# ggplot(plot_df, aes(x = cluster, y = factor(rank), fill = n)) +
#     geom_tile(color = "white", linewidth = 0.2) +
#     scale_fill_viridis_c(
#         option = "plasma",
#         #trans = "sqrt",
#         name = "# Simulation",
#         breaks = seq(0, 10000, 2000),
#         guide = guide_colorbar(barwidth = 2, barheight = 10)
#     ) +
#     labs(
#         x = "cluster",
#         y = "cluster rank (1 = most significant)",
#         title = "Rank of clusters across simulated genes",
#         subtitle = "Lower ranks correspond to lower p-value"
#     ) +
#     theme_minimal(base_size = 14) +
#     theme(
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid = element_blank(),
#         plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
#         plot.subtitle = element_text(size = 12, hjust = 0.5),
#         axis.title.x = element_text(margin = margin(t = 10)),
#         axis.title.y = element_text(margin = margin(r = 10))
#     )
# dev.off()


################################################################################
library(tidyverse)
top_df <- ranked_df %>%
    mutate(rank = as.integer(p_rank))




# expand grid → every seed × gene × cluster combo
df_full <- top_df[, c("seed", "gene", "rank", "cluster")] %>%
    mutate(cluster = as.character(cluster),
           rank    = as.integer(rank)) %>%  
    group_by(seed, gene) %>%
    tidyr::complete(cluster = all_clusters, fill = list(rank = 7L)) %>%
    ungroup() %>%
    mutate(cluster = factor(cluster, levels = all_clusters))

# summarise across genes for each seed × cluster
df_reps <-  df_full %>%
    group_by(seed, cluster) %>%
    summarise(
        mean_rank   = mean(rank),
        median_rank = median(rank),
        .groups = "drop"
    )
df_reps$cluster = factor(df_reps$cluster, levels = all_clusters)
pdf(file.path(fig2, "rank_p_value_heatmap_per_run.pdf"), width=10,height=5)
ggplot(df_reps, aes(x = seed, y = cluster, fill = mean_rank)) +
    geom_tile() +
    scale_fill_viridis_c(
        option    = "plasma",
        direction = -1,  
        name      = "Mean rank",
        breaks    = 1:7,   # labels from 1 to 7
        limits    = c(1, 7), # clamp the scale
        guide     = guide_colorbar(barwidth = 2, barheight = 10)
    )+
    labs(
        x = "100 Simulation runs",
        y = " ",
        title = "Rank of clusters across simulated genes",
        subtitle = "Lower ranks correspond to lower p-value"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title= element_text(size=12), 
        panel.grid = element_blank()
    )

dev.off()

################################################################################
# 
# dff = full_df[, c("gene", "seed", "glm_coef", "cluster", "p_value")]%>%
#     mutate(cluster = as.character(cluster)) %>%  
#     group_by(seed, gene) %>%
#     tidyr::complete(cluster = all_clusters, fill = list(p_value = 1,glm_coef=0)) %>%
#     ungroup() %>%
#     mutate(cluster = factor(cluster, levels = all_clusters))
# dff$cluster = factor(dff$cluster, levels = c(paste("c", 1:13, sep=""), 
#                                              "falsecode","negprobe"))
# dff$glm_coef = as.numeric(dff$glm_coef) 
# summary(dff[!(dff$cluster %in% c("falsecode","negprobe")),"glm_coef"])
# summary(dff[dff$cluster == "falsecode","glm_coef"])
# summary(dff[dff$cluster == "negprobe","glm_coef"])
# 
# avg_dff = dff %>% group_by(seed, cluster)%>% summarise(mean_glm = mean(glm_coef),
#                                                        median_glm= median(glm_coef),
#                                                        sd_glm = sd(glm_coef),
#                                                        mean_p = mean(p_value),
#                                                        median_p= median(p_value),
#                                                        sd_p = sd(p_value),
#                                                        .groups = "drop")
# p1<- ggplot(avg_dff, aes(x = cluster, y = mean_p)) +
#     #geom_point(size = 1) +
#     geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6)+
#     #geom_violin(scale = "width")+
#     theme_classic() +
#     geom_hline(yintercept = 0.05, colour = "grey80",linewidth = 1,
#                linetype = "dashed")+
#     scale_y_continuous(breaks = c(0, 0.05, 0.2, 0.4, 0.6, 0.8,1),
#                        labels= c(0, 0.05, 0.2, 0.4, 0.6, 0.8,1))+
#     labs(title = " ", x = " ", y = "Mean p-value") +
#     theme( axis.line = element_blank(),
#            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
#            axis.text.y = element_text(size=10, vjust=0.5),
#            axis.text.x = element_text(size = 10, angle =45, hjust = 1),
#            axis.title.y= element_text(size=12))
# 
# pdf(file.path(fig2, "figure2_cosmx_hlc_boxplot_glm_pval.pdf"), width=8, height=4)
# p1
# dev.off()
# 
# p2<- ggplot(avg_dff, aes(x = cluster, y = median_glm)) +
#     #geom_point(size = 1) +
#     geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6)+
#     #geom_violin(scale = "width")+
#     theme_classic() +
# 
#     labs(title = " ", x = " ", y = "Median model coefficient") +
#     theme( axis.line = element_blank(),
#            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
#            axis.text.y = element_text(size=10, vjust=0.5),
#            axis.text.x = element_text(size = 10, angle =45, hjust = 1),
#            axis.title.y= element_text(size=12))
# 
# pdf(file.path(fig2, "figure2_cosmx_hlc_boxplot_glm_coef.pdf"), width=8, height=4)
# p2
# dev.off()

