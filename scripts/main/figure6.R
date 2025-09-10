################################################################################
# Figure6: Technical performance
library(dplyr) 
library(knitr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ComplexUpset)
library(xtable)
library(spatstat.geom)
library(ggtext)
library(here)
source(here("scripts/utils.R"))
################################################################################
# Computational complexity 
## Time and memory usage of analysed datasets
cosmx_cancer = read.csv(here("data/dataset_computational_complexity/cosmx_human_liver_cancer_5core.csv"))
cosmx_normal= read.csv(here("data/dataset_computational_complexity/cosmx_human_healthy_liver_5core.csv"))
merscope_hbreast = read.csv(here("data/dataset_computational_complexity/merscope_human_breast_cancer_5core.csv"))
xenium_mbrain =read.csv(here("data/dataset_computational_complexity/xenium_mouse_brain_5core.csv"))
xenium_hbreast=read.csv(here("data/dataset_computational_complexity/xenium_human_breast_cancer_5core.csv"))
xenium_hlung = read.csv(here("data/dataset_computational_complexity/xenium_human_lung_cancer_5core.csv"))
update_resource_usage <- function(df) {
    df_long = reshape2::melt(df, 
                             id.var=c("technology","dataset", "genes_n",
                                      "transcript_n", "cells_n", "clusters_n",
                                      "grid_length"), variable.name="type", 
                             value.name="value")
    df_long$type = as.character(df_long$type)
    df_long$name=  sapply(df_long$type, function(x) strsplit(x, "_")[[1]][1])
    df_long$usage=   sapply(df_long$type, function(x) strsplit(x, "_")[[1]][2])
    df_long = df_long[,c("technology","dataset", "genes_n", "transcript_n",
                         "cells_n", "clusters_n", "grid_length", 
                         "value","name", "usage" )]
    df_wide <- df_long %>% pivot_wider(names_from = usage, 
                                       values_from = value )
    df_wide = as.data.frame(df_wide)
    df_wide$Elapsed = df_wide$Elapsed/60
    
    colnames(df_wide) = c("technology","dataset", "genes_n", "transcript_n",
                          "cells_n", "clusters_n", "grid_length", "name", 
                          "time_min","peak_memory_mb" )
    
    df_wide[df_wide$name=="glm","time_min"]=df_wide[df_wide$name=="glm","time_min"]+df_wide[df_wide$name=="sv","time_min"]
    df_wide[df_wide$name=="glm","peak_memory_mb"]=df_wide[df_wide$name=="glm","peak_memory_mb"]+df_wide[df_wide$name=="sv","peak_memory_mb"]
    df_wide = df_wide[df_wide$name != "sv", ]
    
    return(df_wide)
}

# Apply the function to each data frame
xenium_mbrain <- update_resource_usage(xenium_mbrain)
xenium_hlung <- update_resource_usage(xenium_hlung)
cosmx_cancer <- update_resource_usage(cosmx_cancer)
cosmx_normal <- update_resource_usage(cosmx_normal)

xenium_hbreast <- update_resource_usage(xenium_hbreast)
merscope_hbreast <- update_resource_usage(merscope_hbreast)

dff = rbind(xenium_hlung,cosmx_cancer,merscope_hbreast,
            cosmx_normal,xenium_mbrain, xenium_hbreast)
labels_refs = as.data.frame(cbind(name = c("perm","glm", "limma","fm"), 
                                  label = c("jazzPanda-correlation(5000 permutation)",
                                            "jazzPanda-glm",
                                            "limma","Wilcoxon Rank Sum tests")))
dff = merge(dff, labels_refs, by="name" )
dff$label = factor(dff$label,levels =labels_refs$label)

dff_data_summary <- dff %>%  distinct(dataset,.keep_all = TRUE) %>% 
    select(technology, dataset, genes_n,transcript_n, cells_n,  clusters_n)
dff_data_summary = dff_data_summary[order(dff_data_summary$transcript_n, decreasing = FALSE),]
xtable(dff_data_summary, format = "latex", 
       display = c("s","s","d","d","e","e","d"),
       caption = "tested datasets",  include.rownames=FALSE)

################################################################################
# Memory/ Time complexity

dff$dataset_info = paste0(dff$dataset,"<br>", "*", dff$cells_n, 
                          " cells <br>",dff$transcript_n," transcripts",
                          "*",  sep="")
dff$dataset_info = factor(dff$dataset_info, 
                          levels =paste0(dff_data_summary$dataset,"<br>", "*", 
                                         dff_data_summary$cells_n, 
                                         " cells <br>",dff_data_summary$transcript_n,
                                         " transcripts","*",  sep=""))
dff = dff[order(dff$dataset_info), ]
pdf(file.path(fig6, "time_usage_per_dataset.pdf"), width=16, height=5)
ggplot(dff, aes(x = dataset_info, y = time_min, fill =label )) +
    geom_bar(stat = "identity", position = position_dodge())+
    theme_classic() +
    labs(title = "", x = "", y ="Time (min)", fill = " ") +
    guides(fill = guide_legend(nrow = 1, size=2))+
    scale_y_continuous(expand = c(0.02,0.02)) +
    theme(axis.text.y=element_text( vjust = 0.5, hjust=0.5, size=18),
          axis.text.x=element_markdown(angle=0, vjust = 0.5, hjust=0.5, size=12),
          legend.position = "top",
          axis.title.y =  element_text(size=18),
          legend.title = element_blank(), 
          legend.key.spacing.x= unit(1.5, 'cm'),     
          legend.text = element_text(size=15))

dev.off()


pdf(file.path(fig6, "memory_usage_per_dataset.pdf"), width=16, height=5)
ggplot(dff, aes(x = dataset_info, y = peak_memory_mb/1024, fill = label)) +
    geom_bar(stat = "identity", position = position_dodge())+
    theme_classic() +
    labs(title = "", x = "", y ="Memory (GB)", fill = " ") +
    scale_y_continuous(expand = c(0.02,0.02)) +
    theme(axis.text.y=element_text( vjust = 0.5, hjust=0.5, size=18),
          axis.text.x=element_markdown(angle=0, vjust = 0.5, hjust=0.5, size=12),
          legend.position = "top",
          axis.title.y =  element_text(size=18),
          legend.title = element_blank(), 
          legend.key.spacing.x= unit(1.5, 'cm'),     
          legend.text = element_text(size=15))
dev.off()


################################################################################
# Marker gene comparision vs number of tiles 

# get all slurm array results for marker genes vs ntiles (xenium human breast)
xenium_files <- list.files(path = mg_ntiles_comp, pattern = "^xenium.*\\.csv$", full.names = TRUE)
# Read all files into a list of data frames
data_list <- lapply(xenium_files, read.csv)

keep_dfs <- lapply(data_list, function(df) df[, c(1,4,5)])
xenium_merged_df <- Reduce(function(x, y) merge(x,y,
                                                by = names(x)[1], all = TRUE), keep_dfs)
colnames(xenium_merged_df)[1] = "gene"

## Xenium huamn breast cancer samples
clusters_df <- data.frame(
    cluster = c("c1", "c3", "c2", "c4", "c5", "c7", "c9", "c6", "c8"),
    anno = c("Tumor", "Macrophages", "Stromal", "Myoepithelial", "T", 
             "Endothelial", "Mast", "B", "Dendritic"),
    stringsAsFactors = FALSE
)
plt_lst = list()
for (cl in c("c1", "c5","c8")){
    anno_name = clusters_df[clusters_df$cluster==cl,"anno"]
    gr10 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr10"]==cl,"gene"]
    gr20 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr20"]==cl,"gene"]
    gr30 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr30"]==cl,"gene"]
    gr40 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr40"]==cl,"gene"]
    gr50 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr50"]==cl,"gene"]
    gr60 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr60"]==cl,"gene"]
    gr70 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr70"]==cl,"gene"]
    gr80 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr80"]==cl,"gene"]
    gr90 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr90"]==cl,"gene"]
    gr100 = xenium_merged_df[xenium_merged_df[,"top_cluster_gr100"]==cl,"gene"]
    
    df_mt =as.data.frame(matrix(FALSE,nrow=nrow(xenium_merged_df),ncol=10))
    row.names(df_mt) =xenium_merged_df$gene
    colnames(df_mt)=paste(paste(seq(10,100, by=10), "x", sep=""),
                          seq(10,100, by=10),sep="")
    df_mt[gr10,"10x10"] = TRUE
    df_mt[gr20,"20x20"] = TRUE
    df_mt[gr30,"30x30"] = TRUE
    df_mt[gr40,"40x40"] = TRUE
    df_mt[gr50,"50x50"] = TRUE
    df_mt[gr60,"60x60"] = TRUE
    df_mt[gr70,"70x70"] = TRUE
    df_mt[gr80,"80x80"] = TRUE
    df_mt[gr90,"90x90"] = TRUE
    df_mt[gr100,"100x100"] = TRUE
    p<-upset(df_mt,intersect=colnames(df_mt),
                               wrap=TRUE, keep_empty_groups= FALSE, name="",
                               themes=theme_grey(),
                               stripes=upset_stripes(geom=geom_segment(size=5),
                                                     colors=c('grey95', 'grey95', 'grey95')),
                               sort_intersections_by ="cardinality", sort_sets= FALSE,min_degree=1,
                               set_sizes = FALSE,
                               sort_intersections= "descending", warn_when_converting=FALSE,
                               warn_when_dropping_groups=TRUE,encode_sets=TRUE,
                               matrix=(intersection_matrix()+
                                           theme(axis.text.x=element_blank(),
                                                 panel.background = element_rect(fill="NA"),
                                                 axis.ticks = element_blank(),
                                                 axis.title = element_blank())),
                               base_annotations=list('Intersection size'=(intersection_size(bar_number_threshold=1,color='grey9',fill='grey80',
                                                                                            text = list(size = 3, vjust = -0.1))+
                                                                              theme(axis.text.x = element_blank(),
                                                                                    axis.title.x = element_blank(),
                                                                                    panel.background = element_rect(fill="NA"),
                                                                                    panel.grid = element_line(color="grey90"),
                                                                                    axis.ticks.x = element_blank()))),
                               width_ratio=0.5, height_ratio=0.7)+
                             ggtitle(paste(anno_name,"cells"))+
                             theme(
                                 plot.title = element_text(size = 18, face = "bold")
                             )
    plt_lst[[cl]] <-p
}  
combined_plot <- wrap_plots(plt_lst, ncol = 3)
pdf(file.path(fig6, "xenium_mg_ntiles.pdf"), height=6, width=18)
print(combined_plot)
dev.off()



