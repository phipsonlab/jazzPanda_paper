# define figure and data path
library(here)
data_path <- here::here("data", "dataset_computational_complexity")
overview_PA <- here::here("figures", "supp", "application")
mg_PA <- here::here("figures", "supp", "application", "marker_genes")
comp_PA <- here::here("figures", "supp")
mg_ntiles_comp <- here::here("scripts", "main","discussion_markergenes_vs_ntiles")


fig2 <- here::here("figures", "main","figure2_cosmx_hliver")
fig3 <- here::here("figures", "main","figure3_xenium_hbreast")
fig4 <- here::here("figures", "main","figure4_simulation")
fig5 <- here::here("figures", "main","figure5_compare_methods")
fig6 <- here::here("figures", "main","figure6")



# cluster colors
my_colors <- c(
    "#800000", "#f032e6","#ffe119", "#e6beff",
    "#e6194b", "#008080", "#0082c8", "#3cb44b",
    "#aa6e28", "#fabebe", "#d2f53c", "#911eb4",
    "#fffac8", "#f58231", "#46f0f0","#808080"
)

#############################################################################
# theme setting for figures
defined_theme <- theme(strip.text = element_text(size = rel(2)),
                       axis.line=element_blank(),
                       axis.ticks=element_blank(),
                       axis.text=element_blank(),
                       legend.position = "none",
                       axis.title=element_blank(),
                       panel.background=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.background=element_blank(),
                       panel.border = element_rect(fill=NA, colour = "black"))


#############################################################################
# Function to automatically determine hex bin size based on data density
auto_hex_bin <- function(n, target_points_per_bin = 5) {
    k <- n / target_points_per_bin
    bins <- round(sqrt(k))
    return(max(10, bins)) 
}
#############################################################################
# A help function to load Xenium data
get_xenium_data<-function(path,mtx_name, trans_name="transcript_info.csv.gz",
                          cells_name="cell_info.csv.gz"){
    
    transcript_info <- as.data.frame(fread(paste(path, trans_name,sep="")))
    cell_info <- as.data.frame(fread(paste(path,cells_name,sep="")))
    
    data <- Read10X(data.dir = paste(path,mtx_name, sep=""))
    
    cm <- as.matrix(data$`Gene Expression`)
    r_codeword <- as.matrix(data$`Negative Control Codeword`)
    
    r_probe <- as.matrix(data$`Negative Control Probe`)
    # merge negative control genes and real genes
    cm_neg <- as.data.frame(rbind(r_probe, r_codeword))
    zero_cells <- colnames(cm)[colSums(cm)==0]
    
    transcript_info$x <- as.numeric(transcript_info$x_location)
    transcript_info$y <- as.numeric(transcript_info$y_location)
    
    cell_info$x <- as.numeric(cell_info$x_centroid)
    cell_info$y <- as.numeric(cell_info$y_centroid)
    
    return (list(cm = cm, cm_neg=cm_neg, zero_cells = zero_cells,
                 trans_info=transcript_info, cell_info=cell_info,
                 probe = row.names(r_probe),
                 codeword=row.names(r_codeword)))
    
}

#############################################################################
# help functions for figures 
plot_cluster_props <- function(cluster_info, file_prefix,
                               ct_nm,
                               my_colors,
                               out_dir = overview_PA) {
    sum_df <- cluster_info %>%
        dplyr::group_by(!!sym(ct_nm), sample) %>%
        dplyr::count(name = "n")
    
    # Rename the grouping column to "cellType" for plotting consistency
    colnames(sum_df)[colnames(sum_df) == ct_nm] <- "cellType"
    
    
    n_sample = length(unique(cluster_info$sample))
    # Build output file path
    out_file <- here(out_dir, paste0(file_prefix, "_cluster_prop", ".pdf"))
    
    # Save to PDF
    pdf(out_file, width = 6*min(1.3,n_sample), height = 10)
    print(ggplot(sum_df, aes(x = sample, y = n, fill = cellType)) +
              geom_bar(stat = "identity", position = "fill") +
              scale_y_continuous(labels = scales::percent, 
                                 expand = c(0.01, 0.01)) +
              labs(title = " ", x = " ", y = "proportion of cell types", 
                   fill = "cell types") +
              scale_fill_manual(values = my_colors) +
              theme(
                  legend.position = "right",
                  axis.text.x = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  axis.line = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_blank()
              ))
    dev.off()
}

#############################################################################
# help function to plot spatial coordinates for cells in a dataset
plot_data_sp <- function(cluster_info, file_prefix,
                         ct_nm,
                         my_colors,
                         out_dir = overview_PA,reverse_y = FALSE) {
    
    colnames(cluster_info)[colnames(cluster_info) == ct_nm] <- "cellType"
    n_sample = length(unique(cluster_info$sample))
    # Build output file path
    out_file <- here(out_dir, paste0(file_prefix, "_dataset_sp", ".jpg"))
    
    jpeg(out_file, width = 1000*n_sample, height = 900, res=200)
    print(ggplot(data = cluster_info,
                 aes(x = x, y = y, color=cellType))+
              geom_point(size=0.001, alpha=0.7)+
              facet_wrap(~sample, nrow=1)+
              theme_classic()+
              (if (reverse_y) scale_y_reverse() else NULL) +
              scale_color_manual(values = my_colors)+
              guides(color=guide_legend(title=" ", ncol = 1,
                                        override.aes=list(alpha=1, size=8)))+ 
              defined_theme+
              theme(strip.background = element_rect(color="white"),
                    legend.position = "none",legend.key.width = unit(1, "cm"),
                    legend.title = element_text(size = 12),
                    legend.text = element_text(size = 8),
                    legend.key.height = unit(1, "cm"),
                    strip.text = element_text(size = 12))
    )
    dev.off()
}


#############################################################################
# help function to plot UMAP
plot_umap_seu <- function(cluster_info, seu, file_prefix,
                          ct_nm,
                          my_colors,
                          out_dir = overview_PA, fig_w=1000) {
    
    
    # Build output file path
    out_file <- here(out_dir, paste0(file_prefix, "_UMAP", ".jpg"))
    # Save to PDF
    jpeg(out_file, width = fig_w, height = 1000, res=200)
    print(DimPlot(seu, reduction = "umap",split.by = "sample")+  ggtitle("") +
              labs(title = " ", x = "UMAP1", y = "UMAP2", fill=" ") +
              scale_color_manual(values = my_colors)+
              theme(legend.box.margin = margin(0, 0, 0, 2), 
                    plot.margin = margin(5, 5, 5, 5),     
                    legend.position = "right",
                    panel.border = element_blank(), 
                    panel.grid = element_blank(),
                    legend.text = element_text(size = 12),
                    panel.background = element_blank(),
                    axis.title = element_text(size=12)
              ) )
    dev.off()
}