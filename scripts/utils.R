# define figure and data path
library(here)
data_path <- here::here("data", "dataset_computational_complexity")
overview_PA <- here::here("figures", "supp", "application")
mg_PA <- here::here("figures", "supp", "application", "marker_genes")
comp_PA <- here::here("figures", "supp")
mg_ntiles_comp <- here::here("scripts", "main","discussion_markergenes_vs_ntiles")

MERSCOPE_RAW_DATA <- "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Merscope_human_breast_sample/"
fig2 <- here::here("figures", "main", "figure2_simulation")
fig3 <- here::here("figures", "main", "figure3_cosmx_hliver")
fig4 <- here::here("figures", "main", "figure4_xenium_hbreast")
fig5 <- here::here("figures", "main", "figure5_compare_methods")
fig6 <- here::here("figures", "main", "figure6_sv_extension")
fig7 <- here::here("figures", "main", "figure7_technical_performance")

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
auto_hex_bin <- function(n, target_points_per_bin = 0.5, min_bins = 15, max_bins = 200) {
    if (n < 50) {
        return(min_bins)
    }
    k <- n / target_points_per_bin
    bins <- round(sqrt(k))
    
    # For highly expressed genes, allow more bins to resolve spatial structure
    # Scale up target more aggressively past a threshold
    if (n > 5000) {
        bins <- round(sqrt(n / 3))
    }
    
    return(min(max_bins, max(min_bins, bins)))
}

# Function to calculate the cumulative average correlation 
get_cmr_ma<- function(genes, cor_M, cl){
    curr_corrs = cor_M[genes, cl]
    mv_avgs = cumsum(curr_corrs)/(1:length(genes))
    return(mv_avgs)
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
    
    cm = cm[, setdiff(colnames(cm),zero_cells)]
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
                         my_colors,height = 900,fig_ratio=1,
                         out_dir = overview_PA,reverse_y = FALSE) {
    
    colnames(cluster_info)[colnames(cluster_info) == ct_nm] <- "cellType"
    n_sample = length(unique(cluster_info$sample))
    # Build output file path
    out_file <- here(out_dir, paste0(file_prefix, "_dataset_sp", ".jpg"))
    
    jpeg(out_file, width = 1000*n_sample, height = height, res=200)
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
                    aspect.ratio = fig_ratio,
                    legend.text = element_text(size = 8),
                    legend.key.height = unit(1, "cm"),
                    strip.text = element_blank())
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

#############################################################################

plot_nc <- function(df,
                    sample_name,
                    category,
                    data_nm,
                    overview_PA,
                    bins_map = list("negative control probe" = 2,
                                    "negative control codeword" = 1),
                    reverse_y=FALSE, fig_ratio = 1,
                    width,
                    height) {
    # skip if no data
    if (!category %in% df$category) {
        message("No data for category: ", category, " in ", sample_name)
        return(invisible(NULL))
    }
    
    dff <- df[df$category == category, ]

    # choose target_points_per_bin from bins_map (fallback = 2)
    target_bin <- ifelse(category %in% names(bins_map),
                         bins_map[[category]],
                         2)
    
    out_file <- file.path(
        overview_PA,
        paste0(data_nm, "_", sample_name, "_",
               gsub(" ", "_", category), "_nc_xy.jpg")
    )
    base_size <- 18
    scale_factor <- sqrt(width * height) / 1000
    
    jpeg(out_file, width = width, height = height, res = 200)
    print(
        ggplot(dff, aes(x = x, y = y)) +
            geom_hex(bins = auto_hex_bin(nrow(dff), target_points_per_bin = target_bin)) +
            theme_bw() +
            scale_fill_gradient(low = "white", high = "black") +
            guides(fill = guide_colorbar(barheight = unit(0.06, "npc"),
                                         barwidth  = unit(0.4, "npc"),
                                         )) +
            (if (reverse_y) scale_y_reverse() else NULL) +
            facet_grid(~sample) +
            defined_theme +
            theme(legend.title = element_text(size = base_size * scale_factor,
                                              hjust = 1, vjust=0.8),
                legend.position = "bottom",
                legend.text = element_text(size = 0.8*base_size * scale_factor),
                aspect.ratio = fig_ratio,
                strip.text = element_blank(),
                strip.background = element_blank()
            )
    )
    dev.off()
}
