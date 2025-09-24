
library(SpatialExperimentIO)
library(SpatialExperiment)
library(Seurat)
library(Banksy)
se =  readMerscopeSXE(dirName = "/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/HumanBreastCancerPatient1/",
                      countMatPattern = "cell_by_gene.csv", metaDataPattern = "cell_metadata.csv")
x_avg <- (se@colData$min_x + se@colData$max_x) / 2
y_avg <- (se@colData$min_y + se@colData$max_y) / 2

# create a matrix
coords <- cbind(x = x_avg, y = y_avg)

# assign to spatialCoords
SpatialExperiment::spatialCoords(se) <- coords

blank_genes <- grep("^Blank", rownames(se), value = TRUE)

blank_genes

blank_idx <- grep("^Blank", rownames(se))
se <- se[-blank_idx, ]

lib_size <- Matrix::colSums(assay(se, "counts"))
colData(se)$lib_size <- lib_size

# 3. Filter cells with library size between 0 and 2500
se <- se[, lib_size > 50 & lib_size < 2500]

seu <- as.Seurat(se, data = NULL)
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
seu <- FindVariableFeatures(seu, nfeatures = nrow(seu))

# copy log-normalized data back to se
aname <- "logcounts"
logcounts_mat <- GetAssayData(seu, slot = "data")[, colnames(se)]
assay(se, aname, withDimnames = FALSE) <- logcounts_mat

lambda <- 0.2
k_geom <- 15
use_agf <- TRUE

se <- computeBanksy(se, assay_name = aname, 
                     compute_agf = TRUE, k_geom = k_geom)

se <- runBanksyPCA(se, use_agf = use_agf, lambda = lambda, seed = 1000)
se <- runBanksyUMAP(se, use_agf = use_agf, lambda = lambda, seed = 1000)

cat("Clustering starts\n")
se <- clusterBanksy(se, use_agf = use_agf, lambda = lambda,
                     resolution = c(0.8, 1, 1.2), seed = 1000)

se <- connectClusters(se)


se = readRDS("/vast/projects/xenium_5k/jazzPanda_paper/scripts/main/merscope_hbreast_se.Rds")

cluster_info <- data.frame(
    x = spatialCoords(se)[, 1],
    y = spatialCoords(se)[, 2],
    cell_id = colnames(se), 
    cluster = paste0("c", colData(se)$clust_M1_lam0.2_k50_res0.5),  
    sample = colData(se)$sample_id
)

umap_bkm <- as.data.frame(reducedDim(se, "UMAP_M1_lam0.2"))
umap_bkm$cluster <- paste0("c", colData(se)$clust_M1_lam0.2_k50_res0.5)
umap_bkm$sample_id <- se$sample_id  
# umap_bkm$anno_name <- cluster_info$anno_name[match(row.names(umap_bkm),
#                                                    cluster_info$cell_id)]
# umap_bkm$anno_name =factor(umap_bkm$anno_name, 
#                            levels = anno_df$anno_name)
# umap_bkm = umap_bkm[order(umap_bkm$anno_name), ]

# rm(se)


jpeg("umap.jpg", width=1000, height=1000, res=200)

ggplot(umap_bkm, aes(x = V1, y = V2, color = cluster)) +
    ggrastr::geom_point_rast(size = 0.01, alpha = 0.8) +
    # defined_theme +
    facet_wrap(~sample_id, nrow=1)+
    # scale_color_manual(values = my_colors)+
    # guides(color = guide_legend(override.aes = list(size = 3), 
    #                             ncol = 1, title = "")) +
    labs(title = " ", x = "UMAP1", y = "UMAP2", fill=" ") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 20, face = "plain"),
          legend.position = "right",
          legend.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", 
                                      fill = NA, linewidth = 0.3),
          axis.line = element_blank(),
          legend.justification = "center",
          axis.ticks = element_blank(),
          legend.box.just = "center", 
          axis.text = element_blank(),
          axis.title = element_blank())
dev.off()

jpeg(paste0("xy_persample.jpg"), width=800*3, height=700*3, res=200)

ggplot(cluster_info, aes(x = x, y = y, color = cluster)) +
    #ggrastr::geom_point_rast(size = 0.01, alpha = 0.8) +
    geom_point(size=0.001)+
   # defined_theme +
    facet_wrap(~cluster,ncol=3)+
    #scale_color_manual(values = my_colors)+
    # guides(color = guide_legend(override.aes = list(size = 3),
    #                             ncol = 1, title = "")) +
    labs(title = " ", x = "UMAP1", y = "UMAP2", fill=" ") +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 20, face = "plain"),
          legend.position = "right",
          legend.text = element_text(size = 12),
          panel.border = element_rect(colour = "black",
                                      fill = NA, linewidth = 0.3),
          axis.line = element_blank(),
          legend.justification = "center",
          axis.ticks = element_blank(),
          legend.box.just = "center",
          axis.text = element_blank(),
          axis.title = element_blank())
dev.off()



# define low library cells (example cutoff = 0 or < 300 UMIs)
low_cutoff <- 50
low_cells <- lib_size < low_cutoff

# get coordinates
coords <- spatialCoords(se)

# build plotting dataframe
df <- data.frame(
    x = coords[,1],
    y = coords[,2],
    lib_size = lib_size,
    low = low_cells
)

jpeg(paste0("xy_persample.jpg"), width=1000, height=1000, res=200)
# scatter plot with highlighting
ggplot(df[df$low==TRUE,], aes(x = x, y = y, color = low)) +
    geom_point(size = 0.0001) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    coord_equal() +
    theme_minimal() +
    labs(color = paste0("lib_size <", low_cutoff))

dev.off()