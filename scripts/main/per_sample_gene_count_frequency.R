library(jazzPanda)
library(Seurat)
library(data.table)
library(dplyr)
library(glmnet)
library(caret)
library(limma)
library(edgeR)
library(speckle)
library(peakRAM)
library(dplyr) 
library(scCustomize)
library(SpatialExperimentIO)
library(SpatialExperiment)
library(here)
library(Matrix)
source(here("scripts/utils.R"))
args_all   = commandArgs(trailingOnly = FALSE)
script_arg = grep("^--file=", args_all, value = TRUE)
script_dir = if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())
}
project_root = normalizePath(file.path(script_dir, "..", ".."))


# Per-value frequency table for one count matrix (dense base matrix or sparse Matrix).
# Zero bin computed as prod(dim) - nnz (double) so it never overflows the integer range.
count_freq = function(cm, dataset_nm){
    if (is(cm, "sparseMatrix")){
        cm = drop0(cm)
        nz = cm@x
    } else {
        nz = cm[cm != 0]
    }
    total  = prod(dim(cm))
    counts = tabulate(nz + 1L)
    if (length(counts) == 0L) counts = 0
    counts[1] = total - length(nz)
    names(counts) = 0:(length(counts) - 1L)
    counts = counts[counts > 0]
    prop   = counts / sum(counts)
    data.frame(
        sample = dataset_nm,
        value  = names(counts),
        count  = as.numeric(counts),
        prop   = round(prop, 4),
        label  = paste0(format(counts, big.mark = ",", scientific = FALSE, trim = TRUE),
                        " (", round(100 * prop, 4), "%)"),
        row.names = NULL
    )
}

summary_tab = data.frame()

# ---- Xenium human breast ----
sp1_path = "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_breast_samples/Xenium_hbreast_sample1/"
sp1 = get_xenium_data(sp1_path, mtx_name = "cell_feature_matrix",
                      trans_name = "transcripts.csv.gz", cells_name = "cells.csv.gz")
sp1$trans_info = sp1$trans_info[sp1$trans_info$qv >= 20 &
                                    sp1$trans_info$cell_id != -1 &
                                    !(sp1$trans_info$cell_id %in% sp1$zero_cells), ]

sp2_path = "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_breast_samples/Xenium_hbreast_sample2/"
sp2 = get_xenium_data(sp2_path, mtx_name = "cell_feature_matrix",
                      trans_name = "transcripts.csv.gz", cells_name = "cells.csv.gz")
sp2$trans_info = sp2$trans_info[sp2$trans_info$qv >= 20 &
                                    sp2$trans_info$cell_id != -1 &
                                    !(sp2$trans_info$cell_id %in% sp2$zero_cells), ]

summary_tab = rbind(summary_tab, count_freq(sp1$cm, "xenium_hbreast_sp1"))
summary_tab = rbind(summary_tab, count_freq(sp2$cm, "xenium_hbreast_sp2"))

# ---- MERSCOPE human breast ----
se = readMerscopeSXE(dirName = "/stornext/Bioinf/data/lab_phipson/givanna/merscope_data/HumanBreastCancerPatient1/",
                     countMatPattern = "cell_by_gene.csv", metaDataPattern = "cell_metadata.csv")
x_avg  = (se@colData$min_x + se@colData$max_x) / 2
y_avg  = (se@colData$min_y + se@colData$max_y) / 2
coords = cbind(x = x_avg, y = y_avg)
SpatialExperiment::spatialCoords(se) = coords

se = se[-grep("^Blank", rownames(se)), ]
lib_size = Matrix::colSums(assay(se, "counts"))
colData(se)$lib_size = lib_size
se = se[, lib_size > 30 & lib_size < 2500]

seu = as.Seurat(se, data = NULL)
cm  = GetAssayData(seu, assay = "originalexp", layer = "counts")
summary_tab = rbind(summary_tab, count_freq(cm, "merscope_sp1"))

# ---- CosMx liver (normal + cancer) ----
seu = readRDS("/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/CosMx_normal_and_diseased_liver_samples/LiverDataReleaseSeurat_newUMAP.RDS")
metadata = as.data.frame(seu@meta.data)
qc_cols = c("qcFlagsRNACounts", "qcFlagsCellCounts", "qcFlagsCellPropNeg",
            "qcFlagsCellComplex", "qcFlagsCellArea", "qcFlagsFOV")
metadata = metadata[!apply(metadata[, qc_cols], 1, function(x) any(x == "Fail")), ]
cell_info_cols = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm",
                   "nCount_negprobes", "nFeature_negprobes", "nCount_falsecode",
                   "nFeature_falsecode", "slide_ID_numeric", "Run_Tissue_name",
                   "fov", "cellType", "niche", "cell_id")
cellCoords = metadata[, cell_info_cols]
counts = seu[["RNA"]]@counts

normal_cells = row.names(cellCoords[cellCoords$slide_ID_numeric == 1, ])
cancer_cells = row.names(cellCoords[cellCoords$slide_ID_numeric == 2, ])

summary_tab = rbind(summary_tab, count_freq(counts[, normal_cells], "cosmx_healthy"))
summary_tab = rbind(summary_tab, count_freq(counts[, cancer_cells], "cosmx_cancer"))

# ---- Xenium human lung ----
hl_cancer = get_xenium_data(path = "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_human_lung_cancer/",
                            mtx_name = "cell_feature_matrix/", trans_name = "transcripts.csv.gz",
                            cells_name = "cells.csv.gz")
summary_tab = rbind(summary_tab, count_freq(hl_cancer$cm, "xenium_lung"))

# ---- Xenium mouse brain (3 sections) ----
mbrain_base = "/vast/projects/xenium_5k/data/jazzPanda_paper_dataset/Xenium_mouse_brain/"
mbrain_dirs = c(xenium_mbrain_sp1 = "Xenium_V1_FF_Mouse_Brain_MultiSection_1_outs/",
                xenium_mbrain_sp2 = "Xenium_V1_FF_Mouse_Brain_MultiSection_2_outs/",
                xenium_mbrain_sp3 = "Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/")
for (nm in names(mbrain_dirs)){
    rep = get_xenium_data(path = paste0(mbrain_base, mbrain_dirs[[nm]]),
                          mtx_name = "cell_feature_matrix/", trans_name = "transcripts.csv.gz",
                          cells_name = "cells.csv.gz")
    summary_tab = rbind(summary_tab, count_freq(rep$cm, nm))
}

out_dir = file.path(project_root, "data")
out_dir = normalizePath(out_dir)
write.csv(summary_tab,  file.path(out_dir,"gene_count_frequency.csv"),
          row.names = FALSE)

#################


df = fread("/vast/projects/xenium_5k/jazzPanda_paper/data/gene_count_frequency.csv")
df = as.data.frame(df)
df$value = as.integer(df$value)


keep_vals = 0:5
val_df = df %>% filter(value %in% keep_vals)

p_val = ggplot(val_df, aes(factor(value), prop, fill = factor(value))) +
    geom_col() +
    geom_text(aes(label = percent(prop, accuracy = 0.1)),
              vjust = -0.3, size = 2.2) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.15))) +
    facet_wrap(~ sample, ncol = 5) +
    labs(x = "Gene count", y = "Proportion of gene count", fill="gene count") +
    
    theme_bw()+
    theme(strip.background = element_rect(fill=NA, color="black"))

ggsave("fig_gene_count_proportion_0_5.png", p_val, width = 12, height = 5, dpi = 200)


cap = 10
dist_df = df %>%
    mutate(value_grp = ifelse(value >= cap, cap, value)) %>%
    group_by(sample, value_grp) %>%
    summarise(count = sum(count), .groups = "drop")
lab_fun = function(v){
    v = as.integer(as.character(v))
    ifelse(v >= cap, paste0(cap, "+"), as.character(v))
}

p_dist = ggplot(dist_df, aes(factor(value_grp), count)) +
    geom_col(fill = "steelblue") +
    scale_x_discrete(labels = lab_fun) +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
    facet_wrap(~ sample, scales = "free_y") +
    labs(x = "Count value", y = "Number of entries (log scale)",
         title = "Count values collapse onto 0 and 1") +
    theme_minimal(base_size = 11)

dist_df = df %>%
    mutate(value_grp = ifelse(value >= cap, cap, value)) %>%
    group_by(sample, value_grp) %>%
    summarise(prop = sum(prop), .groups = "drop")

p_dist = ggplot(dist_df, aes(factor(value_grp), prop)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = percent(prop, accuracy = 0.1)), vjust = -0.3, size = 2) +
    scale_x_discrete(labels = lab_fun) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ sample) +
    labs(x = "Gene count", y = "Proportion of gene count", fill="gene count") +
    theme_bw()+
    theme(strip.background = element_rect(fill=NA, color="black"))


ggsave("fig_count_distribution.png", p_dist, width = 11, height = 7, dpi = 200)
# ggsave("fig_nonzero_distribution.png", p_nz, width = 11, height = 7, dpi = 200)










# pdf(paste(dataset_nm,"_gene_expression_frequency.pdf"), width=12, height=8)
# ggplot(summary_tab, aes(x = value_num, y = count)) +
#     geom_col(fill = "steelblue", width = 0.9) +
#     geom_text(data = lab_df, aes(label = short),
#               vjust = -0.4, size = 2.6) +
#     facet_zoom(x = value_num %in% zoom_lab$value) +
#     scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
#     labs(x = "Value", y = "Count", title = "Distribution of count values") +
#     theme_minimal(base_size = 12)
# dev.off()
# 

