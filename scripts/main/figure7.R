# calculate cluster-cluster correlation with spatial vectors 
library(ggplot2)
library(here)
library(data.table)
library(jazzPanda)
library(dplyr)
library(xtable)
source(here("scripts/utils.R"))

data_nm  <- "merscope_hbreast"
cluster_info = readRDS(here(data_path,paste0(data_nm, "_clusters.Rds")))
cluster_names = paste0("c", 1:10)
anno_df <- data.frame(
    cluster = cluster_names,
    major_class = c(
        "Epithelial tumor","Stromal","Epithelial tumor","Epithelial tumor",
        "Myeloid","Myeloid","Lymphoid","Endothelial / Vascular","B lineage","Dendritic cell"
    ),
    sub_class = c(
        "Luminal-like / ERBB2-high","Cancer-associated fibroblast (CAF)","IFN-γ–licensed","Cycling (G2/M)",
        "Tumor-associated macrophage (TAM)","Inflammatory/angiogenic TAM","T cells (activated)","Blood endothelial","Plasma cell / plasmablast","Langerhans / cDC2-like"
    ),
    cell_type = c(
        "Carcinoma—proliferative ERBB2/EGFR⁺",
        "myCAF / iCAF hybrid",
        "Carcinoma—MHC-II⁺ (APC-like)",
        "Carcinoma—G2/M cycling",
        "C1QC⁺ M2-like TAM",
        "SPP1⁺ (osteopontin) TAM",
        "Mixed CD4/CD8 with Treg features",
        "Endothelium—PLVAP/KDR⁺",
        "Plasma cells",
        "CD207⁺ APC"
    ),
    supporting_genes = c(
        "EPCAM, ERBB2, ERBB3, EGFR, CDH1, MYC, FOXM1, MYBL2, MCM2, PCNA, CCNB1, LGR5, DKK1, FZD7, EPHB3",
        "FN1, COL1A1, COL11A1, COL6A3, COL5A1, ACTA2, PDGFRA, TNC, FAP, CXCL12, IL6, SERPINE1, MFAP5, ELN, LOX",
        "CDH1, EPCAM, HLA-DPB1/DRB1, CIITA, IDO1, TNFSF10, TGFB2, CEACAM1",
        "BIRC5, PLK1, CCNB1, AURKA, AURKB, FOXM1, MYBL2, MCM2, TP53",
        "C1QC, LYZ, FCGR3A, CD14, MRC1, CD163, MSR1, CSF1R, TREM2, LGALS9",
        "SPP1, VEGFA, MMP9, MMP12, CXCL8, SERPINE1, TGFBI, ICAM1, ATF3",
        "TRAC, CD3D/E/G, CD2, CD8A, GZMA/H/K, NKG7, CXCR3, TIGIT, ICOS, FOXP3, CCR7, TBX21, EOMES",
        "PECAM1, VWF, CDH5, KDR, FLT1, CLDN5, PLVAP, MMRN2, ANGPT2, PGF, CLEC14A",
        "XBP1, MZB1, IRF4, POU2AF1, CD79A/B, DERL3, FCRL5",
        "CD207, FCER1A, CD1C, CD1B/E, CCL22, CSF1R"
    ),
    notes = c(
        "Luminal-like/HER2-high epithelial tumor with strong proliferation and WNT cues.",
        "Mixed myCAF/iCAF phenotype; ECM deposition and inflammatory signaling coexist.",
        "IFN-γ–induced antigen presentation by tumor cells; epithelial identity preserved.",
        "Mitotic (G2/M) program dominates; classic cycling-tumor signature.",
        "Immunoregulatory TAMs with M2-like polarization; chemokine-rich.",
        "Angiogenic/remodeling TAMs; SPP1/VEGFA/MMP axis.",
        "Activated T cells with cytotoxic and Treg subsets co-enriched.",
        "Activated tumor endothelium; tip/stalk features (PLVAP/ANGPT2).",
        "Canonical plasma cell program; minor contamination unlikely to change call.",
        "Langerin+ dendritic cells bridging cDC2/Langerhans features."
    ),
    confidence = c("High","High","High","High","High","High",
                   "Medium","High","High","High"),
    
    anno_name = c("Tumor_Luminal_ERBB2", 
                  "CAF_Mixed",             
                  "Tumor_IFNg_APC",        
                  "Tumor_Cycling_G2M",     
                  "TAM_C1QC",             
                  "TAM_SPP1",              
                  "Tcell_Activated",       
                  "Endo_Blood",            
                  "Plasma",           
                  "DC_Langerhans_cDC2"),
    stringsAsFactors = FALSE
)


cluster_info = merge(cluster_info, anno_df, by="cluster")

cluster_info$cluster = factor(cluster_info$cluster,
                              levels=cluster_names)

cluster_info$anno_name = factor(cluster_info$anno_name,
                                levels=anno_df$anno_name)
cluster_info = cluster_info[order(cluster_info$anno_name), ]


jazzPanda_res_lst = readRDS(here(data_path,paste0(data_nm, "_jazzPanda_res_lst.Rds")))

perm_lst = readRDS(here(data_path,paste0(data_nm, "_perm_lst.Rds")))
perm_res = get_perm_adjp(perm_lst)
obs_corr = get_cor(perm_lst)


sv_lst = readRDS(here(data_path,paste0(data_nm, "_sq50_vector_lst.Rds")))
nbins = 2500

exp_ord <- anno_df$anno_name
cluster_mt = sv_lst$cluster_mt[, paste0("c", 1:10)]

# cluster_mt <-cluster_mt[,exp_ord]

colnames(cluster_mt) = anno_df$anno_name
cor_cluster_mt <- cor(cluster_mt,cluster_mt, method = "pearson")

col <- grDevices::colorRampPalette(c("#4477AA", "#77AADD", 
                                     "#FFFFFF","#EE9988", "#BB4444"))

pdf(file.path(fig7, "figure7_clustercluster_corr.pdf"), width=10, height=10)
corrplot::corrplot(cor_cluster_mt, method="color", col=col(200), diag=TRUE,
                   addCoef.col = "black",
                   type="upper",cl.cex=1,
                   tl.col="black", tl.srt=45, mar=c(0,0,5,0),sig.level = 0.05, 
                   insig = "blank", 
                   title = " "
)

dev.off()
# 
# selected_mg = get_top_mg(jazzPanda_res_lst) %>%
#     filter(top_cluster != "NoSig") %>%
#     mutate(top_cluster = factor(top_cluster, levels = paste0("c", 1:10))) %>%
#     group_by(top_cluster) %>%
#     slice_max(order_by = glm_coef, n = 3) %>%
#     arrange(top_cluster, desc(glm_coef)) %>%
#     select(gene, top_cluster, glm_coef) %>% data.frame()
# 
# 
# # Calculate pairwise correlations
# cor_gene_mt <- cor(sv_lst$gene_mt[, selected_mg$gene], sv_lst$gene_mt[, selected_mg$gene],
#                    method = "pearson")
# xtable(cor_gene_mt)
# pdf(file.path(fig7, "figure7_genegene_corr_glm.pdf"), width=16, height=16)
# corrplot::corrplot(cor_gene_mt, method="color", col=col(200), diag=TRUE,
#                    addCoef.col = NULL,type="upper",
#                    tl.col="black", tl.srt=45, mar=c(0,0,5,0),sig.level = 0.05,
#                    insig = "blank",
#                    title = " ",cl.cex=1.5
# )
# dev.off()
selected_mg <- unlist(lapply(cluster_names, function(cl) {
    
    obs_cutoff <- quantile(obs_corr[, cl], 0.75)
    
    # genes passing both permutation and observed-corr thresholds
    perm_cl <- intersect(
        rownames(perm_res)[perm_res[, cl] < 0.05],
        rownames(obs_corr)[obs_corr[, cl] > obs_cutoff]
    )
    
    if (length(perm_cl) == 0) return(NULL)
    
    df <- data.frame(
        gene  = perm_cl,
        value = signif(obs_corr[perm_cl, cl], 3),
        stringsAsFactors = FALSE
    )
    
    df <- df[order(df$value, decreasing = TRUE), ]
    
    # return top 3 genes for this cluster
    head(df$gene, 3)
}))


# Calculate pairwise correlations
cor_gene_mt <- cor(sv_lst$gene_mt[, unique(selected_mg)],
                   sv_lst$gene_mt[, unique(selected_mg)],
                   method = "pearson")
xtable(cor_gene_mt)
pdf(file.path(fig7, "figure7_genegene_corr_perm.pdf"), width=16, height=16)
corrplot::corrplot(cor_gene_mt, method="color", col=col(200), diag=TRUE,
                   addCoef.col = NULL,type="upper",
                   tl.col="black", tl.srt=45, mar=c(0,0,5,0),sig.level = 0.05,
                   insig = "blank",
                   title = " ",cl.cex=1.5
)
dev.off()
