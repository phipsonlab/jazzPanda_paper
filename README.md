# jazzPanda: Paper Repository
This repository provides the complete set of **R scripts** used to generate **all main and supplementary figures** for the *jazzPanda* paper.  

- **jazzPanda R package**: [Repository](https://bioconductor.org/packages/jazzPanda)  
- **Preprint**: [bioRxiv link](https://doi.org/10.64898/2026.02.13.705867)  
- **Raw/Processed data**: [Download from Zenodo](https://zenodo.org/records/18149456)  
- **Analysis workflow**: [Workflow documentation](https://phipsonlab.github.io/jazzPanda_workflowr/)  

# jazzPanda: Paper Repository
This repository provides the complete set of **R scripts** used to generate **all main and supplementary figures** for the *jazzPanda* paper.

- **jazzPanda R package**: [Repository](https://bioconductor.org/packages/jazzPanda)
- **Preprint**: [bioRxiv link](https://doi.org/10.64898/2026.02.13.705867)
- **Raw/Processed data**: [Download from Zenodo](https://zenodo.org/records/18149456)
- **Analysis workflow**: [Workflow documentation](https://phipsonlab.github.io/jazzPanda_workflowr/)

```sh
.
├── scripts/
│   ├── main/
│   └── supp/
├── figures/                  # outputs (PDF/JPG) — created by scripts
│   ├── main/
│   └── supp/
├── data/
├── .gitignore
└── README.md                 # you are here
```

## Datasets

Four imaging-based spatial transcriptomics datasets are used throughout the paper:

| Short name | Platform | Tissue |
|---|---|---|
| `xenium_hbreast` | 10x Xenium | Human HER2+ breast cancer (2 samples) |
| `cosmx_hhliver`  | Nanostring CosMx | Human healthy liver |
| `cosmx_hlc`      | Nanostring CosMx | Human liver cancer |
| `merscope_hbreast` | Vizgen MERSCOPE | Human breast cancer |

## Data Availability
Raw experimental data can be accessed via the links below. Processed outputs used in the manuscript are included in the repository under the `data/` directory, with intermediate results stored as `.Rds` files. For very large processed datasets that exceed the repository size limits, we provide access via **Zenodo**.

### Raw data
-1) 10x Xenium Human HER2+ breast cancer data [🔗](https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast)

-2) Nanostring CosMx Human liver healthy and cancer data [🔗](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/)

-3) Vizgen MERSCOPE Human breast cancer data [🔗](https://info.vizgen.com/ffpe-showcase?submissionGuid=c9a25730-3fe5-444b-bef8-1a74d51ddefb)

### Processed data
Processed outputs are organized in the `data/` directory of this repository, with intermediate results saved as `.Rds` files. Very large processed datasets that exceed repository size limits are deposited on [Zenodo](https://zenodo.org/records/18149456) for convenient download.

## Figure map

| Figure | Content | Script | Output |
|---|---|---|---|
| Figure 1 | Data sparsity and performance of single-cell marker analysis methods | `scripts/main/figure_gene_count.R` | `figures/main/figure_intro/` |
| Figure 2 | Method overview (schematic) | — | — |
| Figure 3 | Negative control probes and simulation-based evaluation | `scripts/main/figure_simulation.R` | `figures/main/figure_simulation/` |
| Figure 4 | Application to CosMx healthy liver and Xenium HER2+ breast cancer | `scripts/main/figure_result.R` | `figures/main/figure_result/` |
| Figure 5 | Comparison of marker analysis methods | `scripts/main/figure_compare_methods.R` | `figures/main/figure_compare_methods/` |
| Figure 6 | Extension of the jazzPanda framework (MERSCOPE breast cancer) | `scripts/main/figure_sv_extension.R` | `figures/main/figure_sv_extension/` |
| Figure 7 | Technical performance (bin size, runtime, memory) | `scripts/main/figure_technical_performance.R` | `figures/main/figure_technical_performance/` |

## Script Overview

### Marker gene detection

| Location       | Scripts                                                     | Output                                                      |
| -------------- | ----------------------------------------------------------- | ------------------------------------------------------------|
| `scripts/main` | `run_mg_xenium_human_breast_cancer.sh`<br>`mg_xenium_human_breast_cancer.R` | `data/dataset_computational_complexity/xenium_hbreast_*.Rds` |
| `scripts/main` | `run_mg_cosmx_human_healthy_liver.sh`<br>`mg_cosmx_human_healthy_liver.R` | `data/dataset_computational_complexity/cosmx_hhliver_*.Rds` |
| `scripts/main` | `run_mg_cosmx_human_liver_cancer.sh`<br>`mg_cosmx_human_liver_cancer.R` | `data/dataset_computational_complexity/cosmx_hlc_*.Rds` |
| `scripts/main` | `run_mg_merscope_human_breast_cancer.sh`<br>`mg_merscope_human_breast_cancer.R` | `data/dataset_computational_complexity/merscope_hbreast_*.Rds` |

### Simulation
| **Analysis**                              | **Location**      | **Scripts**                                                                           | **Output**                                   |
|-------------------------------------------|-------------------|---------------------------------------------------------------------------------------|----------------------------------------------|
| Simulation: CosMx human liver cancer      | `scripts/main`    | `cosmx_hlc_simulation_simbg_sa.sh`<br>`cosmx_hlc_simulation_simbg_slurmarray_temp.R`  | `scripts/main/cosmx_hlc_simulation_result/`  |

### Technical performance
#### Tiles on marker genes
| Location                                        | Scripts                                                              | Output                                                                               |
| ----------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| `scripts/main/discussion_markergenes_vs_ntiles` | `xenium_hbreast_ntile_markergenes_sa.sh`<br>`xenium_hbreast_ntile_markergene_result.R` | `scripts/main/discussion_markergenes_vs_ntiles/xenium_hbreast_ntiles_mg_gr{10–100}_*.csv` |

#### Computational complexity on number of cores
| Location                             | Scripts                                                                                | Output                                              |
| ------------------------------------ | -------------------------------------------------------------------------------------- | --------------------------------------------------- |
| `scripts/main/discussion_complexity` | `complexity_ncores_sa.sh`<br>`complexity_ncores_cosmx_hliver_cancer_slurmarray_temp.R` | `scripts/main/discussion_complexity/ncores_result/` |

#### Computational complexity on number of transcripts
| Location                             | Scripts                                                     | Output                                              |
| ------------------------------------ | ----------------------------------------------------------- | --------------------------------------------------- |
| `scripts/main/discussion_complexity` | `complexity_ntr_sim_sa.sh`<br>`complexity_ntr_simulation.R` | `scripts/main/discussion_complexity/ngenes_result/` |

#### Computational complexity on number of tiles
| Location                             | Scripts                                                                                | Output                                                         |
| ------------------------------------ | -------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| `scripts/main/discussion_complexity` | `complexity_ntiles_squarebins_sim_sa.sh`<br>`complexity_ntiles_squarebin_simulation.R` | `scripts/main/discussion_complexity/ntiles_squarebins_result/` |
| `scripts/main/discussion_complexity` | `complexity_ntiles_hexbin_sim_sa.sh`<br>`complexity_ntiles_hexbin_simulation.R` | `scripts/main/discussion_complexity/ntiles_hexbin_result/` |


## Supplementary material

Supplementary figures and tables are generated as R Markdown notebooks in `scripts/supp/`, each rendered to a self-contained `.html` of the same name.

| Notebook | Content |
|---|---|
| `supplementary_simulation.Rmd` | Permutation p-values and correlation of simulated genes across clusters; false discovery rates |
| `supplementary_application_cosmx_healthy_human_liver.Rmd` | Dataset overview; top marker genes per cluster and cluster–gene vector relationships (jazzPanda-correlation and jazzPanda-glm); marker gene overlap and cumulative average correlation across methods |
| `supplementary_application_cosmx_human_liver_cancer.Rmd` | Dataset overview |
| `supplementary_application_xenium_human_breast_cancer.Rmd` | Dataset overview; top marker genes per cluster and cluster–gene vector relationships (jazzPanda-glm); marker gene overlap and cumulative average correlation across methods; marker gene overlap across tile lengths |
| `supplementary_application_merscope_human_breast_cancer.Rmd` | Dataset overview; spatial visualisation of selected cell types|
| `supplementary_technical_performance.Rmd` | Computational complexity of spatial vector construction across square and hex bin lengths; runtime and memory benchmarks |

