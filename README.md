# jazzPanda_paper
This repository provides the complete set of **R scripts** used to generate **all main and supplementary figures** for the *jazzPanda* paper.  

- **jazzPanda R package**: [GitHub Repository](https://github.com/phipsonlab/jazzPanda)  
- **Preprint**: [bioRxiv link](https://doi.org/XXXXXX)  
- **Raw datasets**: [Download from Zenodo](https://zenodo.org/record/XXXXXX)  
- **Analysis workflow**: [Workflow documentation](https://github.com/phipsonlab/jazzPanda_workflow)  

Processed outputs are organized in the `data/` directory, and large intermediate results are saved as `.Rds` files.

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

## Figure map 

| Figure    | Script                      | Output                       |
|-----------|-----------------------------|------------------------------|
| Figure 2  | `scripts/main/figure2.R`       | `figures/main/figure2_cosmx_hliver/`      |
| Figure 3  | `scripts/main/figure3.R`       | `figures/main/figure3_xenium_hbreast/`      |
| Figure 4  | `scripts/main/figure4.R`       | `figures/main/figure4_simulation/`      |
| Figure 5  | `scripts/main/figure5.R`       | `figures/main/figure5_compare_methods/`      |
| Figure 6  | `scripts/main/figure6.R`       | `figures/main/figure6/`      |


## Script Overview

### Marker gene: 

| Location       | Scripts                                                     | Output                                                      |
| -------------- | ----------------------------------------------------------- | ------------------------------------------------------------|
| `scripts/main` | `run_mg_xenium_mouse_brain.sh`<br>`mg_xenium_mouse_brain.R` | `data/dataset_computational_complexity/xenium_mbrain_*.Rds` |
| `scripts/main` | `run_mg_xenium_human_breast_cancer.sh`<br>`mg_xenium_human_breast_cancer.R` | `data/dataset_computational_complexity/xenium_hbreast_*.Rds` |
| `scripts/main` | `run_mg_xenium_human_breast_cancer.sh`<br>`mg_xenium_human_breast_cancer.R` | `data/dataset_computational_complexity/xenium_hlc_*.Rds` |
| `scripts/main` | `run_mg_merscope_human_breast_cancer.sh`<br>`mg_merscope_human_breast_cancer.R` | `data/dataset_computational_complexity/merscope_hbreast_*.Rds` |
| `scripts/main` | `run_mg_cosmx_human_healthy_liver.sh`<br>`mg_cosmx_human_healthy_liver.R` | `data/dataset_computational_complexity/cosmx_hhliver_*.Rds` |
| `scripts/main` | `run_mg_cosmx_human_liver_cancer.sh`<br>`mg_cosmx_human_liver_cancer.R` | `data/dataset_computational_complexity/cosmx_hlc_*.Rds` |

### Simulation
| **Analysis**                              | **Location**      | **Scripts**                                                                           | **Output**                                   |
|-------------------------------------------|-------------------|---------------------------------------------------------------------------------------|----------------------------------------------|
| Simulation: CosMx human liver cancer      | `scripts/main`    | `cosmx_hlc_simulation_simbg_sa.sh`<br>`cosmx_hlc_simulation_simbg_slurmarray_temp.R`  | `scripts/main/cosmx_hlc_simulation_result/`  |


### Technical performance: 
#### Tiles on marker genes 
| Location                                        | Scripts                                                              | Output                                                                               |
| ----------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| `scripts/main/discussion_markergenes_vs_ntiles` | `cosmx_ntile_markergenes_sa.sh`<br>`cosmx_ntile_markergene_result.R` | `scripts/main/discussion_markergenes_vs_ntiles/cosmx_hlc_ntiles_mg_gr{10–100}_*.csv` |
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
