# jazzPanda_paper
This repository contains the **exact R scripts** used to generate **all main figures** and **all supplementary figures** for the jazzPanda paper.

## Figure map 

| Figure    | Script                      | Output                       |
|-----------|-----------------------------|------------------------------|
| Figure 2  | `scripts/main/fig2.R`       | `figures/main/fig2.pdf`      |

```sh
.
├── scripts/
│   ├── main/                 
│   └── supp/                 
├── figures/                  # outputs (PDF/PNG/SVG) — created by scripts
│   ├── main/
│   └── supp/
├── configs/
│   └── paths.yml             # local paths to data dirs (edit me)`
├── data/
│   ├── README.md             # where/how to obtain datasets
│   ├── raw/                  # large downloads (gitignored)
│   └── processed/            # derived caches (gitignored)

├── renv.lock                 # pinned R package versions
├── .gitignore
└── README.md                 # you are here
```
