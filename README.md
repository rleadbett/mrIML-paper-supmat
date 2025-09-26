# MrIML 2.0 Case Studies

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Quarto](https://img.shields.io/badge/Made%20with-Quarto-blue.svg)](https://quarto.org/)

Supplementary code for: **Fountain-Jones, N.M., Appaw, R.C., et al.** (2024). Advancing ecological community analysis with MrIML 2.0: Unravelling taxa associations through interpretable machine learning. *Molecular Ecology Resources*.

This repository contains the code to reproduce the simulation study and the three case studies:

- Avian Parasites
- Forest Communities
- Ostrich Microbiome

## Quick Start

The supplementary material is organised in a quarto document. To reproduce the work, start by cloning the repository and reproducing the R enviroment with [`renv`](https://rstudio.github.io/renv/articles/renv.html).

```bash
git clone https://github.com/yourusername/mrIML-case-studies.git
cd mrIML-case-studies
R -e "renv::restore()"
```

To render the case studies, run

## Case Studies

1. **[Avian Parasites](case-studies/01-avian-parasites/)** - 449 birds, 4 parasites (~30 min)
2. **[Forest Communities](case-studies/02-forest-communities/)** - Multi-taxa forest analysis (~2 hrs) 
3. **[Ostrich Microbiome](case-studies/03-ostrich-microbiome/)** - Gut microbiome development (~4 hrs)
4. **[Data Curation](case-studies/04-microbiome-curation/)** - Microbiome preprocessing (~15 min)

## Requirements

- R ≥ 4.1.0
- See `renv.lock` for package versions

## Contact

nick.fountainjones@utas.edu.au