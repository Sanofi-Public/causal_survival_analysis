# Causal Survival Analysis Project (R + renv)

This project accompanies the paper:

> *Causal survival analysis combines survival analysis and causal inference to evaluate the effect of a treatment or intervention on a time-to-event outcome, such as survival time. It offers an alternative to relying solely on Cox models for assessing these effects. In this paper, we present a comprehensive review of estimators for the average treatment effect measured with the restricted mean survival time, including regression-based methods, weighting approaches, and hybrid techniques. We investigate their theoretical properties and compare their performance through extensive numerical experiments. Our analysis focuses on the finite-sample behavior of these estimators, the influence of nuisance parameter selection, and their robustness and stability under model misspecification. By bridging theoretical insights with practical evaluation, we aim to equip practitioners with both state-of-the-art implementations of these methods and practical guidelines for selecting appropriate estimators for treatment effect estimation. Among the approaches considered, G-formula two-learners, AIPCW-AIPTW, Buckley-James estimators, and causal survival forests emerge as particularly promising.*


## Authors

- [Charlotte Voinot](https://chvoinot.github.io/) – [INRIA PreMeDICaL](https://team.inria.fr/premedical/), [Université de Montpellier](https://www.umontpellier.fr/), [Sanofi R&D](https://www.sanofi.com/en)
- [Clément Berenfeld](https://cberenfeld.github.io/) – [Universität Potsdam](https://www.uni-potsdam.de/en/university-of-potsdam)
- Imke Mayer – [Charité – Universität Berlin](https://www.charite.de/)
- Bernard Sebastien – [Sanofi R&D](https://www.sanofi.com/en)
- [Julie Josse](http://juliejosse.com/) – [INRIA PreMeDICaL](https://team.inria.fr/premedical/), [Université de Montpellier](https://www.umontpellier.fr/)

---

## Quickstart

This project uses [`renv`](https://rstudio.github.io/renv/) to manage R package dependencies in an isolated and reproducible environment.

### 1. Open the project in R or RStudio

Make sure you're in the **root directory** of the project, where the `renv.lock`, `.Rprofile`, and `renv/` folder are located.

### 2. Automatic activation

The `.Rprofile` file contains:

```r
source("renv/activate.R")
``` 

This ensures that the correct project environment is automatically activated when opening the project in R or RStudio. You don’t need to manually call renv::activate().


### 3. Activate the renv environment(if needed)

This tells R to use the project-specific libraries instead of your global ones:

```r
install.packages("renv")
renv::activate()
```

### 4. Restore packages from renv.lock

To install exactly the right versions of all required packages:

```r
renv::restore()
```

This may take a few minutes the first time.

##  Project Structure

```bash
.
├── renv/                               # renv environment folder
├── renv.lock                           # lockfile with package versions
├── Notebook_causal_survival.qmd        # main analysis script
├── utilitary.R                         # main analysis script
├── README.md                           # this file
└── simulations/                        # .RData simulations 
```

