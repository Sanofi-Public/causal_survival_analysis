# Causal Survival Analysis – Supplementary Material

This repository provides all code, data, and instructions necessary to reproduce the results of the manuscript:

> **Treatment Effect Estimation in Causal Survival Analysis: Practical Recommendations**  
> Charlotte Voinot – INRIA PreMeDICaL, Université de Montpellier, Sanofi R&D  
> Clément Berenfeld – Universität Potsdam  
> Imke Mayer – Charité – Universität Berlin  
> Bernard Sebastien – Sanofi R&D  
> Julie Josse – INRIA PreMeDICaL, Université de Montpellier  
> *Submitted to Biometrical Journal, 2025*

---

## Repository Structure

This repository is organized as follows:

```
.
├── code/                      # Analysis and simulation code
│   ├── main_simulations.qmd           # Master script to run all simulations
│   ├── simulate_data.R                # Functions to generate simulated data
│   ├── estimators/                    # Estimator-specific folder
│   │   ├── Estimators.R               # Functions of the estimators
│   │   └── utilitary.R                # utilitary function to compute the estimators
│   ├── results/                       # Precomputed simulation outputs used in the article
│   ├── figures/                       # Figures saved in .pdf
│   └── user_runs/                     # Folder to store user-generated simulations
├── README.md
├── renv.lock                          # R environment snapshot (full reproducibility)
├── LICENSE.txt

```

---

## Reproducing the Results


To support full computational reproducibility, we provide:

- **`sessionInfo.txt`**: output from `sessionInfo()` after loading only the packages needed to render figures. This reflects the minimal runtime environment.

- **`renv.lock`**: full snapshot of the R environment used for the project, including all packages required to run the simulations and fit estimators. Recommended if you wish to reuse the estimator implementations.

- **`installed_packages.txt`**: complete list of installed packages (with dependencies) at the time of rendering.


### Requirements

- **R version**: 4.2.x or higher  
- Recommended: RStudio  
- All packages are managed via `renv`

To restore the environment, run:

```{r}
renv::activate()
renv::restore()
```

If `quarto` is not installed:

```{r}
install.packages("quarto")
install.packages("this.path")
```

### Run the Notebook

**Note**: Running all simulations from scratch can be time-consuming and resource-intensive, due to the large number of estimators to compute — each evaluated over 100 repetitions and in both parametric and non-parametric settings. We strongly recommend using a powerful virtual machine or computing server to reproduce the full set of results.
By default, the notebook loads precomputed simulations due to the very long time of computation (e.g., `load("results/simulation_rct1.RData")`). 

1. To generate the figures and results from precomputed simulations:

```{r}
quarto::quarto_render("code/main_simulations.qmd", execute_params = list(run_simulations = FALSE))
```

This will use the saved results to generate figures and tables


2. To generate the figures and results from scratch: 

```{r}
quarto::quarto_render("code/main_simulations.qmd", execute_params = list(run_simulations = TRUE))
```

This will run simulations under all defined scenarios.

At the end of the process, both HTML and PDF reports will be generated, containing the corresponding results. Additionally, all figures will be saved as individual .pdf files in the figures/ folder.

*Note*: Figures 8 to 14 in the main text are identical to Figures 15 to 21, except that they display a selected subset of estimators for improved readability.

This code generates Figures 15 to 21 only, from which the main figures are derived.

---

## Data Files

This project relies entirely on synthetic data generated from pre-specified data-generating processes (DGPs).

- All DGPs are implemented in `code/simulate_data.R`, defining how covariates, treatment assignment, event times, and censoring are simulated.
