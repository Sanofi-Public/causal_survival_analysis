# Causal Survival Analysis Project (R + renv)

This project uses [`renv`](https://rstudio.github.io/renv/) to manage R package dependencies in an isolated and reproducible environment.

## Quickstart

### 1. Open the project in R or RStudio

Make sure you're working **in the root directory** of the project, where the `renv.lock` and `renv/` folder are located.

You can check your current working directory with:

```r
getwd()
```

### 2. Automatic environment activation

This project includes a `.Rprofile` file containing the line to active renv.
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
