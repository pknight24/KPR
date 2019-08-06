# Kernel Penalized Regression

The `KPR` package provides estimation and inference methods for kernel
penalized regression models, designed for doubly-structured high
dimensional data. An explanation of the theory and usage of kernel
penalized regression can be found [here](https://projecteuclid.org/euclid.aoas/1520564483).

## Installation

Install with `devtools::install_github("pknight24/KPR")`

## Usage

Model fitting is performed with the `KPR()` function.

```{r}
kpr.out <- KPR(designMatrix = X, Y = Y, H = H, Q = Q)
```

The `kpr.out` object has class `KPR` and includes all of the data used to
fit the model, as well as coefficient estimates and p-values based on the GMD inference.

For a more thorough example, see [the vignette](https://pknight24.github.io/KPR).
