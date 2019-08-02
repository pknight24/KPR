# Kernel Penalized Regression

The `KPR` package provides estimation and inference methods for kernel
penalized regression models, designed for doubly-structured high
dimensional data.

## Installation

Install with `devtools::install_github("pknight24/KPR")`

## Usage

Model fitting is performed with the `KPR()` function.

```{r}
kpr.out <- KPR(designMatrix = X, Y = Y, H = H, Q = Q)
```

The `kpr.out` object has class `KPR`, includes all of the data used to
fit the model, as well as coefficient estimates.

Several high-dimensional inference procedures are provided for testing
the significance of regression coefficients. You can run tests on your
fitted model with the `inference()` function, specifying the procedure
to use with the `method` parameter.

```{r}
inference(kpr.out, method = "GMD")
```

For a more thorough example, see [the vignette](pknight24.github.io/KPR).
