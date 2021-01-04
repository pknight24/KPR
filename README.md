# Kernel-penalized regression

The `KPR` package provides estimation and inference methods for kernel-penalized regression models, designed for doubly-structured high
dimensional data. An explanation of the theory and usage of kernel-penalized regression can be found [here](https://projecteuclid.org/euclid.aoas/1520564483).

## Installation

Install with `devtools::install_github("pknight24/KPR")`

## Usage

Model fitting is performed with the `KPR()` function.

```{r}
kpr.out <- KPR(X = X, Y = Y, H = H, Q = Q)
```

The `kpr.out` object has class `KPR` and includes all of the data used to
fit the model, as well as coefficient estimates and tuning parameters.

The package also provides support for variable selection using the GMD inference.

```{r}
kpr.out <- inference(kpr.out)
```

This new `kpr.out` object also contains a list of p-values corresponding to each variable in the model.

For a more detailed example, see the [vignette](https://pknight24.github.io/KPR).
