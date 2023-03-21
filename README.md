# pybalonor
This Python package helps to perform a Bayesian analysis of log-normally
distributed data
(*PYthon package for Bayesian Analysis of the LOg-NORmal distribution*).
Performing a Bayesian analysis of log-normally distributed data requires
care in the prior choice to yield posterior predictive distributions with
finite moments (e.g. Fabrizi & Trivisano, [2012](https://doi.org/10.1214/12-BA733)).

This package uses a simple uniform prior for the log-location and log-variance
parameter. The problem of normalizing the posterior of the mean is solved by
imposing a finite upper bound on the log-variance parameter.

## Installation and Requirements
The following software is required to install pybalonor:
- A modern C++ compiler
- Boost Math
- The Meson build system
- Cython
- NumPy
- [Mebuex](https://github.com/mjziebarth/Mebuex)

The Python package can be built from the repository's root directory using
the setuptools build system. For instance, you may call the following command
from the repository's root directory:
```bash
pip install --user .
```

## Usage
Currently, pybalonor provides one class, `CyLogNormalPosterior`:
```python
class CyLogNormalPosterior:
    def __init__(self, X, l0_min, l0_max, l1_min, l1_max):
        pass

    def log_posterior(self, l0, l1):
        pass

    def log_posterior_predictive(self, x):
        pass

    def posterior_predictive(self, x):
        pass

    def posterior_predictive_cdf(self, x):
        pass

    def log_mean_posterior(self, mu):
        pass
```
The parameters are as follows:
| Parameter | Type  | Purpose |
| --------- | ----- | ---------------------------------------------------------------------- |
| `X`       | dbuf1 | The data set.                                                          |
| `x`       | dbuf1 | Where to evaluate the posterior predictive (same dimension as `X`).    |
| `mu`      | dbuf1 | Log-Normal distribution mean (evaluated as density over the posterior) |
| `l0`      | dbuf1 | Log-location parameter $l_0$ at which to evaluate the posterior.       |
| `l1`      | dbuf1 | Log-variance parameter $l_1$ (like `l0`)                               |
| `l0_min`  | float | Minimum of log-location parameter for prior.                           |
| `l0_max`  | float | Maximum of log-location parameter.                                     |
| `l1_min`  | float | Minimum of log-variance parameter for prior.                           |
| `l1_max`  | float | Maximum of log-variance parameter.                                     |

Note: dbuf1 refers to a C-contiguous buffer of doubles (e.g. a one-dimensional NumPy array).
