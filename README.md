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

## Preface
If you are looking for *an* analysis of the log-normal distribution, you might
likely want to check out the R package
[BayesLN](https://cran.r-project.org/web/packages/BayesLN/index.html) by
Gardini, Fabrizi, and Trivisano. Their conjugate prior is more sophisticated
than the flat prior of pybalonor, and, from limited analysis, seems to lead to
tighter posterior bounds.

If instead you are looking for an analysis based on a flat prior, looking for a
Python solution, or working with a large data set, go ahead!


## Installation and Requirements
The following software is required to install pybalonor:
- A modern C++ compiler
- Boost Math (v1.80.0 or later recommended for numerical stability)
- The Meson build system
- Cython
- NumPy
- Mebuex

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

For more information, visit the [pybalonor documentation](https://mjziebarth.github.io/pybalonor/).

## License
This software is licensed under the European Public License (EUPL) version 1.2
or later (`EUPL-1.2`). See the LICENSE file in this directory.

## Changelog
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### [1.0.0] - 2023-05-04
#### Added
- Initial release.
