# SPDX-License-Identifier: EUPL-1.2
#
# The log-normal posterior.
#
# Authors: Malte J. Ziebarth (mjz.science@fmvkb.de)
#
# Copyright (C) 2023 Malte J. Ziebarth
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.

import numpy as np
from typing import Iterable, Union
from .prior import Prior, FlatPrior
from .bayes import CyLogNormalPosterior

# Types:
Vector = Union[Iterable[float],float]

class Posterior:
    """
    The posterior distribution of the log-normal distribution
    parameters given a sample and prior.

    Parameters
    ----------
    sample : array_like
       The sample for which to compute the posterior distribution.
    prior : FlatPrior
       The prior distribution.
    n_chebyshev : int, optional
       The number of Chebyshev points used in the barycentric Lagrange
       interpolation of the posterior predictive CDF for the posterior
       predictive quantiles. Defaults to 100.
    """
    def __init__(self, sample: Iterable[float], prior: Prior,
                 n_chebyshev: int = 100):
        if not isinstance(prior, Prior):
            raise TypeError("'prior' has to be an instance of a Prior "
                            "subclass.")
        # Sample as a flat numpy array:
        sample = np.array(sample, copy=False).flatten()

        # Number of Chebyshev points for the evaluation of the CDF (and
        # related functions):
        n_chebyshev = int(n_chebyshev)
        if n_chebyshev <= 2:
            raise ValueError("Need at least 2 Chebyshev points.")

        # Compute the posterior object:
        if isinstance(prior, FlatPrior):
            self.posterior = CyLogNormalPosterior(sample, prior.l0_min,
                                                  prior.l0_max, prior.l1_min,
                                                  prior.l1_max, n_chebyshev)
        else:
            raise NotImplementedError("Only 'FlatPrior' priors implemented.")


    def density(self, l0: Vector, l1: Vector) -> np.ndarray:
        """
        Posterior probability density.

        Parameters
        ----------
        l0 : array_like
           Vector of parameters :math:`l_0` at which to evaluate the posterior
           distribution.
        l1 : array_like
           Vector of parameters :math:`l_1` at which to evaluate the posterior
           distribution.

        Returns
        -------
        density : array_like
           NumPy array of the posterior density evaluated at :python:`l0`
           and :python:`l1`.
        """
        ld = self.log_density(l0, l1)
        return np.exp(ld, out=ld)


    def log_density(self, l0: Vector, l1: Vector) -> np.ndarray:
        """
        Logarithm of the posterior probability density.

        Parameters
        ----------
        l0 : array_like
           Vector of parameters :math:`l_0` at which to evaluate the posterior
           distribution.
        l1 : array_like
           Vector of parameters :math:`l_1` at which to evaluate the posterior
           distribution.

        Returns
        -------
        density : array_like
           NumPy array of the logarithm of the posterior density evaluated at
           :python:`l0` and :python:`l1`.
        """
        l0,l1 = np.broadcast_arrays(l0, l1)
        shape = l0.shape
        l0 = np.array(l0, copy=False, order='C').reshape(-1)
        l1 = np.array(l1, copy=False, order='C').reshape(-1)
        lpdf = self.posterior.log_posterior(l0, l1)
        return lpdf.reshape(shape)


    def mean_pdf(self, mu: Vector) -> np.ndarray:
        """
        Posterior density of the distribution mean :math:`\\mu`.

        Parameters
        ----------
        mu : array_like
           Vector of the distribution mean :math:`\\mu` at which to evaluate
           its posterior distribution.

        Returns
        -------
        pdf : array_like
           NumPy array of the posterior density evaluated at :python:`mu`.
        """
        lmd = self.log_mean_pdf(mu)
        return np.exp(lmd, out=lmd)


    def log_mean_pdf(self, mu: Vector) -> np.ndarray:
        """
        Logarithm of the posterior density of the distribution mean
        :math:`\\mu`.

        Parameters
        ----------
        mu : array_like
           Vector of the distribution mean :math:`\\mu` at which to evaluate
           its posterior distribution.

        Returns
        -------
        log_pdf : array_like
           NumPy array of the logarithm of the posterior density evaluated at
           :python:`mu`.
        """
        mu = np.array(mu, copy=False, order='C')
        shape = mu.shape
        lmd = self.posterior.log_mean_posterior(mu.reshape(-1))
        return lmd


    def predictive_pdf(self, x: Vector) -> np.ndarray:
        """
        Posterior predictive density.

        Parameters
        ----------
        x : array_like
           Vector of :math:`x` at which to evaluate the posterior
           predictive distribution.

        Returns
        -------
        pdf : array_like
           NumPy array of the posterior predictive density evaluated at
           :python:`x`.
        """
        lpdf = self.log_predictive_pdf(x)
        return np.exp(lpdf, out=lpdf)


    def log_predictive_pdf(self, x: Vector) -> np.ndarray:
        """
        Logarithm of the posterior predictive density.

        Parameters
        ----------
        x : array_like
           Vector of :math:`x` at which to evaluate the posterior
           predictive distribution.

        Returns
        -------
        pdf : array_like
           NumPy array of the logarithm of the posterior predictive density
           evaluated at :python:`x`.
        """
        x = np.array(x, copy=False, order='C')
        shape = x.shape
        lpdf = self.posterior.log_posterior_predictive(x.reshape(-1))
        return lpdf.reshape(shape)


    def predictive_cdf(self, x: Vector) -> np.ndarray:
        """
        Posterior predictive cumulative distribution function.

        Parameters
        ----------
        x : array_like
           Vector of :math:`x` at which to evaluate the posterior
           predictive distribution.

        Returns
        -------
        cdf : array_like
           NumPy array of the cumulative posterior predictive distribution
           evaluated at :python:`x`.
        """
        x = np.array(x, copy=False, order='C')
        shape = x.shape
        cdf = self.posterior.posterior_predictive_cdf(x.reshape(-1))
        return cdf.reshape(shape)


    def predictive_ccdf(self, x: Vector) -> np.ndarray:
        """
        Posterior predictive complementary distribution function (or survivor
        function).

        Parameters
        ----------
        x : array_like
           Vector of :math:`x` at which to evaluate the posterior
           predictive distribution.

        Returns
        -------
        ccdf : array_like
           NumPy array of the complementary cumulative posterior predictive
           distribution evaluated at :python:`x`.
        """
        x = np.array(x, copy=False, order='C')
        shape = x.shape
        ccdf = self.posterior.posterior_predictive_ccdf(x.reshape(-1))
        return ccdf.reshape(shape)


    def predictive_quantiles(self, q: Vector) -> np.ndarray:
        """
        Quantiles of the posterior predictive distribution.

        Parameters
        ----------
        q : array_like
           Vector of quantiles :math:`q` of the posterior predictive
           distribution to compute.

        Returns
        -------
        x : array_like
           NumPy array of the arguments :math:`x` of the posterior predictive
           distribution corresponding to the quantiles :math:`q`.
        """
        q = np.array(q, copy=False, order='C')
        shape = q.shape
        x = self.posterior.posterior_predictive_quantiles(q.reshape(-1))
        return x.reshape(shape)


    def predictive_tail_quantiles(self, q: Vector) -> np.ndarray:
        """
        Tail quantiles of the posterior predictive distribution (or quantiles
        of the complementary posterior predictive distribution).

        Parameters
        ----------
        q : array_like
           Vector of quantiles :math:`q` of the complementary cumulative
           posterior predictive distribution to compute.

        Returns
        -------
        x : array_like
           NumPy array of the arguments :math:`x` of the posterior predictive
           distribution corresponding to the quantiles :math:`1-q`.
        """
        q = np.array(q, copy=False, order='C')
        shape = q.shape
        x = self.posterior.posterior_predictive_tail_quantiles(q.reshape(-1))
        return x.reshape(shape)
