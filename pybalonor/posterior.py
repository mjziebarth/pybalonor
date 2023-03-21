# The log-normal posterior.

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
        Probability density.
        """
        ld = self.log_density(l0, l1)
        return np.exp(ld, out=ld)


    def log_density(self, l0: Vector, l1: Vector) -> np.ndarray:
        """
        Logarithm of the probability density.
        """
        l0,l1 = np.broadcast_arrays(l0, l1)
        shape = l0.shape
        l0 = np.array(l0, copy=False, order='C').reshape(-1)
        l1 = np.array(l1, copy=False, order='C').reshape(-1)
        lpdf = self.posterior.log_posterior(l0, l1)
        return lpdf.reshape(shape)


    def mean_pdf(self, mu: Vector) -> np.ndarray:
        """
        Density of the distribution mean.
        """
        lmd = self.log_mean_pdf(mu)
        return np.exp(lmd, out=lmd)


    def log_mean_pdf(self, mu: Vector) -> np.ndarray:
        """
        Logarithm of the density of the distribution mean.
        """
        mu = np.array(mu, copy=False, order='C')
        shape = mu.shape
        lmd = self.posterior.log_mean_posterior(mu.reshape(-1))
        return lmd


    def predictive_pdf(self, x: Vector) -> np.ndarray:
        """
        Predictive density.
        """
        lpdf = self.log_predictive_pdf(x)
        return np.exp(lpdf, out=lpdf)


    def log_predictive_pdf(self, x: Vector) -> np.ndarray:
        """
        Logarithm of the predictive density.
        """
        x = np.array(x, copy=False, order='C')
        shape = x.shape
        lpdf = self.posterior.log_posterior_predictive(x.reshape(-1))
        return lpdf.reshape(shape)


    def predictive_cdf(self, x: Vector) -> np.ndarray:
        """
        Predictive cumulative distribution function.
        """
        x = np.array(x, copy=False, order='C')
        shape = x.shape
        cdf = self.posterior.posterior_predictive_cdf(x.reshape(-1))
        return cdf.reshape(shape)


    def predictive_ccdf(self, x: Vector) -> np.ndarray:
        """
        Predictive complementary distribution function (or survivor
        function).
        """
        x = np.array(x, copy=False, order='C')
        shape = x.shape
        ccdf = self.posterior.posterior_predictive_ccdf(x.reshape(-1))
        return ccdf.reshape(shape)


    def predictive_quantiles(self, q: Vector) -> np.ndarray:
        """
        Quantiles of the predictive distribution.
        """
        q = np.array(q, copy=False, order='C')
        shape = q.shape
        x = self.posterior.posterior_predictive_quantiles(q.reshape(-1))
        return x.reshape(shape)


    def predictive_tail_quantiles(self, q: Vector) -> np.ndarray:
        """
        Tail quantiles of the predictive distribution (or quantiles
        of the complementary predictive distribution).
        """
        q = np.array(q, copy=False, order='C')
        shape = q.shape
        x = self.posterior.posterior_predictive_tail_quantiles(q.reshape(-1))
        return x.reshape(shape)
