# Code for Bayesian analysis of the lognormal distribution.

import numpy as np
from libcpp.memory cimport shared_ptr, make_shared
from cython.operator cimport dereference as deref

cdef extern from "lognormal.hpp" namespace "indiapaleale" nogil:
    cdef cppclass LogNormalPosterior:
        LogNormalPosterior(const size_t N, const double* X,
                           const double l0_min, const double l0_max,
                           const double l1_min, const double l1_max,
                           size_t n_chebyshev) except+

        @staticmethod
        int log_posterior(const size_t M, const double* l0,
                          const double* l1, double* log_post) except+

        @staticmethod
        int posterior(const size_t M, const double* l0, const double* l1,
                      double* post) except+

        @staticmethod
        int log_mean_posterior(const size_t M, const double* mu,
                               double* log_posterior) except+

        @staticmethod
        int log_posterior_predictive(const size_t M, const double* x,
                                     double* log_post_pred) except+

        @staticmethod
        int posterior_predictive_cdf(const size_t M, const double* x,
                                     double* log_post_pred) except+


cdef class CyLogNormalPosterior:
    """

    """
    cdef shared_ptr[LogNormalPosterior] _posterior

    def __init__(self, const double[::1] X, double l0_min, double l0_max,
                 double l1_min, double l1_max, size_t n_chebyshev):
        """

        """
        cdef size_t N = X.shape[0]
        if N == 0:
            raise RuntimeError("No data given.")
        if l0_max <= l0_min:
            raise RuntimeError("l0_max has to be larger than l0_min.")
        if l1_max <= l1_min:
            raise RuntimeError("l1_max has to be larger than l1_min.")
        if l1_min < 0:
            raise RuntimeError("l0_min has to be non-negative.")
        self._posterior = make_shared[LogNormalPosterior](N, &X[0], l0_min,
                                                          l0_max, l1_min,
                                                          l1_max, n_chebyshev)

    def log_posterior(self, const double[::1] l0, const double[::1] l1):
        """
        Compute the logarithm of the posterior.
        """
        if l0.shape[0] != l1.shape[0]:
            raise RuntimeError("Shapes of l0 and l1 must be equal!")
        cdef size_t M = l0.shape[0]
        if M == 0:
            return np.empty(0)

        if not self._posterior:
            raise RuntimeError("CyLogNormalPosterior not properly initialized.")

        cdef double[::1] log_post = np.empty(M)
        with nogil:
            deref(self._posterior).log_posterior(M, &l0[0], &l1[0],
                                                 &log_post[0])

        return log_post.base

    def log_mean_posterior(self, const double[::1] mu):
        """
        Compute the logarithm of the posterior.
        """
        cdef size_t M = mu.shape[0]
        if M == 0:
            return np.empty(0)

        if not self._posterior:
            raise RuntimeError("CyLogNormalPosterior not properly initialized.")

        cdef double[::1] log_mean_post = np.empty(M)
        with nogil:
            deref(self._posterior).log_mean_posterior(M, &mu[0],
                                                      &log_mean_post[0])

        return log_mean_post.base

    def log_posterior_predictive(self, const double[::1] x):
        """
        Compute the logarithm of the posterior predictive.
        """
        cdef size_t M = x.shape[0]
        if M == 0:
            return np.empty(0)

        if not self._posterior:
            raise RuntimeError("CyLogNormalPosterior not properly initialized.")

        cdef double[::1] log_post_pred = np.empty(M)
        with nogil:
            deref(self._posterior).log_posterior_predictive(M, &x[0],
                                                            &log_post_pred[0])

        return log_post_pred.base

    def posterior_predictive(self, const double[::1] x):
        """
        Compute the posterior predictive.
        """
        log_post_pred = self.log_posterior_predictive(x)
        np.exp(log_post_pred, out=log_post_pred)
        return log_post_pred


    def posterior_predictive_cdf(self, const double[::1] x):
        """
        Compute the posterior predictive CDF.
        """
        cdef size_t M = x.shape[0]
        if M == 0:
            return np.empty(0)

        if not self._posterior:
            raise RuntimeError("CyLogNormalPosterior not properly initialized.")

        cdef double[::1] post_pred_cdf = np.empty(M)
        with nogil:
            deref(self._posterior).posterior_predictive_cdf(M, &x[0],
                                                            &post_pred_cdf[0])

        return post_pred_cdf.base
