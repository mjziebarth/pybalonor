# Code for Bayesian analysis of the lognormal distribution.

import numpy as np

cdef extern from "lognormal.hpp" namespace "indiapaleale" nogil:
    cdef cppclass LogNormalPosterior:
        @staticmethod
        int log_posterior(const size_t M, const double* l0,
                          const double* l1, double* log_post,
                          const size_t N, const double* X,
                          const double l0_min, const double l0_max,
                          const double l1_min, const double l1_max) except+

        @staticmethod
        int posterior(const size_t M, const double* l0, const double* l1,
                      double* post,
                      const size_t N, const double* X,
                      const double l0_min, const double l0_max,
                      const double l1_min, const double l1_max) except+

        @staticmethod
        int log_mean_posterior(const size_t M, const double* mu,
                               double* log_posterior,
                               const size_t N, const double* X,
                               const double l0_min, const double l0_max,
                               const double l1_min, const double l1_max
                              ) except+


def log_normal_log_posterior(const double[::1] l0, const double[::1] l1,
                             const double[::1] X, double l0_min, double l0_max,
                             double l1_min, double l1_max):
    """
    Compute the logarithm of the posterior.
    """
    if l0.shape[0] != l1.shape[0]:
        raise RuntimeError("Shapes of l0 and l1 must be equal!")
    cdef size_t M = l0.shape[0]
    cdef size_t N = X.shape[0]
    if M == 0:
        return np.empty(0)
    if N == 0:
        raise RuntimeError("No data given.")
    if l0_max <= l0_min:
        raise RuntimeError("l0_max has to be larger than l0_min.")
    if l1_max <= l1_min:
        raise RuntimeError("l1_max has to be larger than l1_min.")
    if l1_min < 0:
        raise RuntimeError("l0_min has to be non-negative.")

    cdef double[::1] log_post = np.empty(M)
    with nogil:
        LogNormalPosterior.log_posterior(M, &l0[0], &l1[0], &log_post[0],
                                         N, &X[0], l0_min, l0_max, l1_min,
                                         l1_max)

    return log_post.base


def log_normal_log_mean_posterior(const double[::1] mu, const double[::1] X,
                                  double l0_min, double l0_max, double l1_min,
                                  double l1_max):
    """
    Compute the logarithm of the posterior.
    """
    cdef size_t M = mu.shape[0]
    cdef size_t N = X.shape[0]
    if M == 0:
        return np.empty(0)
    if N == 0:
        raise RuntimeError("No data given.")
    if l0_max <= l0_min:
        raise RuntimeError("l0_max has to be larger than l0_min.")
    if l1_max <= l1_min:
        raise RuntimeError("l1_max has to be larger than l1_min.")
    if l1_min < 0:
        raise RuntimeError("l0_min has to be non-negative.")

    cdef double[::1] log_mean_post = np.empty(M)
    with nogil:
        LogNormalPosterior.log_mean_posterior(M, &mu[0], &log_mean_post[0],
                                              N, &X[0], l0_min, l0_max, l1_min,
                                              l1_max)

    return log_mean_post.base