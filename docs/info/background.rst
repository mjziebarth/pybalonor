.. _Scientific Background:

=====================
Scientific Background
=====================

Preface
-------
If you are looking for *an* analysis of the log-normal distribution, you might
likely want to check out the R package
`BayesLN <https://cran.r-project.org/web/packages/BayesLN/index.html>`_ by
Gardini, Fabrizi, and Trivisano. Their conjugate prior is more sophisticated
than the flat prior of pybalonor, and, from limited analysis, seems to lead to
tighter posterior bounds.

If instead you are looking for an analysis based on a flat prior, looking for a
Python solution, or working with a large data set, go ahead!

Details about the maths of pybalonor that go beyond this doc page are given in
[Z2023]_.


Bayesian Analysis of the Log-Normal Distribution with a Flat Prior
------------------------------------------------------------------
Performing a Bayesian analysis of log-normally distributed data requires
care in the prior choice to yield posterior predictive distributions with
finite moments (e.g. [F2012]_).

This package uses a simple uniform prior

.. math ::

   \phi(l_0, l_1 \,|\, l_0^\text{min}, l_0^\text{max}, l_1^\text{min},
       l_1^\text{max}) = \frac{1}{(l_0^\text{max} - l_0^\text{min})
                         (l_1^\text{max} - l_1^\text{min})}

for the log-location parameter
:math:`l_0` and log-variance parameter :math:`l_1` of the log-normal
distribution

.. math ::

   p\big(x\,|\, l_0, l_1 \big) = \frac{1}{\sqrt{2\pi} l_1 x }
   \exp\!\left(-\frac{\big(\ln x - l_0\big)^2}{2 l_1^2}\right)

within the bounds

.. math ::

   l_0^\text{min}>0\,, \quad l_0^\text{max}>l_0^\text{min},\quad
   l_1^\text{min}>0\,,\quad \text{and} \quad l_1^\text{max}>l_1^\text{min}\,.

The problem of normalizing the posterior of the mean is solved by
imposing a finite upper bound on the log-variance parameter. The prior is
represented by the :py:class:`~pybalonor.prior.FlatPrior` class.

The pybalonor package allows to compute a number of posterior distributions
given a sample

.. math ::

   \mathcal{N} = \{x_i \,|\, i = 1..n\}\,.

The computation of these posterior quantities are performed using the
:py:class:`~pybalonor.Posterior` class. Objects of this class are intialized
with a prior and a sample, and they represent the frozen posterior state given
the prior and evidence.


Posterior Distribution
""""""""""""""""""""""
The posterior distribution of the parameters :math:`l_0` and :math:`l_1` is

.. math ::

   p\big(l_0, l_1 \,|\,\mathcal{N}\big)
      \sim \mathcal{L}\big(l_0, l_1\,|\, \mathcal{N}\big)

where

.. math ::

   \mathcal{L}\big(l_0, l_1 \,|\, \mathcal{N}\big)
       = \phi(l_0, l_1)\prod\limits_{i=1}^n \frac{1}{\sqrt{2\pi l_1 x_i}}
         \exp\!\left(-\frac{\big(\ln x_i - l_0\big)^2}{2 l_1^2} \right)

is the likelihood of :math:`l_0` and :math:`l_1` given the evidence
:math:`\mathcal{N}` and prior :math:`\phi`. Details about the computation of the
normalization constant for :math:`p\big(l_0, l_1 \,|\,\mathcal{N}\big)` are
given in [Z2023]_.


Posterior Predictive Distribution
"""""""""""""""""""""""""""""""""
The posterior predictive distribution :math:`p\big(x\,|\,\mathcal{N}\big)`
gives the posterior probability density of the random variable :math:`X` given
a sample :math:`\mathcal{N}`. It is computed using the posterior density of
:math:`l_0` and :math:`l_1`:

.. math ::

   p\big(x\,|\,\mathcal{N}\big)
      = \int\limits_{l_0^\text{min}}^{l_0^\text{max}}\!\!\!\mathrm{d}l_0\!\!
        \int\limits_{l_1^\text{min}}^{l_1^\text{max}}\!\!\!\mathrm{d}l_1\;
        p\big(x\,|\,l_0, l_1\big) p\big(l_0, l_1\,|\,\mathcal{N}\big)\,.

Details about the computation are given in [Z2023]_.

Posterior Distribution of the Mean
""""""""""""""""""""""""""""""""""
The mean of the log-normal distribution :math:`p\big(x\,|\,l_0, l_1\big)`
parameterized by :math:`l_0` and :math:`l_1` is

.. math ::

   \mu = \exp\!\left(l_0 + \frac{1}{2}l_1^2\right)\,.

The posterior distribution of :math:`\mu` given a sample :math:`\mathcal{N}`
follows by accumulating all posterior mass
:math:`p\big(l_0,l_1\,|\,\mathcal{N}\big)` for those parameters :math:`l_0`
and :math:`l_1` that amount to :math:`\mu`. This amounts to computing
:math:`p(\mu)` so that

.. math ::

   p(\mu\,|\,\mathcal{N})\mathrm{d}\mu
      = \!\!\iint\limits_{A(\mu,\mathrm{d}\mu)} \!\!\mathrm{d}l_0\,
        \mathrm{d}l_1\; p(l_0, l_1\,|\,\mathcal{N})\,,

where :math:`A(\mu,\mathrm{d}\mu)` is the subset of the :math:`(l_0,l_1)`
parameter space in which the distribution mean is in the interval
:math:`[\mu,\mu+\mathrm{d}\mu]`. Details about the computation are given in
[Z2023]_.




Bibliography
------------

.. [F2012]  Fabrizi, E. & Trivisano, C. (2012). Bayesian Estimation of
            Log-Normal Means with Finite Quadratic Expected Loss. Bayesian
            Analysis 7(4), 975-996.
            https://doi.org/10.1214/12-BA733

.. [Z2023] Ziebarth, Malte J. (2023). Bayesian Analysis of the Log-Normal
           Distribution with a Flat Prior: The pybalonor package.