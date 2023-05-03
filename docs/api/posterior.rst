=========
Posterior
=========

Given a :ref:`Prior <prior>` and a data set
:math:`\mathcal{N}=\{x_i\,:\,i=1..n\}`, pybalonor comes with the
:py:class:`~pybalonor.Posterior` class that allows to evaluate the posterior
predictive distribution, its quantiles, and the posterior distribution of the
distribution mean.

The formulae evaluated for the posterior ar given in
:ref:`Scientific Background <Scientific Background>`, and [Z2023]_ details how
these equations are evaluated.

|

:py:class:`~pybalonor.Posterior`
--------------------------------------

.. role:: python(code)
   :language: python

.. autoclass:: pybalonor.Posterior
   :members: