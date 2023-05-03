=====
Prior
=====

The pybalonor package supports a prior that is uniform in the mode of logs
:math:`l_0` and the standard deviation of logs :math:`l_1`,

.. math ::

   \phi(l_0, l_1 \,|\, l_0^\text{min}, l_0^\text{max}, l_1^\text{min},
        l_1^\text{max}) = \frac{1}{(l_0^\text{max} - l_0^\text{min})
                          (l_1^\text{max} - l_1^\text{min})}\,,

where the rectangular support of :math:`l_0` and :math:`l_1` is defined by

.. math ::

   l_0^\text{min}>0\,, \quad l_0^\text{max}>l_0^\text{min},\quad
   l_1^\text{min}>0\,,\quad \text{and} \quad l_1^\text{max}>l_1^\text{min}\,.

The uniform prior weights all parameter combinations :math:`(l_0,l_1)` within
the bounds

.. math ::

   l_0^\text{min} \leq l_0 \leq l_0^\text{max}

   l_1^\text{min} \leq l_1 \leq l_1^\text{max}

equally.

|

:py:class:`~pybalonor.prior.FlatPrior`
--------------------------------------

.. role:: python(code)
   :language: python

.. autoclass:: pybalonor.prior.FlatPrior
   :members: