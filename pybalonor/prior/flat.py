# Flat prior.

from .prior import Prior

class FlatPrior(Prior):
    """
    Flat prior of the log-normal distribution.

    Parameters
    ----------
    l0_min : float
        Minimum :math:`l_0` within the prior support.
    l0_max : float
        Maximum :math:`l_0` within the prior support. Has to fulfill
        :math:`\,l_0^\\text{max} > l_0^\\text{min}`.
    l1_min : float
        Minimum :math:`l_1` within the prior support. Has to fulfill
        :math:`\,l_1^\\text{min} \geq 0`.
    l1_max : float
        Maximum :math:`l_1` within the prior support. Has to fulfill
        :math:`\,l_1^\\text{max} > l_1^\\text{min}`.
    """
    l0_min: float
    l0_max: float
    l1_min: float
    l1_max: float

    def __init__(self, l0_min: float, l0_max: float, l1_min: float,
                 l1_max: float):
        self.l0_min = float(l0_min)
        self.l0_max = float(l0_max)
        self.l1_min = float(l1_min)
        self.l1_max = float(l1_max)
