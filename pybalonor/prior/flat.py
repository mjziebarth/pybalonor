# Flat prior.

from .prior import Prior

class FlatPrior(Prior):
    """
    Flat prior of the log-normal distribution.
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
