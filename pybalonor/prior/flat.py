# SPDX-License-Identifier: EUPL-1.2
#
# Flat prior.
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
