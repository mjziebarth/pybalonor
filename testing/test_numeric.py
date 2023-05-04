# SPDX-License-Identifier: EUPL-1.2
#
# Test some numerical properties of the posterior distributions using
# pre-evaluated test cases.
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
from pybalonor import Posterior, FlatPrior

def test_mean_posterior():
    """
    Ensures that the posterior of the mean is correctly normalized
    for a test data set.
    """
    # Normally-distributed random values with mean 0 and variance 1:
    X = np.array([-1.7992800466709749,  1.2487606638593516, -0.084337761710497,
                   0.6510473809233797,  0.8433388434280636,  0.6817916295962333,
                   0.5758152867770215, -1.135789421044811,  -0.3399192133810009,
                  -0.0162180438823278, -1.218121832512448,  -0.1862849433706261,
                   1.0157978378575432])

    # Transform to log-normal:
    l0 = 1.0
    l1 = 0.5
    X = np.exp(l1 * X + l0)

    # Prior parameters the resulting limits on µ:
    L0_MIN = -10.0
    L0_MAX = 2.0
    L1_MIN = 0.0
    L1_MAX = 3.0
    prior = FlatPrior(L0_MIN, L0_MAX, L1_MIN, L1_MAX)
    mu_min = np.exp(L0_MIN + 0.5*L1_MIN**2)
    mu_max = np.exp(L0_MAX + 0.5*L1_MAX**2)

    # Evaluate the posterior within these limits:
    mu = np.linspace(mu_min, mu_max, 10000)
    posterior = Posterior(X, prior)
    y = posterior.mean_pdf(mu)

    # Ensure that the posterior is normed:
    assert abs(0.5 * (np.sum(y[1:]) + np.sum(y[:-1])) * (mu[1] - mu[0])
               - 1.0) < 1e-4

    # Ensure that the mode is close to the actual µ:
    assert abs(mu[np.argmax(y)] - np.exp(l0 + 0.5*l1**2)) < 0.1
