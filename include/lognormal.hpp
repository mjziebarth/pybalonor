/*
 * SPDX-License-Identifier: EUPL-1.2
 *
 * Bayesian inference of the log-normal distribution.
 *
 * Authors: Malte J. Ziebarth (mjz.science@fmvkb.de)
 *
 * Copyright (C) 2023 Malte J. Ziebarth
 *
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
 * the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and
 * limitations under the Licence.
 */

#include <cstddef>
#include <vector>
#include <limits>

#include <cdfbli.hpp>

#ifndef INDIAPALEALE_LOGNORMAL_HPP
#define INDIAPALEALE_LOGNORMAL_HPP

namespace indiapaleale {

class LogNormalPosterior {
public:
	LogNormalPosterior(const size_t N, const double* X, const double l0_min,
	                   const double l0_max, const double l1_min,
	                   const double l1_max, size_t n_chebyshev);

	void log_posterior(const size_t M, const double* l0, const double* l1,
	                   double* log_post) const;

	void posterior(const size_t M, const double* l0, const double* l1,
	               double* post) const;

	void log_mean_posterior(const size_t M, const double* mu,
	                        double* log_posterior) const;

	void log_posterior_predictive(const size_t M, const double* x,
	                              double* log_post_pred) const;

	void posterior_predictive_cdf(const size_t M, const double* x,
	                              double* pred_cdf) const;

	void posterior_predictive_ccdf(const size_t M, const double* x,
	                               double* pred_ccdf) const;

	void posterior_predictive_quantiles(const size_t M, const double* q,
	                                    double* pred_quantiles) const;

	void posterior_predictive_tail_quantiles(const size_t M, const double* q,
	                                         double* pred_tail_quantiles) const;

private:
	struct prior_t {
		double l0_min;
		double l0_max;
		double l1_min;
		double l1_max;
	};

	const prior_t prior;
	const std::vector<long double> lX;
	const long double x_sum;
	const long double lx_sum;
	const long double lga;

	mutable long double lI = -std::numeric_limits<long double>::infinity();
	void compute_lI() const;

	/*
	 * Variables for the log_mean_posterior:
	 */
	mutable long double lmp_log_norm = -1.0;

	/*
	 * Variables for the posterior predictive:
	 */

	struct post_pred_t {
		/* The following vector contains the joint set of data x_i and
		 * the evaluation point x (which can always be set to the */
		std::vector<long double> lX_xi;
		const long double dlga;

		post_pred_t(size_t N, long double dlga);
	};

	post_pred_t predictive_parameters() const;

	/*
	 * Variables for the posterior predictive CDF and cCDF:
	 */
	typedef CenteredExpTanTransform<long double> pred_trans_t;
	typedef BarycentricLagrangeCDFInterpolator<pred_trans_t> pred_cdf_interp_t;
	typedef BarycentricLagrangeCDFInterpolator<pred_trans_t,true>
	        pred_ccdf_interp_t;
	size_t n_chebyshev;
	mutable std::optional<pred_cdf_interp_t> predictive_cdf_interpolator;
	mutable std::optional<pred_ccdf_interp_t> predictive_ccdf_interpolator;

	enum cumulative_t {
		CDF,
		CCDF
	};
	void init_predictive_cumulative_interpolator(cumulative_t which) const;


	static prior_t sanity_check(double l0_min, double l0_max, double l1_min,
	                            double l1_max);
};

}

#endif