/*
 * Bayes of lognormal
 */

#include <cstddef>
#include <vector>
#include <limits>

#ifndef INDIAPALEALE_LOGNORMAL_HPP
#define INDIAPALEALE_LOGNORMAL_HPP

namespace indiapaleale {

class LogNormalPosterior {
public:
	LogNormalPosterior(const size_t N, const double* X, const double l0_min,
	                   const double l0_max, const double l1_min,
	                   const double l1_max);

	void log_posterior(const size_t M, const double* l0, const double* l1,
	                   double* log_post) const;

	void posterior(const size_t M, const double* l0, const double* l1,
	               double* post) const;

	void log_mean_posterior(const size_t M, const double* mu,
	                        double* log_posterior) const;

	void log_posterior_predictive(const size_t M, const double* x,
	                              double* log_post_pred) const;

	void posterior_predictive_cdf(const size_t M, const double* x,
	                              double* log_post_pred) const;

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

	struct post_pred_t {
		/* The following vector contains the joint set of data x_i and
		 * the evaluation point x (which can always be set to the */
		std::vector<long double> lX_xi;
		const long double dlga;

		post_pred_t(size_t N, long double dlga);
	};

	post_pred_t predictive_parameters() const;

	static prior_t sanity_check(double l0_min, double l0_max, double l1_min,
	                            double l1_max);
};

}

#endif