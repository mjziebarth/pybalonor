/*
 * Bayes of lognormal
 */

#include <cstddef>

#ifndef INDIAPALEALE_LOGNORMAL_HPP
#define INDIAPALEALE_LOGNORMAL_HPP

namespace indiapaleale {

class LogNormalPosterior {
public:
	static int log_posterior(const size_t M, const double* l0, const double* l1,
	                         double* log_post,
	                         const size_t N, const double* X,
	                         const double l0_min, const double l0_max,
	                         const double l1_min, const double l1_max);

	static int posterior(const size_t M, const double* l0, const double* l1,
	                     double* post,
	                     const size_t N, const double* X,
	                     const double l0_min, const double l0_max,
	                     const double l1_min, const double l1_max);

	static int log_mean_posterior(const size_t M, const double* mu,
	                              double* log_posterior,
	                              const size_t N, const double* X,
	                              const double l0_min, const double l0_max,
	                              const double l1_min, const double l1_max
	                             );

	static int log_posterior_predictive(const size_t M, const double* x,
	                                    double* log_post_pred,
	                                    const size_t N, const double* X,
	                                    const double l0_min,
	                                    const double l0_max,
	                                    const double l1_min,
	                                    const double l1_max);

};

}

#endif