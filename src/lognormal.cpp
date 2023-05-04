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

#include <lognormal.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <iostream>

#define BOOST_ENABLE_ASSERT_HANDLER
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>

/*
 * Handling asserts:
 */
void boost::assertion_failed(char const* expr, char const* function,
                             char const* file, long line)
{
	std::string msg(expr);
	msg.append("\n   in function ");
	msg.append(function);
	msg.append("\n   in file ");
	msg.append(file);
	msg.append("\n   in line");
	msg.append(std::to_string(line));
	throw std::runtime_error(msg);
}



using boost::math::quadrature::tanh_sinh;
using boost::math::quadrature::exp_sinh;
using boost::math::quadrature::gauss_kronrod;
using boost::math::tools::toms748_solve;
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::brent_find_minima;
using boost::math::tools::eps_tolerance;
using indiapaleale::LogNormalPosterior;
using indiapaleale::math;

template<typename real>
struct lognormal_params_t {
	real l0;
	real l1;
};

/*
 * log-normal maximum likelihood.
 */
template<typename real>
lognormal_params_t<real> log_normal_mle(const std::vector<long double>& lX)
{
	lognormal_params_t<real> p({0.0, 0.0});
	for (double lxi : lX){
		p.l0 += lxi;
		p.l1 += lxi * lxi;
	}
	p.l0 /= lX.size();
	p.l1 /= lX.size();
	p.l1 -= p.l0 * p.l0;

	return p;
}




/*
 * A function that computes the logarithm of the difference between two
 * regularized incomplete gamma functions with equal `a`.
 */
template<typename real>
real log_delta_gamma_1(real alpha, real x0, real x1)
{
	constexpr real tol = 1e-10 * std::numeric_limits<real>::epsilon();
	/* This function is made for large x0 & x1 where the difference
	 * between the two incomplete gamma functions is small.
	 */
	const real lx1 = math<real>::log(x1);
	real uk_x1 = 1.0;
	real uk_x0 = math<real>::exp((alpha - 1.0) * (math<real>::log(x0) - lx1)
	                             - (x0 - x1));
	real duk_x10 = uk_x1 - uk_x0;
	real S = duk_x10;
	size_t k = 1;
	while (duk_x10 > tol * S){
		uk_x1 *= (alpha - k) / x1;
		uk_x0 *= (alpha - k) / x0;
		duk_x10 = uk_x1 - uk_x0;
		S += duk_x10;
		++k;
	}
	return (alpha - 1.0) * lx1 - x1 + math<real>::log(S);
}

template<typename real>
real log_delta_gamma_div_loggamma(real alpha, real x0, real x1)
{
	/*
	 * This function computes the logarithm of the difference
	 * between two incomplete gamma functions with equal alpha:
	 *    log(gamma_p(alpha, x0) - gamma_p(alpha, x1)) - loggamma(a)
	 * which requires x0 > x1.
	 */
	constexpr real tol
	    = std::sqrt(std::numeric_limits<real>::epsilon());
	const real g0 = boost::math::gamma_p(alpha, x0);
	const real g1 = boost::math::gamma_p(alpha, x1);
	if (std::abs<real>(g0 - g1) > tol * g0)
		return math<real>::log(g0 - g1);
	return log_delta_gamma_1<real>(alpha, x0, x1) - math<real>::lgamma(alpha);
}


template<typename real>
real compute_log_integral_I(real a, real b, real c, real d,
                            const std::vector<long double>& lx)
{
	real alpha = 0.5 * lx.size() - 0.5;

	/* This scale will be used to keep the integrand in a numerically
	 * feasible range:
	 */
	real log_scale = 0.0;
	real l0_mle = 0.0;
	real I = 0.0;
	if (c == 0.0 && std::isinf(d)){
		/*
		 * Infinite interval with c == 0.
		 */
		{
			auto log_integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1);
				return log_integrand;
			};

			lognormal_params_t<real> p(log_normal_mle<real>(lx));
			l0_mle = p.l0;
			if (l0_mle < a || l0_mle > b)
				log_scale = std::max<real>(log_integrand(a), log_integrand(b));
			else
				log_scale = log_integrand(l0_mle);
		}

		if (std::isinf(log_scale))
			I = 0.0;
		else {
			auto integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1);
				return math<real>::exp(log_integrand - log_scale);
			};

			tanh_sinh<real> integrator;
			if (a < l0_mle && b > l0_mle)
				I =   integrator.integrate(integrand, a, l0_mle)
				    + integrator.integrate(integrand, l0_mle, b);
			else {
				I = integrator.integrate(integrand, a, b);
			}
		}
	} else if (c == 0.0){
		/*
		 * Finite interval with c == 0.
		 */
		real d2 = d * d;
		{
			auto log_integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1)
				    + math<real>::log(boost::math::gamma_q(alpha, C1/d2));
				return log_integrand;
			};

			lognormal_params_t<real> p(log_normal_mle<real>(lx));
			l0_mle = p.l0;
			if (l0_mle < a || l0_mle > b)
				log_scale = std::max<real>(log_integrand(a), log_integrand(b));
			else
				log_scale = log_integrand(l0_mle);
		}

		if (std::isinf(log_scale))
			I = 0.0;
		else {
			auto integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1)
				    + math<real>::log(boost::math::gamma_q(alpha, C1/d2));
				real I = math<real>::exp(log_integrand - log_scale);
				if (std::isinf(I) || std::isnan(I)){
					std::cout << "Integrand infinite!\nC1 = "
					          << C1
					          << "\nlog_integrand = " << log_integrand << "\nlog_scale = "
					          << log_scale << "\n" << std::flush;
				}
				return I;
			};

			tanh_sinh<real> integrator;
			if (a < l0_mle && b > l0_mle)
				I =   integrator.integrate(integrand, a, l0_mle)
				    + integrator.integrate(integrand, l0_mle, b);
			else {
				I = integrator.integrate(integrand, a, b);
			}
		}

	} else if (std::isinf(d)) {
		/*
		 * Infinite interval with c > 0.
		 */
		real c2 = c * c;
		{
			auto log_integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1)
				    + math<real>::log(boost::math::gamma_p(alpha, C1/c2));
				return log_integrand;
			};

			lognormal_params_t<real> p(log_normal_mle<real>(lx));
			l0_mle = p.l0;
			if (l0_mle < a || l0_mle > b)
				log_scale = std::max<real>(log_integrand(a), log_integrand(b));
			else
				log_scale = log_integrand(l0_mle);
		}

		if (std::isinf(log_scale))
			I = 0.0;
		else {
			auto integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1)
				    + math<real>::log(boost::math::gamma_p(alpha, C1/c2));
				return math<real>::exp(log_integrand - log_scale);
			};

			tanh_sinh<real> integrator;
			if (a < l0_mle && b > l0_mle)
				I =   integrator.integrate(integrand, a, l0_mle)
				    + integrator.integrate(integrand, l0_mle, b);
			else {
				I = integrator.integrate(integrand, a, b);
			}
		}
	} else {
		/*
		 * Finite interval with c > 0.
		 */
		real c2 = c * c;
		real d2 = d * d;
		{
			auto log_integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1)
					+ log_delta_gamma_div_loggamma(alpha, C1/c2, C1/d2);
				return log_integrand;
			};

			lognormal_params_t<real> p(log_normal_mle<real>(lx));
			l0_mle = p.l0;
			if (l0_mle < a || l0_mle > b)
				log_scale = std::max<real>(log_integrand(a), log_integrand(b));
			else
				log_scale = log_integrand(l0_mle);
		}

		if (std::isinf(I))
			I = 0.0;
		else {
			auto integrand = [&](real l0) -> real {
				real C1 = 0.0;
				for (double lxi : lx){
					const real dlx = lxi - l0;
					C1 += dlx * dlx;
				}
				C1 *= 0.5;
				real log_integrand = -alpha * math<real>::log(C1)
				    + log_delta_gamma_div_loggamma(alpha, C1/c2, C1/d2);
				real I = math<real>::exp(log_integrand - log_scale);
				return I;
			};

			tanh_sinh<real> integrator;
			if (a < l0_mle && b > l0_mle)
				I =   integrator.integrate(integrand, a, l0_mle)
				    + integrator.integrate(integrand, l0_mle, b);
			else {
				I = integrator.integrate(integrand, a, b);
			}
		}
	}

	return math<real>::log(I) + log_scale;
}

LogNormalPosterior::prior_t
LogNormalPosterior::sanity_check(double l0_min, double l0_max, double l1_min,
                                 const double l1_max)
{
	/* Sanity check: */
	if (l0_min > l0_max){
		throw std::runtime_error("l0_min > l0_max not allowed.");
	}
	if (l1_min > l1_max){
		throw std::runtime_error("l1_min > l1_max not allowed.");
	} else if (l1_min < 0){
		throw std::runtime_error("l1_min < 0 not allowed.");
	}
	return prior_t({.l0_min=l0_min, .l0_max=l0_max, .l1_min=l1_min,
	                .l1_max=l1_max});
}

static std::vector<long double> compute_lX(const size_t N, const double* X)
{
	if (N <= 1)
		throw std::runtime_error("Only samples with more than one data point "
		                         "supported.");
	std::vector<long double> lX(N);
	for (size_t i=0; i<N; ++i){
		const double xi = X[i];
		if (xi <= 0.0)
			throw std::runtime_error("Encountered x < 0.");
		lX[i] = math<long double>::log(xi);
	}
	/* Sort lX for later median computation: */
	std::sort(lX.begin(), lX.end());
	return lX;
}

LogNormalPosterior::LogNormalPosterior(const size_t N, const double* X,
                       const double l0_min, const double l0_max,
                       const double l1_min, const double l1_max,
                       size_t n_chebyshev)
   : prior(sanity_check(l0_min, l0_max, l1_min, l1_max)),
     lX(compute_lX(N,X)),
     x_sum(std::accumulate(X, X+N, 0.0)),
     lx_sum(std::accumulate(lX.cbegin(), lX.cend(), 0.0)),
     lga(math<long double>::lgamma(0.5 * N - 0.5)),
     n_chebyshev(n_chebyshev)
{
}

void LogNormalPosterior::compute_lI() const
{
	if (std::isinf(lI) && lI < 0.0)
		lI = compute_log_integral_I<long double>(prior.l0_min, prior.l0_max,
		                                         prior.l1_min, prior.l1_max,
		                                         lX);
}

void LogNormalPosterior::log_posterior(const size_t M, const double* l0,
                                       const double* l1, double* log_post) const
{
	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 not implemented in log_posterior.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 not implemented in log_posterior.");
	}

	/* Normalization constant: */
	constexpr long double ln2 = std::log((long double)2.0);
	compute_lI();

	/* Evaluate posterior: */
	//#pragma omp parallel for
	const size_t N = lX.size();
	for (size_t i=0; i<M; ++i){
		const double l0i = l0[i];
		const double l1i = l1[i];
		if (l0i < prior.l0_min || l0i > prior.l0_max || l1i < prior.l1_min
		    || l1i > prior.l1_max || l1i == 0.0)
			log_post[i] = -std::numeric_limits<double>::infinity();
		else {
			long double C1 = 0.0;
			for (double lxi : lX){
				const long double dlx = lxi - l0i;
				C1 += dlx * dlx;
			}
			C1 *= 0.5;
			log_post[i] = ln2 - lI - lga - N * math<long double>::log(l1i)
			              - C1 / (l1i * l1i);
		}
	}
}

void LogNormalPosterior::posterior(const size_t M, const double* l0,
                                   const double* l1, double* post) const
{
	log_posterior(M, l0, l1, post);
	for (size_t i=0; i<M; ++i)
		post[i] = std::exp(post[i]);
}

/*
 * This function evaluates the posterior of the mean for one particular
 * value of the mean.
 */
template<typename real>
real log_mean_posterior_eval(const real mu_i,
                             const std::vector<long double>& lX,
                             const real lx_sum,
                             const real l0_min, const real l0_max,
                             const real l1_min, const real l1_max,
                             const real log_norm
                          )
{
	typedef math<real> mth;
	const size_t N = lX.size();

	const real ln_mu = mth::log(mu_i);
	if (ln_mu <= l0_min){
		return -std::numeric_limits<real>::infinity();
	}

	const real lambda0
	   = std::min(mth::sqrt(std::max<real>(2.0 * (ln_mu - l0_max), 0.0)),
	              l1_min);
	const real lambda1
	   = std::min(mth::sqrt(std::max<real>(2.0 * (ln_mu - l0_min), 0.0)),
	              l1_max);
	if (lambda1 <= lambda0){
		return -std::numeric_limits<real>::infinity();
	}

	/*
	 * Compute the maximum of the integrand:
	 */
	real log_scale = 0.0;
	real l1_peak = 0.0;
	bool max_at_boundary = false;
	{
		auto log_integrand = [&](real l1) -> real {
			if (l1 == 0.0)
				return -std::numeric_limits<real>::infinity();
			const real dlx = ln_mu - 0.5 * l1 * l1;
			real S = 0.0;
			for (double lxi : lX){
				S += (lxi - dlx) * (lxi - dlx);
			}
			return -(N * mth::log(l1)) - 0.5 * S / (l1 * l1);
		};
		auto fun0 = [&](real l1) -> real {
			const real dlx = (ln_mu - 0.5 * l1 * l1);
			real S = 0.0;
			for (double lxi : lX){
				S += (lxi - dlx) * (lxi - dlx);
			}
			return N * (ln_mu - 1.0) - 0.5 * N * l1 * l1 - lx_sum
			       + S / (l1 * l1);
		};
		if ((fun0(lambda0) > 0) == (fun0(lambda1) > 0)){
			/* No change in sign -> maximum at the boundary */
			const real li_l = log_integrand(lambda0);
			const real li_r = log_integrand(lambda1);
			max_at_boundary = true;
			log_scale = std::max(li_l, li_r);
		} else {
			/*
			 * Use Newton-Raphson to compute the maximum of the
			 * log integrand.
			 */
			auto fun1 = [&](real l1) -> real {
				const real dlx = (ln_mu - 0.5 * l1 * l1);
				real S0 = 0.0;
				real S1 = 0.0;
				for (double lxi : lX){
					real di = (lxi - dlx) / l1;
					S0 += di * di;
					S1 += di;
				}
				return -(N * l1) - 2 * S0 / l1  + 2 * S1;
			};
			/* Initial guess: */
			if (2.0 * ln_mu - 2.0/N * lx_sum < 0.0)
				l1_peak = 0.5 * (lambda1 + lambda0);
			else
				l1_peak
				   = std::max(
				        std::min(
				           mth::sqrt(
				             std::max<real>(2 * ln_mu - 2.0/N * lx_sum,
				                            0.0)),
				           lambda1),
				        lambda0);
			const real tol
			  = mth::sqrt(std::numeric_limits<real>::epsilon());
			for (size_t j=0; j<100; ++j){
				real dl1 = - fun0(l1_peak) / fun1(l1_peak);
				dl1 = std::min(std::max(dl1, -0.9 * (l1_peak - lambda0)),
				               0.9 * (lambda1 - l1_peak));
				l1_peak += dl1;
				if (std::abs(dl1) < tol * l1_peak)
					break;
			}
			log_scale = log_integrand(l1_peak);
		}
	}
	/* Now integrate: */
	auto integrand = [&](real l1) -> real {
		if (l1 == 0.0)
			return 0.0;
		const real dlx = ln_mu - 0.5 * l1 * l1;
		real S = 0.0;
		for (double lxi : lX){
			S += (lxi - dlx) * (lxi - dlx);
		}
		real log = -(N * mth::log(l1)) - 0.5 * S / (l1 * l1);
		return mth::exp(log - log_scale);
	};
	tanh_sinh<real> integrator;
	real I_l1 = 0.0;
	/* On numerical difficulties, return NaN: */
	try {
		if (max_at_boundary){
			I_l1 = integrator.integrate(integrand, lambda0, lambda1);
		} else {
			I_l1 =   integrator.integrate(integrand, lambda0, l1_peak)
				   + integrator.integrate(integrand, l1_peak, lambda1);
		}
		return -log_norm + mth::log(I_l1) + log_scale - ln_mu;
	} catch (...) {
		std::cout << "quiet_NaN\n";
		return std::numeric_limits<double>::quiet_NaN();
	}
}

template<typename real>
std::pair<real,real>
maximum_posterior_mu(const std::vector<long double>& lX, const real lx_sum,
                     const real x_sum, const real l0_min, const real l0_max,
                     const real l1_min, const real l1_max)
{
	/*
	 * This function aims to find the maximum of the posterior in mu.
	 * It assumes that there is exactly one maximum.
	 */
	const size_t N = lX.size();

	auto fun = [&](real mu) -> real {
		return -log_mean_posterior_eval<real>(mu, lX, lx_sum, l0_min, l0_max,
		                                      l1_min, l1_max, 0);
	};

	/* Evaluate the posterior at the mean: */
	const real mu0 = x_sum / N;
	const real y_mu0 = -fun(mu0);

	/* Try to find the mu at which the posterior has fallen significantly: */
	real mu_r = 10.0 * mu0;
	while (-fun(mu_r) > y_mu0 - 20.0)
	{
		mu_r *= 10.0;
	}

	/* Within the bracket [0, mu_r], find the maximum: */
	constexpr int bits = std::numeric_limits<real>::digits;
	std::uintmax_t max_iter = 50;
	std::pair<real,real>
	   maxpost = brent_find_minima(fun, static_cast<real>(0.0), mu_r, bits,
	                               max_iter);

	/* Return the log of the posterior: */
	maxpost.second = -maxpost.second;
	return maxpost;
}


void LogNormalPosterior::log_mean_posterior(const size_t M, const double* mu,
                                            double* log_posterior) const
{
	typedef math<long double> mth;
	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 (l0_min == l0_max) not "
		                       "implemented.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 (l1_min == l1_max) not "
		                       "implemented.");
	}

	/* Compute the normalization constant: */
	compute_lI();
	if (lmp_log_norm == -1.0){
		lmp_log_norm = lI + lga - mth::log(2.0);
	}

	/* Now evaluate posterior: */
	#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		log_posterior[i]
		   = log_mean_posterior_eval<long double>(mu[i], lX, lx_sum,
		                                          prior.l0_min, prior.l0_max,
		                                          prior.l1_min, prior.l1_max,
		                                          lmp_log_norm);
	}

}


LogNormalPosterior::post_pred_t::post_pred_t(size_t N, long double dlga)
   : lX_xi(N+1), dlga(dlga)
{}

LogNormalPosterior::post_pred_t
LogNormalPosterior::predictive_parameters() const
{
	const size_t N = lX.size();
	post_pred_t params(N, math<long double>::lgamma(0.5 * (N + 1.0) - 0.5)
	                      - math<long double>::lgamma(0.5 * N - 0.5));
	std::copy(lX.cbegin(), lX.cend(), params.lX_xi.begin());

	return params;
}


void LogNormalPosterior::log_posterior_predictive(const size_t M,
                                                  const double* x,
                                                  double* log_post_pred) const
{
	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 (l0_min == l0_max) not "
		                       "implemented.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 (l1_min == l1_max) not "
		                       "implemented.");
	}

	/* Normalization constant: */
	constexpr long double ln_sqrt_2pi
	   = 0.5 * std::log((long double)2.0 * std::numbers::pi_v<long double>);
	compute_lI();
	post_pred_t pred_params = predictive_parameters();

	/* Evaluate posterior: */
	#pragma omp parallel for firstprivate(pred_params)
	for (size_t i=0; i<M; ++i){
		const double xi = x[i];
		if (xi <= 0.0)
			log_post_pred[i] = -std::numeric_limits<double>::infinity();
		else {
			try {
				/* Create the joint set of lX and the logarithm of the
				 * evaluation point: */
				const double lxi = std::log(xi);
				pred_params.lX_xi.back() = lxi;

				const long double lI_xi
				   = compute_log_integral_I<long double>(prior.l0_min,
				                                         prior.l0_max,
				                                         prior.l1_min,
				                                         prior.l1_max,
				                                         pred_params.lX_xi);
				log_post_pred[i] = -ln_sqrt_2pi - lxi + pred_params.dlga
				                   + lI_xi - lI;
			} catch (...) {
				log_post_pred[i] = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}
}


void
LogNormalPosterior::init_predictive_cumulative_interpolator(
    LogNormalPosterior::cumulative_t which
) const
{
	typedef math<long double> mth;

	if ((!predictive_cdf_interpolator && which == CDF) ||
	    (!predictive_ccdf_interpolator && which == CCDF))
	{
		/*
		 * Normalization constant:
		 */
		constexpr long double ln_sqrt_2pi
		   = 0.5 * std::log((long double)2.0 * std::numbers::pi_v<long double>);
		compute_lI();

		/* Compute the median (easy since lX is sorted): */
		long double center = mth::exp(lX[lX.size() / 2]);

		/* Setup a transform that spreads the posterior predictive CDF
		 * quite evenly in the interval [-1,1]:
		 */
		pred_trans_t transform(center);

		/* The PDF we integrate: */
		auto pdf = [&](const long double x) -> long double {
			if (x <= 0.0)
				return 0.0;
			else if (std::isinf(x))
				return 0.0;
			else {
				post_pred_t pred_params = predictive_parameters();
				try {
					/* Create the joint set of lX and the logarithm of the
					 * evaluation point: */
					const long double lx = mth::log(x);
					pred_params.lX_xi.back() = lx;


					const long double lI_xi
					   = compute_log_integral_I<long double>(prior.l0_min,
					                                         prior.l0_max,
					                                         prior.l1_min,
					                                         prior.l1_max,
					                                         pred_params.lX_xi);
					return mth::exp(-ln_sqrt_2pi - lx + pred_params.dlga + lI_xi
						            - lI);
				} catch (...) {
					return std::numeric_limits<long double>::quiet_NaN();
				}
			}
		};

		/* Compute the CDF interpolator: */
		if (which == CDF)
			predictive_cdf_interpolator = pred_cdf_interp_t(pdf, transform,
			                                                n_chebyshev);
		else
			predictive_ccdf_interpolator = pred_ccdf_interp_t(pdf, transform,
			                                                  n_chebyshev);
	}
}


void LogNormalPosterior::posterior_predictive_cdf(const size_t M,
                                                  const double* x,
                                                  double* post_pred_cdf) const
{
	typedef math<long double> mth;

	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 (l0_min == l0_max) not "
		                       "implemented.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 (l1_min == l1_max) not "
		                       "implemented.");
	}

	/* Ensure that we have the CDF interpolator: */
	init_predictive_cumulative_interpolator(CDF);

	/* Interpolate: */
	std::exception except;
	bool have_exception = false;
	const pred_cdf_interp_t& cdf(*predictive_cdf_interpolator);
	#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		if (have_exception)
			continue;
		try {
			post_pred_cdf[i] = cdf(x[i]);
		} catch (const std::exception& e) {
			have_exception = true;
			except = e;
		}
	}
	if (have_exception)
		throw except;
}


void LogNormalPosterior::posterior_predictive_ccdf(const size_t M,
                                                   const double* x,
                                                   double* post_pred_ccdf) const
{
	typedef math<long double> mth;

	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 (l0_min == l0_max) not "
		                       "implemented.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 (l1_min == l1_max) not "
		                       "implemented.");
	}

	/* Ensure that we have the CDF interpolator: */
	init_predictive_cumulative_interpolator(CCDF);

	/* Interpolate: */
	std::exception except;
	bool have_exception = false;
	const pred_ccdf_interp_t& ccdf(*predictive_ccdf_interpolator);
	#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		if (have_exception)
			continue;
		try {
			post_pred_ccdf[i] = ccdf(x[i]);
		} catch (const std::exception& e) {
			have_exception = true;
			except = e;
		}
	}
	if (have_exception)
		throw except;
}


void LogNormalPosterior::posterior_predictive_quantiles(const size_t M,
                                                  const double* q,
                                                  double* post_pred_q) const
{
	typedef math<long double> mth;

	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 (l0_min == l0_max) not "
		                       "implemented.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 (l1_min == l1_max) not "
		                       "implemented.");
	}

	/* Ensure that we have the CDF interpolator: */
	init_predictive_cumulative_interpolator(CDF);

	/* Interpolate: */
	std::exception except;
	bool have_exception = false;
	const pred_cdf_interp_t& cdf(*predictive_cdf_interpolator);
	#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		if (have_exception)
			continue;
		try {
			post_pred_q[i] = cdf.quantile(q[i]);
		} catch (const std::exception& e) {
			std::cout << "exception: " << e.what() << "\n" << std::flush;
			have_exception = true;
			except = e;
		}
	}
	if (have_exception)
		throw except;
}


void LogNormalPosterior::posterior_predictive_tail_quantiles(const size_t M,
                                                  const double* q,
                                                  double* post_pred_tq) const
{
	typedef math<long double> mth;

	/* Sanity check: */
	if (prior.l0_min == prior.l0_max){
		/* This is the case of known l0. Not implemented. */
		throw std::logic_error("Known l0 (l0_min == l0_max) not "
		                       "implemented.");
	}
	if (prior.l1_min == prior.l1_max){
		/* This is the case of known l1. Not implemented. */
		throw std::logic_error("Known l1 (l1_min == l1_max) not "
		                       "implemented.");
	}

	/* Ensure that we have the CDF interpolator: */
	init_predictive_cumulative_interpolator(CCDF);

	/* Interpolate: */
	std::exception except;
	bool have_exception = false;
	const pred_ccdf_interp_t& ccdf(*predictive_ccdf_interpolator);
	#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		if (have_exception)
			continue;
		try {
			post_pred_tq[i] = ccdf.quantile(q[i]);
		} catch (const std::exception& e) {
			have_exception = true;
			except = e;
		}
	}
	if (have_exception)
		throw except;
}
