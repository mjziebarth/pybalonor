#include <lognormal.hpp>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <numbers>
#include <stdexcept>
#include <iostream>

#define BOOST_ENABLE_ASSERT_HANDLER
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
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
using boost::math::tools::toms748_solve;
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::brent_find_minima;
using boost::math::tools::eps_tolerance;
using indiapaleale::LogNormalPosterior;

template<typename real>
struct math;

template<>
struct math<long double>
{
	static long double exp(long double x){
		return expl(x);
	}
	static long double log(long double x){
		return logl(x);
	}
	static long double sqrt(long double x){
		return sqrtl(x);
	}
	static long double lgamma(long double x){
		return lgammal(x);
	}
};

template<typename real>
struct lognormal_params_t {
	real l0;
	real l1;
};

/*
 * log-normal maximum likelihood.
 */
template<typename real>
lognormal_params_t<real> log_normal_mle(const std::vector<double>& lX)
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
                            const std::vector<double>& lx)
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
	} else  if (c == 0.0){
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

		auto integrand = [&](real l0) -> real {
			real C1 = 0.0;
			for (double lxi : lx){
				const real dlx = lxi - l0;
				C1 += dlx * dlx;
			}
			C1 *= 0.5;
			real log_integrand = -alpha * math<real>::log(C1)
			    + math<real>::log(boost::math::gamma_q(alpha, C1/d2));
			return math<real>::exp(log_integrand - log_scale);
		};

		tanh_sinh<real> integrator;
		if (a < l0_mle && b > l0_mle)
			I =   integrator.integrate(integrand, a, l0_mle)
			    + integrator.integrate(integrand, l0_mle, b);
		else {
			I = integrator.integrate(integrand, a, b);
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

		auto integrand = [&](real l0) -> real {
			real C1 = 0.0;
			for (double lxi : lx){
				const real dlx = lxi - l0;
				C1 += dlx * dlx;
			}
			C1 *= 0.5;
			std::cout << "l0 = " << l0 << "\n";
			std::cout << "   C1 = " << C1 << "\n";
			real log_integrand = -alpha * math<real>::log(C1)
			    + log_delta_gamma_div_loggamma(alpha, C1/c2, C1/d2);
			std::cout << "   log_integrand = " << log_integrand << "\n" << std::flush;
			std::cout << "   scaled =        " << log_integrand - log_scale << "\n";
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

	return math<real>::log(I) + log_scale;
}

static void sanity_check(const double l0_min, const double l0_max,
                         const double l1_min, const double l1_max)
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
}

int LogNormalPosterior::log_posterior(const size_t M, const double* l0,
                             const double* l1, double* log_post,
                             const size_t N, const double* X,
                             const double l0_min, const double l0_max,
                             const double l1_min, const double l1_max)
{
	/* Sanity check: */
	sanity_check(l0_min, l0_max, l1_min, l1_max);
	if (l0_min == l0_max){
		/* This is the case of known l0. Not implemented. */
		return -1;
	}
	if (l1_min == l1_max){
		/* This is the case of known l1. Not implemented. */
		return -1;
	}

	/* Logarithm of the data: */
	std::vector<double> lX(N);
	for (size_t i=0; i<N; ++i){
		lX[i] = std::log(X[i]);
	}

	/* Normalization constant: */
	constexpr long double ln2 = std::log((long double)2.0);
	const long double lI
	   = compute_log_integral_I<long double>(l0_min, l0_max, l1_min,
	                                         l1_max, lX);

	/* Evaluate posterior: */
	const double lga = std::lgamma(0.5 * N - 0.5);
	//#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		const double l0i = l0[i];
		const double l1i = l1[i];
		if (l0i < l0_min || l0i > l0_max || l1i < l1_min || l1i > l1_max)
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
	return 0;
}

int LogNormalPosterior::posterior(const size_t M, const double* l0,
                         const double* l1, double* post,
                         const size_t N, const double* X,
                         const double l0_min, const double l0_max,
                         const double l1_min, const double l1_max)
{
	int rc = log_posterior(M, l0, l1, post, N, X, l0_min, l0_max, l1_min,
	                       l1_max);
	if (rc)
		return rc;
	for (size_t i=0; i<M; ++i)
		post[i] = std::exp(post[i]);

	return 0;
}


//template<typename fun_t, typename T>
//T find_maximum(fun_t fun, T xmin, T xmax)
//{
////	constexpr double tol = std::sqrt(std::numeric_limits<T>::epsilon());
//	constexpr int bits = std::numeric_limits<T>::digits;
//	T fun_xmax, fun_xmin;
//	if (std::isinf(xmax)){
//		std::array<T,3> xprev;
//		std::array<T,3> fprev;
//		xprev[0] = xmin;
//		fprev[0] = fun(xmin);
//		xprev[1] = xmin;
//		fprev[1] = fun(xmin);
//		xprev[2] = (xmin == 0) ? 1.0 : 2 * xmin;
//		fprev[2] = fun(xprev[2]);
//		while (fprev[2] > fprev[1]){
//			xprev[0] = xprev[1];
//			fprev[0] = fprev[1];
//			xprev[1] = xprev[2];
//			fprev[1] = fprev[2];
//			xprev[2] *= 2;
//			fprev[2] = fun(xprev[2]);
//		}
//		xmin = xprev[0];
//		xmax = xprev[2];
//		fun_xmin = fprev[0];
//		fun_xmax = fprev[2];
//	} else {
//		fun_xmax = fun(xmax);
//		fun_xmin = fun(xmin);
//	}
//	std::cout << "xmin = " << xmin << "\n"
//	          << "xmax = " << xmax << "\n" << std::flush;
//	std::cout << std::setprecision(20);
//	std::pair<T,T> res
//	   = brent_find_minima([&](T x) -> T {return -fun(x);}, xmin, xmax, bits);
//	return res.first;
//}


int LogNormalPosterior::log_mean_posterior(const size_t M, const double* mu,
                          double* log_posterior,
                          const size_t N, const double* X,
                          const double l0_min, const double l0_max,
                          const double l1_min, const double l1_max
                          )
{
	typedef math<long double> mth;
	/* Sanity check: */
	sanity_check(l0_min, l0_max, l1_min, l1_max);
	if (l0_min == l0_max){
		/* This is the case of known l0. Not implemented. */
		return -1;
	}
	if (l1_min == l1_max){
		/* This is the case of known l1. Not implemented. */
		return -1;
	}

	/* Logarithm of the data: */
	std::vector<double> lX(N);
	long double lx_sum = 0.0;
	for (size_t i=0; i<N; ++i){
		lX[i] = std::log(X[i]);
		lx_sum += lX[i];
	}

	/* Normalization constant: */
	constexpr long double ln2 = std::log((long double)2.0);
	const long double lI
	   = compute_log_integral_I<long double>(l0_min, l0_max, l1_min,
	                                         l1_max, lX);

	/* Evaluate posterior: */
	const double lga = std::lgamma(0.5 * N - 0.5);
	//#pragma omp parallel for
	for (size_t i=0; i<M; ++i){
		const double mu_i = mu[i];
		const double ln_mu = std::log(mu_i);
		if (ln_mu <= l0_min){
			log_posterior[i] = -std::numeric_limits<double>::infinity();
			continue;
		}
		const long double dstar
		   = std::min(mth::sqrt(2.0 * (ln_mu - l0_min)), (long double)l1_max);
		if (dstar <= l1_min){
			log_posterior[i] = -std::numeric_limits<double>::infinity();
			continue;
		}
		/*
		 * Compute the maximum of the integrand:
		 */
		long double log_scale = 0.0;
		long double l1_peak = 0.0;
		bool max_at_boundary = false;
		{
			auto log_integrand = [&](long double l1) -> long double {
				if (l1 == 0.0)
					return -std::numeric_limits<long double>::infinity();
				const long double dlx = ln_mu - 0.5 * l1 * l1;
				long double S = 0.0;
				for (double lxi : lX){
					S += (lxi - dlx) * (lxi - dlx);
				}
				return -(N * mth::log(l1)) - 0.5 * S / (l1 * l1);
			};
			auto fun0 = [&](long double l1) -> long double {
				const long double dlx = (ln_mu - 0.5 * l1 * l1);
				long double S = 0.0;
				for (double lxi : lX){
					S += (lxi - dlx) * (lxi - dlx);
				}
				return N * (ln_mu - 1.0) - 0.5 * N * l1 * l1 - lx_sum
				       + S / (l1 * l1);
			};
			if ((fun0(l1_min) > 0) == (fun0(dstar) > 0)){
				/* No change in sign -> maximum at the boundary */
				const long double li_l = log_integrand(l1_min);
				const long double li_r = log_integrand(dstar);
				max_at_boundary = true;
				log_scale = std::max(li_l, li_r);
			} else {
				/*
				 * Use Newton-Raphson to compute the maximum of the
				 * log integrand.
				 */
				auto fun1 = [&](long double l1) -> long double {
					const long double dlx = (ln_mu - 0.5 * l1 * l1);
					long double S0 = 0.0;
					long double S1 = 0.0;
					for (double lxi : lX){
						long double di = (lxi - dlx) / l1;
						S0 += di * di;
						S1 += di;
					}
					return -(N * l1) - 2 * S0 / l1  + 2 * S1;
				};
				/* Initial guess: */
				if (2.0 * ln_mu - 2.0/N * lx_sum < 0.0)
					l1_peak = l1_min + 0.5 * (dstar - l1_min);
				else
					l1_peak = std::max(std::min(
					             mth::sqrt(2 * ln_mu - 2.0/N * lx_sum),
					             dstar), (long double)l1_min);
				const long double tol
				  = mth::sqrt(std::numeric_limits<long double>::epsilon());
				for (size_t j=0; j<100; ++j){
					long double dl1 = - fun0(l1_peak) / fun1(l1_peak);
					dl1 = std::min(std::max(dl1, -0.9 * (l1_peak - l1_min)),
					               0.9 * (dstar - l1_peak));
					l1_peak += dl1;
					if (std::abs(dl1) < tol * l1_peak)
						break;
				}
				log_scale = log_integrand(l1_peak);
			}
		}
		/* Now integrate: */
		auto integrand = [&](long double l1) -> long double {
			if (l1 == 0.0)
				return 0.0;
			const long double dlx = ln_mu - 0.5 * l1 * l1;
			long double S = 0.0;
			for (double lxi : lX){
				S += (lxi - dlx) * (lxi - dlx);
			}
			long double log = -(N * mth::log(l1)) - 0.5 * S / (l1 * l1);
			return mth::exp(log - log_scale);
		};
		tanh_sinh<long double> integrator;
		long double I_l1 = 0.0;
		try {
			if (max_at_boundary){
				I_l1 = integrator.integrate(integrand, l1_min, dstar);
			} else {
				I_l1 =   integrator.integrate(integrand, l1_min, l1_peak)
					   + integrator.integrate(integrand, l1_peak, dstar);
			}
			log_posterior[i] = ln2 - lI - lga + mth::log(I_l1) + log_scale;
		} catch (...) {
			log_posterior[i] = std::numeric_limits<double>::quiet_NaN();
		}
	}
	return 0;

}



int LogNormalPosterior::log_posterior_predictive(const size_t M,
                          const double* x, double* log_post_pred,
                          const size_t N, const double* X,
                          const double l0_min, const double l0_max,
                          const double l1_min, const double l1_max
                          )
{
	/* Sanity check: */
	sanity_check(l0_min, l0_max, l1_min, l1_max);
	if (l0_min == l0_max){
		/* This is the case of known l0. Not implemented. */
		return -1;
	}
	if (l1_min == l1_max){
		/* This is the case of known l1. Not implemented. */
		return -1;
	}

	/* Logarithm of the data: */
	std::vector<double> lX(N);
	for (size_t i=0; i<N; ++i){
		lX[i] = std::log(X[i]);
	}

	/* The following vector contains the joint set of data x_i and
	 * the evaluation point x (which can always be set to the */
	std::vector<double> lX_xi(N+1);
	std::copy(lX.cbegin(), lX.cend(), lX_xi.begin());

	/* Normalization constant: */
	constexpr long double ln_sqrt_2pi
	   = 0.5 * std::log((long double)2.0 * std::numbers::pi_v<long double>);
	const long double lI
	   = compute_log_integral_I<long double>(l0_min, l0_max, l1_min,
	                                         l1_max, lX);

	/* Evaluate posterior: */
	const long double dlga = math<long double>::lgamma(0.5 * (N + 1.0) - 0.5)
	                       - math<long double>::lgamma(0.5 * N - 0.5);
	#pragma omp parallel for firstprivate(lX_xi)
	for (size_t i=0; i<M; ++i){
		const double xi = x[i];
		if (xi <= 0.0)
			log_post_pred[i] = -std::numeric_limits<double>::infinity();
		else {
			/* Create the joint set of lX and the logarithm of the
			 * evaluation point: */
			const double lxi = std::log(xi);
			lX_xi.back() = lxi;

			const long double lI_xi 
			   = compute_log_integral_I<long double>(l0_min, l0_max, l1_min,
			                                         l1_max, lX_xi);
			log_post_pred[i] = -ln_sqrt_2pi - lxi + dlga + lI_xi - lI;
		}
	}
	return 0;
}
