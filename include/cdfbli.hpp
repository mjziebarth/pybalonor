/*
 * Barycentric Lagrance interpolation of the CDF.
 */

#include <stdexcept>
#include <numbers>
#include <vector>
#include <optional>
#include <math.hpp>
#include <iostream>

#define BOOST_ENABLE_ASSERT_HANDLER
#include <boost/math/quadrature/gauss_kronrod.hpp>


#ifndef INDIAPALEALE_CDFBLI_HPP
#define INDIAPALEALE_CDFBLI_HPP

namespace indiapaleale {


template<typename real>
class CenteredExpTanTransform {
public:
	typedef real real_t;
	typedef math<real> mth;

	CenteredExpTanTransform(real center) : center(center)
	{};

	real operator()(real x) const {
		constexpr real half_pi = std::numbers::pi_v<real> / 2;
		if (x < -1.0 || x > 1.0)
			throw std::domain_error("Input coordinate out of domain in "
			                        "CenteredExpTanTransform.");
		if (x == -1.0)
			return 0.0;
		else if (x == 1.0)
			return std::numeric_limits<real>::infinity();
		return center * mth::exp(mth::tan(half_pi * x));
	};

	/*
	 * Inverse transform:
	 */
	real inv(real x) const {
		constexpr real half_pi = std::numbers::pi_v<real> / 2;
		if (x < 0){
			std::string msg("Transformed coordinate ");
			msg.append(std::to_string(x));
			msg.append(" out of transformed domain.");
			throw std::domain_error(msg);
		} else if (x == 0)
			return -1.0;
		else if (std::isinf(x))
			return 1.0;
		return std::max<real>(std::min<real>(
		          mth::atan(mth::log(x / center)) / half_pi, 1.0), -1.0);
	};

	/*
	 * The input domain of this transform:
	 */
	constexpr real domain_min() const {
		return -1.0;
	};
	constexpr real domain_max() const {
		return 1.0;
	};

private:
	real center;
};


template<typename transform_t, bool complementary=false>
class BarycentricLagrangeCDFInterpolator {
public:
	/*
	 * Obtain the floating point type from the transform type:
	 */
	typedef typename transform_t::real_t real;

	template<typename F>
	BarycentricLagrangeCDFInterpolator(F pdf, transform_t transform,
	                                   size_t n_chebyshev)
	   : a(transform.domain_min()), b(transform.domain_max()), trans(transform),
	     support(n_chebyshev, zf_t({0.0, 0.0}))
	{
		typedef boost::math::quadrature::gauss_kronrod<real, 7> GK7;
		typedef boost::math::quadrature::gauss_kronrod<real, 15> GK15;

		/* Sanity checks: */
		if (n_chebyshev <= 1)
			throw std::runtime_error("Need at least 2 Chebyshev points.");
		if (a >= b)
			throw std::runtime_error("a < b is required.");

		/* Prepare the interpolation points: */
		std::vector<real> x(n_chebyshev);
		x[0] = trans(b); // Ensure correct start of the interval.
		std::vector<real> z(n_chebyshev);
		z[0] = 1.0;
		for (size_t i=1; i<n_chebyshev-1; ++i){
			constexpr real pi = std::numbers::pi_v<real>;
			const real zi = std::cos(i * pi / (n_chebyshev-1));
			/*
			 * Generate Chebyshev points in the input domain of the
			 * transform, and transform them to coordinates of the PDF:
			 */
			x[i] = trans(std::min(std::max<real>(  0.5 * (1.0 - zi) * a
			                                     + 0.5 * (1.0 + zi) * b, a),
			                      b));
			z[i] = zi;
		}
		x.back() = trans(a); // Ensure correct end of the interval.
		z.back() = -1.0;

		/* Evaluate the CDF.
		 * Integrate the given PDF in the intervals given by the interval
		 * [a,b] split by the Chebyshev nodes.
		 * The integrals are performed in parallel using OpenMP.
		 * The estimation process is two-step:
		 *   1) Get a first quick estimate of the CDF increments using
		 *      a GK7 rule.
		 *   2) For all intervals in which the resulting error estimate
		 *      is larger than sqrt(epsilon), refine the integral using
		 *      an adaptive Gauss-Kronrod rule.
		 */
		struct result_t {
			real inc;
			real err;
		};
		std::vector<result_t> increments(n_chebyshev-1);
		std::exception except;
		bool have_exception = false;
		/*
		 * Step 1.
		 */
		#pragma omp parallel for
		for (size_t i=0; i<n_chebyshev-1; ++i){
			/* Skip if exception occurred: */
			if (have_exception)
				continue;
			try {
				/* Note: Chebyshev nodes above decrease with i, hence
				 *       x[i+1] is left and x[i] is right boundary.
				 */
				increments[i].inc =
				    GK7::integrate(
				       pdf,
				       x[i+1], // left integration boundary
				       x[i],   // right boundary
				       0,      // max number of descent levels
				       0.0,    // tolerance (not used)
				       &increments[i].err
				);
			} catch (const std::exception& e) {
				have_exception = true;
				except = e;
			}
		}
		if (have_exception)
			throw except;

		/*
		 * Step 2: Refine intervals of large error:
		 */
		constexpr real root_eps
		   = boost::math::tools::root_epsilon<real>();
		std::vector<size_t> to_refine;
		for (size_t i=0; i<n_chebyshev-1; ++i){
			if (increments[i].err > root_eps)
				to_refine.push_back(i);
		}
		#pragma omp parallel for schedule(dynamic,1)
		for (size_t i : to_refine){
			try {
				const real tol = std::min<real>(root_eps / increments[i].inc,
				                                1.0);
				increments[i].inc =
				    GK15::integrate(
				        pdf,
				        x[i+1], // left integration boundary
				        x[i],   // right boundary
				        7,      // max number of descent levels
				        tol,    // tolerance
				        &increments[i].err
				);
			} catch (const std::exception& e) {
				have_exception = true;
				except = e;
			}
		}
		if (have_exception)
			throw except;



		/* Cumulative sum and transfer to support vector.
		 * Note here again the descending property of Chebyshev nodes with
		 * index i.
		 */
		real cdf = 0.0;
		real cumul_err = 0.0;
		if (complementary){
			/* Cumulative sum from the back: */
			support[0].z = z[0];
			support[0].f = 0.0;
			for (size_t i=0; i<n_chebyshev-1; ++i){
				cdf += increments[i].inc;
				cumul_err += increments[i].err;
				support[i+1].z = z[i];
				support[i+1].f = cdf;
			}
		} else {
			/* Cumulative sum from the front: */
			support.back().z = z.back();
			support.back().f = 0.0;
			for (size_t i=0; i<n_chebyshev-1; ++i){
				const size_t j = n_chebyshev-i-2;
				cdf += increments[j].inc;
				cumul_err += increments[j].err;
				support[j].z = z[j];
				support[j].f = cdf;
			}
		}
	};

	real operator()(real x) const {
		/* Transform the coordinate to the input domain of the transform: */
		real z = trans.inv(x);

		/* Perform Barycentric Lagrange interpolation: */
		return bli(z);
	}


	real quantile(real q) const {
		/*
		 * This method computes quantiles.
		 */

		/* Sanity: */
		if (q < 0.0 || q > 1.0)
			throw std::domain_error("Requested quantile outside [0,1].");

		/* Now solve for the quantile: */
		if (q == support.front().f)
			return support.front().z;
		else if (q == support.back().f)
			return support.back().z;
		else {
			/* The typical case. Use TOMS 748 on a quantile
			 * function to find the quantile.
			 */
			auto quantile_function = [&](real z) -> real {
				return bli(z) - q;
			};
			std::uintmax_t max_iter(100);
			boost::math::tools::eps_tolerance<real>
			   eps_tol(std::numeric_limits<real>::digits - 2);
			const real left_qf = (complementary) ? 1.0 - q : -q;
			const real right_qf = (complementary) ? -q : 1.0 - q;
			std::pair<real,real> bracket
			   = toms748_solve(quantile_function, a, b, left_qf, right_qf,
			                   eps_tol, max_iter);

			/* Bracket is in the input domain of the transform.
			 * Need to transform back to the input domain of the PDF
			 * (the output domain of the transform):
			 */
			return trans(0.5*(bracket.first + bracket.second));
		}
	}

private:
	/*
	 * Types:
	 */
	struct zf_t {
		real z;
		real f;
	};

	/*
	 * Members:
	 */
	real a;
	real b;
	transform_t trans;
	std::vector<zf_t> support;


	real bli(real z) const {
		/*
		 * This function evaluates the CDF using Barycentric Lagrange
		 * interpolation of the support vector.
		 * Reference:
		 *    Berrut, J.-P. and Trefethen, Lloyd N. (2004): Barycentric
		 *    Lagrange Interpolation. SIAM Review 46(3), 501-517.
		 *    https://dx.doi.org/10.1137/S0036144502417715
		 */
		auto xfit = support.cbegin();
		if (z == xfit->z)
			return xfit->f;
		real nom = 0.0;
		real denom = 0.0;
		real wi = 0.5 / (z - xfit->z);
		nom += wi * xfit->f;
		denom += wi;
		int8_t sign = -1;
		++xfit;
		for (size_t i=1; i<support.size()-1; ++i){
			if (z == xfit->z)
				return xfit->f;
			wi = sign * 1.0 / (z - xfit->z);
			nom += wi * xfit->f;
			denom += wi;
			sign = -sign;
			++xfit;
		}
		if (z == xfit->z)
			return xfit->f;
		wi = sign * 0.5 / (z - xfit->z);
		nom += wi * xfit->f;
		denom += wi;
		return std::max<real>(std::min<real>(nom / denom, 1.0), 0.0);
	}
};

} // end namespace

#endif