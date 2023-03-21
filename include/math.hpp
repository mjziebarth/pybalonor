/*
 * Template math functions.
 */

#ifndef INDIAPALEALE_MATH_HPP
#define INDIAPLAEALE_MATH_HPP

namespace indiapaleale{

template<typename real>
struct math;


template<>
struct math<long double>
{
	static long double exp(long double x);
	static long double log(long double x);
	static long double tan(long double x);
	static long double atan(long double x);
	static long double sqrt(long double x);
	static long double lgamma(long double x);
};

}

#endif