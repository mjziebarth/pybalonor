/*
 * Template math functions.
 */

#include <math.hpp>
#include <math.h>

using indiapaleale::math;

long double math<long double>::exp(long double x){
	return expl(x);
}

long double math<long double>::log(long double x){
	return logl(x);
}

long double math<long double>::tan(long double x){
	return tanl(x);
}

long double math<long double>::atan(long double x){
	return atanl(x);
}

long double math<long double>::sqrt(long double x){
	return sqrtl(x);
}

long double math<long double>::lgamma(long double x){
	return lgammal(x);
}

