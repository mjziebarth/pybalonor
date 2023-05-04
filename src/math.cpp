/*
 * SPDX-License-Identifier: EUPL-1.2
 *
 * Template math functions.
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

