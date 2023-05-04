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