/*
 * Fit_Polynom.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 29.10.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef FIT_POLYNOM_HH_
#define FIT_POLYNOM_HH_


void Fit_Polynom(double *x, double *y, int Startindex, int Endindex, double x0,
		int Polynomgrad, double *Parametervektor);
void x_zu_Minimum_von_y_in_Intervall(double *x, double *y, int Startindex,
		int Endindex, double &x_min, double &y_min, int &Indexmin);

#endif /* FIT_POLYNOM_HH_ */
