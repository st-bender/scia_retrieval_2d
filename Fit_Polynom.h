/*
 * Fit_Polynom.h
 *
 *  Created on: 29.10.2010
 *      Author: martin
 */

#ifndef FIT_POLYNOM_HH_
#define FIT_POLYNOM_HH_


void Fit_Polynom(double *x, double *y, int Startindex, int Endindex, double x0,
		int Polynomgrad, double *Parametervektor);
int x_zu_Minimum_von_y_in_Intervall(double *x, double *y, int Startindex,
		int Endindex, double &x_min, double &y_min, int &Indexmin);

#endif /* FIT_POLYNOM_HH_ */
