/*
 * Koordinatentransformation.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 08.09.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include <cmath>
#include<iostream>



////////////////////////////////////////////////////////////////////////////////
// Kugelkoordinaten
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////
// Funktionsstart int  Umwandlung_Kugel_in_Karthesisch(double r,double phi,
// double theta,double& x,double& y,double& z);
////////////////////////////////////////////////
void Umwandlung_Kugel_in_Karthesisch(double r, double phi, double theta,
		double &x, double &y, double &z)
{
	// Winkel in Grad
	const double deg = M_PI / 180.0;
	const double cth = cos(deg * theta), sth = sin(deg * theta);
	const double cph = cos(deg * phi), sph = sin(deg * phi);

	x = r * cth * cph;
	y = r * cth * sph;
	z = r * sth;
}
////////////////////////////////////////////////
// ENDE int  Umwandlung_Kugel_in_Karthesisch
////////////////////////////////////////////////

////////////////////////////////////////////////
// Funktionsstart int Umwandlung_Karthesisch_in_Kugel(double x,double y,
// double z,double& r,double& phi,double& theta);
////////////////////////////////////////////////
void Umwandlung_Karthesisch_in_Kugel(double x, double y, double z, double &r,
		double &phi, double &theta)
{
	// Winkel in Grad
	// Der Punkt 0,0,0 ist ausgeschlossen
	const double rad = 180.0 * M_1_PI;

	r = sqrt(x * x + y * y + z * z);
	if (r == 0.0) {
		std::cout << " Umwandlung in Kugelkoordinaten von (0,0,0) nicht sinnvoll\n";
		phi = 0;
		theta = 0;
		return;
	}
	theta = rad * asin(z / r);
	if ((theta >= 90.0) || (theta <= -90.0)) {
		// Phi quasi beliebig
		phi = 0;
		return;
	}

	phi = rad * atan2(y, x);
}
////////////////////////////////////////////////
// ENDE int Umwandlung_Karthesisch_in_Kugel
////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// ENDE  Kugelkoordinaten
////////////////////////////////////////////////////////////////////////////////
