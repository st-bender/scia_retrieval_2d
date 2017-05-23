/*
 * Koordinatentransformation.h
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

#ifndef KOORDINATENTRANSFORMATION_HH_
#define KOORDINATENTRANSFORMATION_HH_

// Funktionen zur Umwandlung von karthesichen Koordinaten
// in andere Systeme und umgekehrt

// Kugelkoordinaten
void Umwandlung_Kugel_in_Karthesisch(double r, double phi, double theta,
		double &x, double &y, double &z);
void Umwandlung_Karthesisch_in_Kugel(double x, double y, double z, double &r,
		double &phi, double &theta);

#endif /* KOORDINATENTRANSFORMATION_HH_ */
