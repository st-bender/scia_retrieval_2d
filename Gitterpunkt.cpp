/*
 * Gitterpunkt.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 27.05.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
#include "Gitterpunkt.h"

// Konstruktoren
Gitterpunkt::Gitterpunkt() :
	m_vorderer_Durchstosspunkt(3), m_hinterer_Durchstosspunkt(3)
{
	//leer
	// sets only the length of the intersection point vectors
}
//Methoden
bool Gitterpunkt::Punkt_in_Gitterpunkt(double Lat, double Hoehe)
{
	return ((Lat <= m_Max_Breite) && (Lat > m_Min_Breite)
			&& (Hoehe <= m_Max_Hoehe) && (Hoehe > m_Min_Hoehe));
}
