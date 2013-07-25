/*
 * Gitterpunkt.cpp
 *
 *  Created on: 27.05.2010
 *      Author: martin
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
