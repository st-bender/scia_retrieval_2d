/*
 * Gitterpunkt.h
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

#ifndef GITTERPUNKT_HH_
#define GITTERPUNKT_HH_

#include "MPL_Vektor.h"

class Gitterpunkt
{
public:
	// Konstruktoren
	Gitterpunkt();
	// Methoden
	bool Punkt_in_Gitterpunkt(double Lat, double Hoehe); //liefert true falls ja, false falls nein
	//int Sonnenzenitwinkel_der_Gitterpunkte_berechnen( diverse andere Winkel);
	// räumliche Lokalisierung des Gitterpunkts
	double m_Max_Hoehe;
	double m_Hoehe;
	double m_Min_Hoehe;
	double m_Max_Breite;
	double m_Breite;
	double m_Min_Breite;
	// longitude
	double longitude;
	// Indizierung
	int m_eigener_Index;
	// Indizierung Nachbarpunkte
	int m_Index_oberer_Nachbar ; // bzw. Hoeher
	int m_Index_unterer_Nachbar; // bzw. Niedriger
	int m_Index_Nord_Nachbar;
	int m_Index_Sued_Nachbar;
	// Die Nachbarn auf den 4 Diagonalen Punkten
	// (später Praktisch für das Raytracing)
	int m_Index_oberer_Nord_Nachbar ;
	int m_Index_unterer_Nord_Nachbar;
	int m_Index_oberer_Sued_Nachbar ;
	int m_Index_unterer_Sued_Nachbar;

	// Die beiden Durchstoßpunkte des LOS Raytracings sollte man sich merken....
	// im LFS Retrieval für Longitude des Startpunkts wichtig
	MPL_Vektor m_vorderer_Durchstosspunkt;  // vorne beim satelliten
	MPL_Vektor m_hinterer_Durchstosspunkt;  //hinten weiter weg vom satelliten

	// Diese Winkel werden bei der Erzeugung der AMF Matrix
	// als Zwischenlösungen benötigt..... Todo Winkel in Raytracing
	//double m_SZA;                       // Sonnenzenitwinkel des Gitterpunktes
	//double m_Streuwinkel;               // Streuwinkel des Gitterpunkts
};
#endif /* GITTERPUNKT_HH_ */
