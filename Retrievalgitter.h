/*
 * Retrievalgitter.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 26.05.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef RETRIEVALGITTER_HH_
#define RETRIEVALGITTER_HH_

#include <vector>
#include "Gitterpunkt.h"

class Retrievalgitter
{
public:
	//Methoden
	void Retrievalgitter_erzeugen(std::vector<class Ausgewertete_Messung_Limb> &AM_Limb,
			double Epsilon, class Konfiguration &Konf);
	void Alle_Durchstosspunkte_Null_setzen();
	void In_Datei_Ausgeben(std::string Dateiname);

// //TODO löschen nach Implementierung im Raytracing
// // wird bei der Berechnung der AMF Matrix für Zwischenrechnungen benötigt
// void SZA_und_Streuwinkel_NULL_setzen();
// //TODO löschen nach Implementierung im Raytracing
// void SZA_und_Streuwinkel_bestimmen(double Sat_Lat,
//   double Sat_Lon, double Sat_Hoehe,
//   //Richtpunkt (Tangentenpunkt oder Grundpunkt)
//   double RP_Lat,double RP_Lon, double RP_Hoehe,
//   double Erdradius, double Deklination, double Stunde, double Minute);
//   // für jeden Strahlengang neu
// //TODO löschen nach test
// void Retrievalgitter_erzeugen_und_Messrichtung_herausfinden(
//   vector<Ausgewertete_Messung_Limb>& AM_Limb);

	// Hilfsfunktionen
	//Membervariablen
	int m_Anzahl_Hoehen;
	int m_Anzahl_Breiten;
	int m_Anzahl_Punkte;
	std::vector<Gitterpunkt> m_Gitter;
};
//Hilfsfunktionen
int Get_Index_of_Maximum(std::vector<double> A);
int Get_Index_of_Minimum(std::vector<double> A, int Startindex);


#endif /* RETRIEVALGITTER_HH_ */
