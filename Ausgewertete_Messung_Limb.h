/*
 * Ausgewertete_Datei_Limb.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 19.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
// Dieses Objekt enthält alle Informationen einer Limbmessung, die nach der
// Bestimmung der Zeilendichte noch benötigt werden
// Für den Aufbau des Gitters der Teilintegranden AMF
// geographische Breite des Tangentenpunkts
// geographische Länge des Tangentenpunkts
// Höhe am Tangentenpunkt
// Erdradius am Tangentenpunkt
// geographische Breite des Satelliten, für Nord-Süd Unterscheidung

// Für den Deklinationswinkel
// Jahr, Monat, Tag, Stunde, Minute
// Funktion zur Brechnung von Deklinationswinkel, Variable Deklinationswinkel

// Die ausgewertete Zeilendichte selbst
#ifndef AUSGEWERTETE_MESSUNG_LIMB_HH_
#define AUSGEWERTETE_MESSUNG_LIMB_HH_

#include <string>
#include <cstdio>
#include <iostream>

class Ausgewertete_Messung_Limb
{
public:
	inline void Ausgabe_auf_Bildschirm();
	// Membervariablen
	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	// total number density at measurement point
	double total_number_density;
	// Zwischenergebnisse
	double m_Deklination;
	double m_Sonnen_Longitude;
	double m_LocalSolarTime;
	// Wellenlänge des Übergangs
	double m_Wellenlaenge;
	double m_Wellenlaenge_abs;
	//Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;        //Uhrzeit ist wichtig für Längengrad der Sonne
	int m_Stunde;
	int m_Minute;
	int m_Sekunde;
	// Geolocation
	double m_Latitude_Sat;
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Latitude_TP;
	double m_Longitude_TP;
	double m_Hoehe_TP;
	double m_Erdradius;
	// phase of orbit (0...1)
	double m_orbit_phase;
	double center_lat, center_lon;
};

// inline Implementierungen
inline void Ausgewertete_Messung_Limb::Ausgabe_auf_Bildschirm()
{
	std::cout << "m_Zeilendichte           : " << "\t" << m_Zeilendichte << std::endl;
	std::cout << "m_Fehler_Zeilendichten   : " << "\t" << m_Fehler_Zeilendichten << std::endl;
	std::cout << "total_number_density     : " << "\t" << total_number_density << std::endl;
	std::cout << "m_Deklination            : " << "\t" << m_Deklination << std::endl;
	std::cout << "m_Sonnen_Longitude       : " << "\t" << m_Sonnen_Longitude << std::endl;
	std::cout << "m_LocalSolarTime         : " << "\t" << m_LocalSolarTime << std::endl;
	std::cout << "Wellenlänge des Übergangs: " << "\t" << m_Wellenlaenge << std::endl;
	std::cout << "Jahr                     : " << "\t" << m_Jahr << std::endl;
	std::cout << "Monat                    : " << "\t" << m_Monat << std::endl;
	std::cout << "Tag                      : " << "\t" << m_Tag << std::endl;
	std::cout << "Stunde                   : " << "\t" << m_Stunde << std::endl;
	std::cout << "Minute                   : " << "\t" << m_Minute << std::endl;
	std::cout << "Sekunde                  : " << "\t" << m_Sekunde << std::endl;
	std::cout << "m_Latitude_Sat           : " << "\t" << m_Latitude_Sat << std::endl;
	std::cout << "m_Longitude_Sat          : " << "\t" << m_Longitude_Sat << std::endl;
	std::cout << "m_Hoehe_Sat              : " << "\t" << m_Hoehe_Sat << std::endl;
	std::cout << "m_Latitude_TP            : " << "\t" << m_Latitude_TP << std::endl;
	std::cout << "m_Longitude_TP           : " << "\t" << m_Longitude_TP << std::endl;
	std::cout << "m_Hoehe_TP               : " << "\t" << m_Hoehe_TP << std::endl;
	std::cout << "m_Erdradius              : " << "\t" << m_Erdradius << std::endl;
	std::cout << "m_orbit_phase            : " << "\t" << m_orbit_phase << std::endl;
	std::cout << "center_lat               : " << "\t" << center_lat << std::endl;
	std::cout << "center_lon               : " << "\t" << center_lon << std::endl;
}


#endif /* AUSGEWERTETE_MESSUNG_LIMB_HH_ */
