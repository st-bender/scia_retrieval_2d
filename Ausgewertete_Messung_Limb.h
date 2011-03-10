/*
 * Ausgewertete_Datei_Limb.h
 *
 *  Created on: 19.04.2010
 *      Author: martin
 */
// Dieses Objekt enthält alle Informationen einer Limbmessung, die nach der Bestimmung der Zeilendichte noch
// benötigt werden
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
#include <stdio.h>
#include <iostream>

using namespace std;

class Ausgewertete_Messung_Limb
{
public:
	inline void Ausgabe_auf_Bildschirm();
	// Membervariablen
	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	double m_Deklination;
	double m_Sonnen_Longitude;
	// Wellenlänge des Übergangs
	double m_Wellenlaenge;
	//Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;        //Uhrzeit ist wichtig für Längengrad der Sonne
	int m_Stunde;
	int m_Minute;
	// Geolocation
	double m_Lattidude_Sat;
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Lattidude_TP;
	double m_Longitude_TP;
	double m_Hoehe_TP;
	double m_Erdradius;
};

// inline Implementierungen
inline void Ausgewertete_Messung_Limb::Ausgabe_auf_Bildschirm()
{
	cout << "m_Zeilendichte                    : " << "\t" << m_Zeilendichte << "\n";
	cout << "m_Fehler_Zeilendichten       : " << "\t" << m_Fehler_Zeilendichten << "\n";
	cout << "m_Deklination                     : " << "\t" << m_Deklination << "\n";
	cout << "m_Sonnen_Longitude                     : " << "\t" << m_Sonnen_Longitude << "\n";
	cout << "Wellenlänge des Übergangs: " << "\t" << m_Wellenlaenge << "\n";
	cout << "Jahr                                    : " << "\t" << m_Jahr << "\n";
	cout << "Monat                                    : " << "\t" << m_Monat << "\n";
	cout << "Tag                                    : " << "\t" << m_Tag << "\n";
	cout << "Stunde                               : " << "\t" << m_Stunde << "\n";
	cout << "Minute                                : " << "\t" << m_Minute << "\n";
	cout << "m_Lattidude_Sat                 : " << "\t" << m_Lattidude_Sat << "\n";
	cout << "m_Longitude_Sat                : " << "\t" << m_Longitude_Sat << "\n";
	cout << "m_Hoehe_Sat                     : " << "\t" << m_Hoehe_Sat << "\n";
	cout << "m_Lattidude_TP                  : " << "\t" << m_Lattidude_TP << "\n";
	cout << "m_Longitude_TP                 : " << "\t" << m_Longitude_TP << "\n";
	cout << "m_Hoehe_TP                      : " << "\t" << m_Hoehe_TP << "\n";
	cout << "m_Erdradius                       : " << "\t" << m_Erdradius << "\n";
}


#endif /* AUSGEWERTETE_MESSUNG_LIMB_HH_ */
