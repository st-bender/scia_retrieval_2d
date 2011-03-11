/*
 * Ausgewertete_Messung_Nadir.h
 *
 *  Created on: 30.04.2010
 *      Author: martin
 */

#ifndef AUSGEWERTETE_MESSUNG_NADIR_HH_
#define AUSGEWERTETE_MESSUNG_NADIR_HH_

#include <string>
#include <cstdio>
#include <iostream>


using namespace std;

class Ausgewertete_Messung_Nadir
{
public:
	inline void Ausgabe_auf_Bildschirm();
	// Membervariablen
	//Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	double m_Deklination;
	double m_Sonnen_Longitude;
	// Zusatzinformationen
	double m_Wellenlaenge; //Übergangswellenlänge
	// Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	// Geolokationen
	double m_Lattitude_Sat;
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Lattitude_Ground;
	double m_Longitude_Ground;
	double m_Erdradius;
};
inline void Ausgewertete_Messung_Nadir::Ausgabe_auf_Bildschirm()
{
	// Ergebnisse
	cout << "m_Zeilendichte           : " << "\t" << m_Zeilendichte << "\n";
	cout << "m_Fehler_Zeilendichten   : " << "\t" << m_Fehler_Zeilendichten << "\n";
	// Zwischenergebnisse
	cout << "m_Deklination            : " << "\t" << m_Deklination << "\n";
	cout << "m_Sonnen_Longitude       : " << "\t" << m_Sonnen_Longitude << "\n";
	// Zusatzinfo
	cout << "Wellenlänge des Übergangs: " << "\t" << m_Wellenlaenge << "\n";
	// Datum
	cout << "m_Jahr: " << "\t" << m_Jahr << "\n";
	cout << "m_Monat: " << "\t" << m_Monat << "\n";
	cout << "m_Tag: " << "\t" << m_Tag << "\n";
	cout << "m_Stunde: " << "\t" << m_Stunde << "\n";
	cout << "m_Minute: " << "\t" << m_Minute << "\n";
	// Geolokationen
	cout << "m_Lattitude_Sat: " << "\t" << m_Lattitude_Sat << "\n";
	cout << "m_Longitude_Sat: " << "\t" << m_Longitude_Sat << "\n";
	cout << "m_Hoehe_Sat: " << "\t" << m_Hoehe_Sat << "\n";
	cout << "m_Lattidude_Ground: " << "\t" << m_Lattitude_Ground << "\n";
	cout << "m_Longitude_Ground: " << "\t" << m_Longitude_Ground << "\n";
	cout << "m_Erdradius: " << "\t" << m_Erdradius << "\n";
}
#endif /* AUSGEWERTETE_MESSUNG_NADIR_HH_ */
