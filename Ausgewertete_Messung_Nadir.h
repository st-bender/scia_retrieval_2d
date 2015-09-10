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
	double m_Wellenlaenge_abs; //Absorptionswellenlänge
	// Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	// Geolokationen
	double m_Latitude_Sat;
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Latitude_Ground;
	double m_Longitude_Ground;
	double m_Erdradius;
};
inline void Ausgewertete_Messung_Nadir::Ausgabe_auf_Bildschirm()
{
	// Ergebnisse
	std::cout << "m_Zeilendichte           : " << "\t" << m_Zeilendichte << std::endl;
	std::cout << "m_Fehler_Zeilendichten   : " << "\t" << m_Fehler_Zeilendichten << std::endl;
	// Zwischenergebnisse
	std::cout << "m_Deklination            : " << "\t" << m_Deklination << std::endl;
	std::cout << "m_Sonnen_Longitude       : " << "\t" << m_Sonnen_Longitude << std::endl;
	// Zusatzinfo
	std::cout << "Wellenlänge des Übergangs: " << "\t" << m_Wellenlaenge << std::endl;
	// Datum
	std::cout << "m_Jahr: " << "\t" << m_Jahr << std::endl;
	std::cout << "m_Monat: " << "\t" << m_Monat << std::endl;
	std::cout << "m_Tag: " << "\t" << m_Tag << std::endl;
	std::cout << "m_Stunde: " << "\t" << m_Stunde << std::endl;
	std::cout << "m_Minute: " << "\t" << m_Minute << std::endl;
	// Geolokationen
	std::cout << "m_Latitude_Sat: " << "\t" << m_Latitude_Sat << std::endl;
	std::cout << "m_Longitude_Sat: " << "\t" << m_Longitude_Sat << std::endl;
	std::cout << "m_Hoehe_Sat: " << "\t" << m_Hoehe_Sat << std::endl;
	std::cout << "m_Latitude_Ground: " << "\t" << m_Latitude_Ground << std::endl;
	std::cout << "m_Longitude_Ground: " << "\t" << m_Longitude_Ground << std::endl;
	std::cout << "m_Erdradius: " << "\t" << m_Erdradius << std::endl;
}
#endif /* AUSGEWERTETE_MESSUNG_NADIR_HH_ */
