/*
 * Messung_Nadir.h
 *
 *  Created on: 28.04.2010
 *      Author: martin
 */

#ifndef MESSUNG_NADIR_HH_
#define MESSUNG_NADIR_HH_

#include <vector>
#include <string>
#include "Messung.h"
#include "Ausgewertete_Messung_Nadir.h"

class Messung_Nadir : public Messung
{
	/****************************************************
	 Der einzige große Unterschied in den Funktionen
	 zwischen Nadir und Limb besteht hier darin,
	 dass ich bei Limb noch Vektors verwendet hab
	 ****************************************************/
public:
	explicit Messung_Nadir(std::string filename = "dummy");
	// Assignmentoperator Overload
	Messung_Nadir &operator =(const Messung_Nadir &rhs);
	//Methoden
	void Zeilendichte_Bestimmen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots, int MessungsNr);
	Ausgewertete_Messung_Nadir Ergebnis_Zusammenfassen();

	//Wartungsfunktionen
	//zum Testen und debuggen und überprüfen, ob der fit halbwegs passt
	void Ausgabe_in_Datei(std::string Dateiname);

	// neue Membervariablen

	// Herkunftsmerkmale
	int m_Messung_ID;
	// Geolokationen für Raytrace
	double m_Latitude_Ground;
	double m_Longitude_Ground;
};

#endif /* MESSUNG_NADIR_HH_ */
