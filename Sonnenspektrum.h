/*
 * Sonnenspektrum.h
 *
 *  Created on: 20.04.2010
 *      Author: martin
 */

#ifndef SONNENSPEKTRUM_HH_
#define SONNENSPEKTRUM_HH_

#include <string>
#include "Messung_Limb.h"

using namespace std;

class Sonnenspektrum
{
public:
	// Diese Funktion ist für die Sonnenspektren von GOME von
	// Mark Weber
	int Laden_GOME(string Dateiname, string Fallback_Dateiname);
	//Sciamachy Sonnenspektrum des Orbits
	int Laden_SCIA(string Dateiname, string Fallback_Dateiname);

	int moving_average(int window_size);
	int savitzky_golay(int window_size);
	int saoref_to_schiamachy();
	int Interpolieren(Messung_Limb &Messung_Erdschein);
	int nicht_interpolieren();
	vector<double> m_Wellenlaengen; // wie lang sind die eigentlich
	vector<double> m_Intensitaeten;
	// Auf Erdschein-Messungs-Wellenlängen Interpolierte Intensität
	// (Wellenlängen sind dann gleich dem ErdscheinWL)
	vector<double> m_Int_interpoliert;
	vector<double> m_WL_interpoliert; //für Ausgabe
	int m_Anzahl_WL;
	//Wartungsfunktion
	int Speichern(string Dateiname);  //zur Kontrolle
	int Speichern_was_geladen_wurde(string Dateiname);
};


#endif /* SONNENSPEKTRUM_HH_ */
