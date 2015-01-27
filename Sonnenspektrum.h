/*
 * Sonnenspektrum.h
 *
 *  Created on: 20.04.2010
 *      Author: martin
 */

#ifndef SONNENSPEKTRUM_HH_
#define SONNENSPEKTRUM_HH_

#include <string>
#include <vector>

class Sonnenspektrum
{
public:
	// Diese Funktion ist für die Sonnenspektren von GOME von
	// Mark Weber
	int Laden_GOME(std::string Dateiname, std::string Fallback_Dateiname);
	//Sciamachy Sonnenspektrum des Orbits
	int Laden_SCIA(std::string Dateiname, std::string Fallback_Dateiname);

	void moving_average(int window_size);
	void savitzky_golay(int window_size);
	void saoref_to_sciamachy();
	double get_rad_at_wl(double wl);
	void Interpolieren(class Messung_Limb &Messung_Erdschein);
	void nicht_interpolieren();
	std::vector<double> m_Wellenlaengen; // wie lang sind die eigentlich
	std::vector<double> m_Intensitaeten;
	// Auf Erdschein-Messungs-Wellenlängen Interpolierte Intensität
	// (Wellenlängen sind dann gleich dem ErdscheinWL)
	std::vector<double> m_Int_interpoliert;
	std::vector<double> m_WL_interpoliert; //für Ausgabe
	int m_Anzahl_WL;
	double int_ref1, int_ref2;
	//Wartungsfunktion
	int Speichern(std::string Dateiname);  //zur Kontrolle
	int Speichern_was_geladen_wurde(std::string Dateiname);
};


#endif /* SONNENSPEKTRUM_HH_ */
