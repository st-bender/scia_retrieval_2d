/*
 * Messung_Limb.h
 *
 *  Created on: 12.04.2010
 *      Author: martin langowski
 */
#ifndef MESSUNG_LIMB_HH_
#define MESSUNG_LIMB_HH_

#include <vector>
#include <string>
#include "Messung.h"
#include "Ausgewertete_Messung_Limb.h"

class Messung_Limb : public Messung
{
public:
	explicit Messung_Limb(std::string filename = "dummy");
	//Methoden
	void Zeilendichte_Bestimmen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots);
	void Saeulendichte_Bestimmen_MgI285nm(class Speziesfenster &Spezfenst,
			int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots,
			double *mean_10_20);
	void Plots_der_Spektren_erzeugen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots,
			double *mean_10_20);
	void slant_column_NO(class NO_emiss &NO, std::string mache_Fit_Plots,
			class Sonnenspektrum &sol_spec, int index,
			class Speziesfenster &Spezfenst, std::string Arbeitsverzeichnis,
			bool debug = true);
	void moving_average(int window_size);
	void savitzky_golay(int window_size);
	double msise_temperature(class Konfiguration &Konf);

	Ausgewertete_Messung_Limb Ergebnis_Zusammenfassen();

	void Ausgabe_in_Datei(std::string Dateiname);
};

#endif /* MESSUNG_LIMB_HH_ */
