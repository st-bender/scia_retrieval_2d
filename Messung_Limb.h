/*
 * Messung_Limb.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 12.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
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
	void moving_average(int window_size);
	void savitzky_golay(int window_size);
	double msise_temperature(class Konfiguration &Konf);

	Ausgewertete_Messung_Limb Ergebnis_Zusammenfassen();

	void Ausgabe_in_Datei(std::string Dateiname);
};

#endif /* MESSUNG_LIMB_HH_ */
