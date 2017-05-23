/*
 * Messung_Nadir.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 28.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
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
	//Methoden
	void Zeilendichte_Bestimmen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots);
	Ausgewertete_Messung_Nadir Ergebnis_Zusammenfassen();

	//Wartungsfunktionen
	//zum Testen und debuggen und überprüfen, ob der fit halbwegs passt
	void Ausgabe_in_Datei(std::string Dateiname);
};

#endif /* MESSUNG_NADIR_HH_ */
