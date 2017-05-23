/*
 * Limbauswertung.h
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
#ifndef LIMBAUSWERTUNG_HH_
#define LIMBAUSWERTUNG_HH_

// Diese Funktion dient lediglich dazu den Programmcode im Hauptprogramm kurz zu
// halten wichtig ist dabei, dass die nötigen Argumente übergeben werden und
// auch die Rückgabe komplett ist

#include <vector>

//Argumente von vorne nach hinten nach erster Nutzung sortiert
int Limb_Auswertung(class Orbitliste &Orbitlist,
					int l,
					class Sonnenspektrum &Solspec,
					class Sonnenspektrum &sol_ref,
					std::vector<class Speziesfenster> &Spezies_Fenster,
					int &counter_Nachtmessungen,
					int &counter_NLC_detektiert,
					int &counter_Richtungsvektor_nicht_ok,
					std::string Arbeitsverzeichnis,
					std::string mache_Fit_Plots,
					std::string limb_meso_thermo,       // "ja" oder "nein"
					std::vector<class Ausgewertete_Messung_Limb> &Ausgewertete_Limbmessung_MgI,
					std::vector<class Ausgewertete_Messung_Limb> &Ausgewertete_Limbmessung_MgII,
					std::vector<class Ausgewertete_Messung_Limb> &Ausgewertete_Limbmessung_unknown,
					std::vector<class Ausgewertete_Messung_Limb> &Ausgewertete_Limbmessung_FeI,
					std::vector<class Ausgewertete_Messung_Limb> &Ausgewertete_Limbmessung_NO,
					Konfiguration &Konf);


#endif /* LIMBAUSWERTUNG_HH_ */
