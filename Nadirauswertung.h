/*
 * Nadirauswertung.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 30.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef NADIRAUSWERTUNG_HH_
#define NADIRAUSWERTUNG_HH_

// Diese Funktion dient lediglich dazu den Programmcode im Hauptprogramm kurz
// zu halten wichtig ist dabei, dass die nötigen Argumente übergeben werden und
// auch die Rückgabe komplett ist

#include <vector>

//Argumente von vorne nach hinten nach erster Nutzung sortiert
int Nadir_Auswertung(class Orbitliste &Orbitlist,
					 int l,
					 class Sonnenspektrum &Solspec,
					 std::vector<class Speziesfenster> &Spezies_Fenster,
					 int &counter_Nachtmessungen_Nadir,
					 int &counter_Nadir_Nacht_Dateien,
					 std::string Arbeitsverzeichnis, std::string mache_Fit_Plots,
					 std::vector<class Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_MgI,
					 std::vector<class Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_MgII,
					 std::vector<class Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_unknown,
					 std::vector<class Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_FeI,
					 std::vector<class Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_NO,
					 class Konfiguration &Konf);

#endif /* NADIRAUSWERTUNG_HH_ */
