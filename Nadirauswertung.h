/*
 * Nadirauswertung.h
 *
 *  Created on: 30.04.2010
 *      Author: martin
 */

#ifndef NADIRAUSWERTUNG_HH_
#define NADIRAUSWERTUNG_HH_

// Diese Funktion dient lediglich dazu den Programmcode im Hauptprogramm kurz
// zu halten wichtig ist dabei, dass die nötigen Argumente übergeben werden und
// auch die Rückgabe komplett ist

#include<vector>
#include "Orbitliste.h"
#include "Sonnenspektrum.h"
#include "Speziesfenster.h"
#include "Ausgewertete_Messung_Nadir.h"

//Argumente von vorne nach hinten nach erster Nutzung sortiert
int Nadir_Auswertung(Orbitliste Orbitlist,
					 int l,
					 Sonnenspektrum &Solspec,
					 std::vector<Speziesfenster>& Spezies_Fenster,
					 int &counter_Nachtmessungen_Nadir,
					 int &counter_Nadir_Nacht_Dateien,
					 std::string Arbeitsverzeichnis, std::string mache_Fit_Plots,
					 std::vector<Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_MgI,
					 std::vector<Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_MgII,
					 std::vector<Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_unknown,
					 std::vector<Ausgewertete_Messung_Nadir> &Ausgewertete_Nadirmessung_FeI);

#endif /* NADIRAUSWERTUNG_HH_ */
