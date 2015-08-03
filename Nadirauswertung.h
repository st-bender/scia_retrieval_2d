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
