/*
 * Nadirauswertung.h
 *
 *  Created on: 30.04.2010
 *      Author: martin
 */

#ifndef NADIRAUSWERTUNG_HH_
#define NADIRAUSWERTUNG_HH_

// Diese Funktion dient lediglich dazu den Programmcode im Hauptprogramm kurz zu halten
// wichtig ist dabei, dass die nötigen Argumente übergeben werden und auch die Rückgabe komplett ist

#include<vector>
#include "Orbitliste.h"
#include "Sonnenspektrum.h"
#include "Speziesfenster.h"
#include "Ausgewertete_Messung_Nadir.h"


using namespace std;

//Argumente von vorne nach hinten nach erster Nutzung sortiert
int Nadir_Auswertung(Orbitliste Orbitlist,
                                  int l,
                                  Sonnenspektrum Solspec,
                                  vector<Speziesfenster>& Spezies_Fenster,
                                  int* counter_Nachtmessungen_Nadir,
                                  int* counter_Nadir_Nacht_Dateien,
                                  string Arbeitsverzeichnis, string mache_Fit_Plots,
                                  vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_MgI,
                                  vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_MgII,
                                  vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_unknown,
                                  vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_FeI);

#endif /* NADIRAUSWERTUNG_HH_ */
