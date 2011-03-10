/*
 * Limbauswertung.h
 *
 *  Created on: 28.04.2010
 *      Author: martin
 */
#ifndef LIMBAUSWERTUNG_HH_
#define LIMBAUSWERTUNG_HH_

// Diese Funktion dient lediglich dazu den Programmcode im Hauptprogramm kurz zu
// halten wichtig ist dabei, dass die nötigen Argumente übergeben werden und
// auch die Rückgabe komplett ist

#include<vector>
#include "Orbitliste.h"
#include "Sonnenspektrum.h"
#include "Speziesfenster.h"
#include "Ausgewertete_Messung_Limb.h"


using namespace std;

//Argumente von vorne nach hinten nach erster Nutzung sortiert
int Limb_Auswertung(Orbitliste Orbitlist,
					int l,
					Sonnenspektrum &Solspec,
					vector<Speziesfenster>& Spezies_Fenster,
					int *counter_Nachtmessungen,
					int *counter_NLC_detektiert,
					int *counter_Richtungsvektor_nicht_ok,
					string Arbeitsverzeichnis,
					string mache_Fit_Plots,
					string limb_meso_thermo,       // "ja" oder "nein"
					vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_MgI,
					vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_MgII,
					vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_unknown,
					vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_FeI);


#endif /* LIMBAUSWERTUNG_HH_ */
