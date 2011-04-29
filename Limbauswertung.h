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

#include <vector>

//Argumente von vorne nach hinten nach erster Nutzung sortiert
int Limb_Auswertung(class Orbitliste &Orbitlist,
					int l,
					class Sonnenspektrum &Solspec,
					std::vector<class Speziesfenster> &Spezies_Fenster,
					class NO_emiss &NO,
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
					Konfiguration &Konf);


#endif /* LIMBAUSWERTUNG_HH_ */
