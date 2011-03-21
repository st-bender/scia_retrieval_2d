/*
 * Nadirauswertung.cpp
 *
 *  Created on: 30.04.2010
 *      Author: martin
 */

#include<vector>
#include "Orbitliste.h"
#include "Sonnenspektrum.h"
#include "Speziesfenster.h"
#include "Ausgewertete_Messung_Nadir.h"
// TODO Das kann weg, komisch, dass das ohne Messung_Nadir funktioniert hat
//#include "Messung_Limb.h"
#include "Messung_Nadir.h"
#include "Datei_IO.h"  //ReadL1C_Nadir_mpl_binary
#include"Messungs_ausschliessende_Tests.h"

using namespace std;

int Nadir_Auswertung(Orbitliste Orbitlist,
					 int l,
					 Sonnenspektrum Solspec,
					 vector<Speziesfenster>& Spezies_Fenster,
					 int &counter_Nachtmessungen_Nadir,
					 int &counter_Nadir_Nacht_Dateien,
					 string Arbeitsverzeichnis, string mache_Fit_Plots,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_MgI,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_MgII,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_unknown,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_FeI)
{
	//cout<<"Start_Nadirauswertung\n";
	unsigned int j, k;
	Messung_Nadir *Rohdaten = 0;
	int Anzahl_Messungen = 0;
	Rohdaten = ReadL1C_Nadir_mpl_binary(Orbitlist.m_Dateinamen[l], Anzahl_Messungen);
	//cerr<<Orbitlist.m_Dateinamen[l]<<" wird bearbeitet\n";
	//cout<<Rohdaten[0].m_Wellenlaengen[0]<<"\t"
	//  <<Rohdaten[0].m_Intensitaeten[0]<<"\n";

	// Hier erste Messungen aussortieren
	bool ist_Nachtmessung;
	ist_Nachtmessung = Test_auf_Nachtmessung_Nadir(Rohdaten, Anzahl_Messungen);
	if (ist_Nachtmessung == true) {
		counter_Nachtmessungen_Nadir += Anzahl_Messungen; //counter setzen
		counter_Nadir_Nacht_Dateien++;
		// Speicher freimachen /////////////
		for (int i = 0; i < Anzahl_Messungen; i++) {
			// das muss sein..unschön..evtl ein Objekt entwerfen, was nur aus
			// so einem Array besteht dann kann man das über den Destruktor
			// erledigen
			Rohdaten[i].save_delete_all_memory();
		}
		if (Rohdaten != 0) {
			delete[] Rohdaten;
		}
		/////////////////////////////////////////////
		return 1;  //Nachtmessung 1
	}


	//cout<<Anzahl_Messungen<<"\n";
	for (int i = 0; i < Anzahl_Messungen; i++) { //Schleife über alle Rohdaten
		Messung_Nadir Messung = Rohdaten[i];
		Messung.Deklinationswinkel_bestimmen();
		Messung.Sonnen_Longitude_bestimmen();
		Messung.Intensitaeten_normieren(Solspec.m_Int_interpoliert);

		for (j = 0; j < Spezies_Fenster.size(); j++) {
			//Schleife über alle Spezies wie z.b. Mg oder Mg+

			//Speziesfenster  Spezfenst=Spezies_Fenster[j];
			for (k = 0; k < Spezies_Fenster[j].m_Wellenlaengen.size(); k++) {
				//Schleife über alle Linien dieser Spezies
				//Streuwinkel schon beim einlesen bestimmt
				//Spezfenst.m_Liniendaten[k].m_theta=Messung.m_Streuwinkel;
				//Streuwinkel muss woanders ermittelt werden
				Spezies_Fenster[j].m_Liniendaten[k].Emissivitaet_ermitteln();
				//Spezfenst.m_Liniendaten[k].Auf_Bildschirm_Ausgeben();

				Messung.Intensitaeten_durch_piF_Gamma_berechnen(Spezies_Fenster[j], k);

				// Jetzt Zeilendichte und Fehler bestimmen
				Messung.Zeilendichte_Bestimmen(Spezies_Fenster[j], k,
						Arbeitsverzeichnis, mache_Fit_Plots, i);

				// Zu Testzwecken fertige Messung in Datei Speichern
				if ((k == 0) && (j == 0) && (i == 0)) {
					Messung.Ausgabe_in_Datei("CHECKDATA/Messung_Nadir_Fenster0_Hoehe_74km_0teLinie.txt");
				}
				// Ergebnis zusammenfassen
				Ausgewertete_Messung_Nadir Ergebnis = Messung.Ergebnis_Zusammenfassen();
				// Die braucht man später für die Luftmassenmatrix
				Ergebnis.m_Wellenlaenge = Spezies_Fenster[j].m_Wellenlaengen[k];
				//Ergebnis.Ausgabe_auf_Bildschirm();
				// Zusammenfassung der Zwischenresultate dem Vektor
				// für die jeweilige Spezies zuordnen
				//cout<<Spezies_Fenster[j].m_Spezies_Name<<"\n";
				if (Spezies_Fenster[j].m_Spezies_Name == "MgI") {
					//cout<<"Ausgewertete_Nadirmessung_MgI.push_back(Ergebnis)\n";
					Ausgewertete_Nadirmessung_MgI.push_back(Ergebnis);
				}
				if (Spezies_Fenster[j].m_Spezies_Name == "MgII") {
					Ausgewertete_Nadirmessung_MgII.push_back(Ergebnis);
				}
				if (Spezies_Fenster[j].m_Spezies_Name == "unknown") {
					Ausgewertete_Nadirmessung_unknown.push_back(Ergebnis);
				}
				if (Spezies_Fenster[j].m_Spezies_Name == "FeI") {
					Ausgewertete_Nadirmessung_FeI.push_back(Ergebnis);
				}
			}//ende k Linie
		}//ende j Spezies_Fenster
		//Speicher der Messung löschen
		// das muss eigentlich nicht hier stehn,
		// da Messung gleich seinen Destruktor aufruft
		Messung.save_delete_all_memory();
	}//Ende i Rohdaten

	for (int i = 0; i < Anzahl_Messungen; i++) {
		// das muss sein..unschön..
		// evtl ein Objekt entwerfen, was nur aus so einem Array besteht
		// dann kann man das über den Destruktor erledigen
		Rohdaten[i].save_delete_all_memory();
	}
	if (Rohdaten != 0) {
		delete[] Rohdaten;
	}

	return 0;
}//Ende Nadir_Auswertung
