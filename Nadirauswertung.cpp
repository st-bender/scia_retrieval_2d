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

int Nadir_Auswertung(Orbitliste &Orbitlist,
					 int l,
					 Sonnenspektrum &Solspec,
					 vector<Speziesfenster>& Spezies_Fenster,
					 int &counter_Nachtmessungen_Nadir,
					 int &counter_Nadir_Nacht_Dateien,
					 string Arbeitsverzeichnis, string mache_Fit_Plots,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_MgI,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_MgII,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_unknown,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_FeI,
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_NO)
{
	//cout<<"Start_Nadirauswertung\n";
	unsigned int i, j, k;
	vector<Messung_Nadir> Rohdaten;
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
		/////////////////////////////////////////////
		return 1;  //Nachtmessung 1
	}

	//cout<<Anzahl_Messungen<<"\n";
	vector<Messung_Nadir>::iterator mnit;
	vector<Speziesfenster>::iterator sfit;
	vector<Liniendaten>::iterator ldit;
	//Schleife über alle Rohdaten
	for (i = 0, mnit = Rohdaten.begin(); mnit != Rohdaten.end(); i++, ++mnit) {
		mnit->Deklinationswinkel_bestimmen();
		mnit->Sonnen_Longitude_bestimmen();
		mnit->Intensitaeten_normieren(Solspec.m_Int_interpoliert);

		//Schleife über alle Spezies wie z.b. Mg oder Mg+
		for (j = 0, sfit = Spezies_Fenster.begin();
				sfit != Spezies_Fenster.end(); j++, ++sfit) {

			//Schleife über alle Linien dieser Spezies
			for (k = 0, ldit = sfit->m_Liniendaten.begin();
					ldit != sfit->m_Liniendaten.end(); k++, ++ldit) {
				//Schleife über alle Linien dieser Spezies
				//Streuwinkel schon beim einlesen bestimmt
				//Spezfenst.m_Liniendaten[k].m_theta=Messung.m_Streuwinkel;
				//Streuwinkel muss woanders ermittelt werden
				ldit->Emissivitaet_ermitteln();
				//Spezfenst.m_Liniendaten[k].Auf_Bildschirm_Ausgeben();

				mnit->Intensitaeten_durch_piF_Gamma_berechnen((*sfit), k);

				// Jetzt Zeilendichte und Fehler bestimmen
				mnit->Zeilendichte_Bestimmen((*sfit), k,
						Arbeitsverzeichnis, mache_Fit_Plots, i);

				// Zu Testzwecken fertige Messung in Datei Speichern
				if ((k == 0) && (j == 0) && (i == 0)) {
					mnit->Ausgabe_in_Datei("CHECKDATA/Messung_Nadir_Fenster0_Hoehe_74km_0teLinie.txt");
				}
				// Ergebnis zusammenfassen
				Ausgewertete_Messung_Nadir Ergebnis = mnit->Ergebnis_Zusammenfassen();
				// Die braucht man später für die Luftmassenmatrix
				Ergebnis.m_Wellenlaenge = sfit->m_Wellenlaengen[k];
				//Ergebnis.Ausgabe_auf_Bildschirm();
				// Zusammenfassung der Zwischenresultate dem Vektor
				// für die jeweilige Spezies zuordnen
				//cout<<Spezies_Fenster[j].m_Spezies_Name<<"\n";
				if (sfit->m_Spezies_Name == "MgI") {
					//cout<<"Ausgewertete_Nadirmessung_MgI.push_back(Ergebnis)\n";
					Ausgewertete_Nadirmessung_MgI.push_back(Ergebnis);
				}
				if (sfit->m_Spezies_Name == "MgII") {
					Ausgewertete_Nadirmessung_MgII.push_back(Ergebnis);
				}
				if (sfit->m_Spezies_Name == "unknown") {
					Ausgewertete_Nadirmessung_unknown.push_back(Ergebnis);
				}
				if (sfit->m_Spezies_Name == "FeI") {
					Ausgewertete_Nadirmessung_FeI.push_back(Ergebnis);
				}
			}//ende k Linie
		}//ende j Spezies_Fenster
	}//Ende i Rohdaten

	return 0;
}//Ende Nadir_Auswertung
