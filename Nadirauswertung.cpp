/*
 * Nadirauswertung.cpp
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
#include "Glaetten.h"
#include "Konfiguration.h"

using std::string;
using std::vector;

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
					 vector<Ausgewertete_Messung_Nadir>& Ausgewertete_Nadirmessung_NO,
					 Konfiguration &Konf)
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

				mnit->Intensitaeten_durch_piF_Gamma_berechnen((*sfit), ldit->m_Gamma);
				// In der Formel ist piF in W/(m^2*Wellenlänge) verlangt..
				// also muss noch mit der Kanalbreite multipliziert werden
				mnit->Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen((*sfit));

				// Jetzt Zeilendichte und Fehler bestimmen
				mnit->Zeilendichte_Bestimmen((*sfit), k,
						Arbeitsverzeichnis, mache_Fit_Plots);

				// Ergebnis zusammenfassen
				Ausgewertete_Messung_Nadir Ergebnis = mnit->Ergebnis_Zusammenfassen();
				// Die braucht man später für die Luftmassenmatrix
				Ergebnis.m_Wellenlaenge
					= Ergebnis.m_Wellenlaenge_abs
					= ldit->m_Wellenlaenge;
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
				if (sfit->m_Spezies_Name == "NO") {
					double sol_fac = spidr_value_from_file(mnit->m_Jahr,
							mnit->m_Monat, mnit->m_Tag,
							Konf.m_Pfad_Solar_Correction_Factors, 1.0);
					std::cerr << "# solar factor = " << sol_fac << std::endl;
					// create new object, same transition but modelled temperature
					// fixed for now at 200 K since we are not sure which altitude
					// we would have to use for it.
					double temp = 200.0;
					int vu = sfit->NO_vec.at(k).get_vu();
					int vl = sfit->NO_vec.at(k).get_vl();
					int vl_abs = sfit->NO_vec.at(k).get_vl_abs();
					std::cerr << "# atmo temperature = " << temp
						<< ", vl_abs (v) = " << vl_abs
						<< ", vu (v') = " << vu
						<< ", vl (v'') = " << vl
						<< std::endl;
					NO_emiss NO_new(vu, vl, vl_abs, temp);
					NO_new.solar = sfit->NO_vec.at(k).solar * sol_fac;
					NO_new.read_luque_data_from_file(Konf.m_Pfad_NO_parameters);
					NO_new.calc_excitation();
					NO_new.calc_line_emissivities();
					//NO_new.pol_corr(mlit->m_TP_SZA, mlit->m_TP_rel_SAA, 0.17, -0.2);
					NO_new.scia_convolve(Rohdaten.at(0));
					double wl_abs = NO_new.get_wl_abs_vu_0();
					double wl_emiss = NO_new.get_wl_emiss_vu_vl();
					std::cerr << "# wls: abs = " << wl_abs << ", emiss = "
						<< wl_emiss << std::endl;
					mnit->slant_column_NO(NO_new, mache_Fit_Plots, Solspec, k,
							*sfit, Arbeitsverzeichnis, Konf);
					Ergebnis = mnit->Ergebnis_Zusammenfassen();
					Ergebnis.m_Wellenlaenge
						= ldit->m_Wellenlaenge
						= sfit->m_Wellenlaengen.at(k)
						= wl_emiss;
					Ergebnis.m_Wellenlaenge_abs = wl_abs;
					Ausgewertete_Nadirmessung_NO.push_back(Ergebnis);
				}
			}//ende k Linie
		}//ende j Spezies_Fenster
	}//Ende i Rohdaten

	return 0;
}//Ende Nadir_Auswertung
