/*
 * Datei_IO.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 19.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include <cmath>
#include "Messung_Limb.h"
#include "Messung_Nadir.h"
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream>
#include <cstdlib> //system
#include "Nachricht_Schreiben.h"
#include "MPL_Matrix.h"
#include "LimbNadir_IO.h"
#include "Datei_IO.h"
#include "Retrievalgitter.h"
#include "Glaetten.h"

using std::string;
using std::vector;

extern int Prioritylevel;

// calculates the intensity average over a range of wavelengths.
// when median is set to true, it returns the median in the range,
// otherwise the arithmetic mean is returned (the default case)
double average_over_wl_range(std::vector<float> rad, std::vector<float> wl,
		double wl_start, double wl_end, bool median = false)
{
	double avg = 0.;
	// copy to a vector for <algorithm>
	vector<float>::iterator wl_low, wl_up;

	if (wl_start > wl_end) {
		// swap
		double tmp = wl_start;
		wl_start = wl_end;
		wl_end = tmp;
	}

	wl_low = lower_bound(wl.begin(), wl.end(), wl_start);
	if (wl_low == wl.end()) return 0.; // not in range

	wl_up = upper_bound(wl.begin(), wl.end(), wl_end);
	if (wl_up != wl.begin()) --wl_up;

	if (wl_up == wl_low) return 0.; // nothing to average

	if (median == false) {
		// build the arithmetic mean
		avg = accumulate(rad.begin(), rad.end(), 0.);
		avg /= rad.size();
	} else {
		// find the median
		nth_element(rad.begin(), rad.begin() + rad.size() / 2, rad.end());
		avg = rad.at(rad.size() / 2);
	}

	return avg;
}

// helper function to copy Limb_Datensatz *Limbdaten and
// float *Wellenlaengen into a vector<Messung_Limb>
vector<Messung_Limb> make_messung_limb_vector(string Dateiname,
		std::vector<Limb_Datensatz> &Limbdaten, std::vector<float> &Wellenlaengen,
		int no_of_alt, int no_of_pix, int Datum[6], float cent_lat_lon[10],
		float orbit_phase, int no_of_heights, int offset, int direction,
		double dark_sig)
{
	bool has_straylight = false;
	// dark signal and error
	// constant dark signal (default and fall-back)
	//dark_sig = 2.731e9;
	//double dark_sig = 3.9e9;
	//double dark_sig = 0.0;
	double dark_err = 0.0;

	if (dark_sig == -1.) {
		// normal average or median for the dark signal correction
		dark_sig = average_over_wl_range(Limbdaten[no_of_alt - 1].m_radiance,
				Wellenlaengen, 238.0, 282.0, true);
		dark_err = average_over_wl_range(Limbdaten[no_of_alt - 1].m_error,
				Wellenlaengen, 238.0, 282.0, true);

		if (dark_sig > 6.e9)
			has_straylight = true;
	}
	std::cerr << "dark signal: " << dark_sig << std::endl;

	// 4. Erstellung des Übergabevektors
	vector<Messung_Limb> Ergebnisvektor;

	for (int i = 0; i < no_of_heights; i++) {
		Messung_Limb ml(Dateiname);
		ml.m_Jahr = Datum[0];
		ml.m_Monat = Datum[1];
		ml.m_Tag = Datum[2];
		ml.m_Stunde = Datum[3];
		ml.m_Minute = Datum[4];
		ml.m_Sekunde = Datum[5];
		ml.m_orbit_phase = orbit_phase;
		ml.center_lat = cent_lat_lon[0];
		ml.center_lon = cent_lat_lon[1];
		ml.m_Latitude_Sat = Limbdaten[offset + direction * i].m_Sub_Sat_Lat; //achtung geodätische Koordinaten
		ml.m_Longitude_Sat = Limbdaten[offset + direction * i].m_Sub_Sat_Lon;
		ml.m_Hoehe_Sat = Limbdaten[offset + direction * i].m_Sat_Hoehe;
		ml.m_Latitude_TP = Limbdaten[offset + direction * i].m_TP_Lat;
		ml.m_Longitude_TP = Limbdaten[offset + direction * i].m_TP_Lon;
		ml.m_Hoehe_TP = Limbdaten[offset + direction * i].m_Tangentenhoehe;
		ml.m_Erdradius = Limbdaten[offset + direction * i].m_Erdradius;
		ml.m_TP_SZA = Limbdaten[offset + direction * i].m_TP_SZA;
		ml.m_TP_rel_SAA = Limbdaten[offset + direction * i].m_TP_SAA;
		ml.m_Number_of_Wavelength = no_of_pix;

		int idx_282nm = 0;
		for (int j = 0; j < no_of_pix; j++) {
			ml.m_Wellenlaengen.push_back(Wellenlaengen[j]);
			if (has_straylight) {
				// the old corrections
				dark_sig = Limbdaten[no_of_alt - 1].m_radiance[j];
				dark_err = Limbdaten[no_of_alt - 1].m_error[j];
			}
			double signal = Limbdaten[offset + direction * i].m_radiance[j];
			double relerr = Limbdaten[offset + direction * i].m_error[j];
			ml.m_Intensitaeten.push_back(signal - dark_sig);
			ml.m_Intensitaeten_relativer_Fehler.push_back(
					std::sqrt(((relerr * signal) * (relerr * signal) +
					(dark_err * dark_sig) * (dark_err * dark_sig))) /
					ml.m_Intensitaeten.back());
			ml.m_Sonne.push_back(0.);
			ml.m_Intensitaeten_durch_piF.push_back(0.);
			ml.m_Intensitaeten_durch_piF_Gamma.push_back(0.);
			ml.m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand.push_back(0.);
			// save the index for 282.03 nm for later
			if (Wellenlaengen[j] < 282.04) idx_282nm = j;
		}

		// Der Pixel 552 (282,03nm zeigt bei nadir( und nur dort, einen Peak)
		// ....interpretation dead pixel
		// use the saved index for 282.03 nm and interpolate the value
		ml.m_Intensitaeten.at(idx_282nm) = 0.5 *
			(ml.m_Intensitaeten.at(idx_282nm - 1) +
			 ml.m_Intensitaeten.at(idx_282nm + 1));
		Ergebnisvektor.push_back(ml);
	}

	return Ergebnisvektor;
}

////////////////////////////////////////////////////////////////////////////////
//
// Funktionsstart ReadL1C_Limb_mpl_binary
//
////////////////////////////////////////////////////////////////////////////////
vector<Messung_Limb> ReadL1C_Limb_mpl_binary(string Dateiname,
		Messung_Limb &Troposphaerische_Saeule, Messung_Limb &mean_10_20,
		int Anzahl_Hoehen, double dark_bg)
{
	//binärdateien sind nicht gepackt(das wär einfach nicht effizient)...
	//ansonsten hier packen und später entpacken
	//..siehe alte versionen
	// erst wird alles geladen, dann analyse und nachberarbeitung durchgeführt

	// 1. Zur Verfügung Stellung der Speicherstrukturen zur Aufnahme der Datei
	// 2. Laden der Datei
	// 3. Nachbearbeitung/Ausschlusskriterien
	//    -> auf seperate Funktion nach laden verschoben worden
	// 4. Erstellung des Übergabevektors
	// //ACHTUNG für Troposphaerische_Saeule nur Intensitäten
	// 5. Speicherfreigabe
	// 6. Rückgabe

	// 1. Zur Verfügung Stellung der Speicherstrukturen zur Aufnahme der Datei
	//int lang_textheader=31;
	string textheader[31];
	int no_of_alt = 0;
	int no_of_pix = 0;
	int Orbitstate[5];
	int Datum[6];
	float Center_Lat_Lon[10];
	float orbit_phase;
	std::vector<float> Wellenlaengen;
	std::vector<Limb_Datensatz> Limbdaten;
	std::vector<Messung_Limb> Ergebnisvektor;
	// 2. Laden der Datei
	//cerr<<" 2. Laden der Datei\n";
	int err = Load_Limb_l_mpl_binary(Dateiname,
						   textheader, no_of_alt, no_of_pix, Orbitstate, Datum,
						   Center_Lat_Lon, orbit_phase, Wellenlaengen,
						   Limbdaten);
	// return empty vector on error
	if (err) return Ergebnisvektor;
	// 3. Nachbearbeitung/Ausschlusskriterien
	// Wird jetzt nach dem Laden durchgeführt
	// 4. Erstellung des Übergabevektors
	// check Anzahl_Hoehen for correctness and change it to
	// the largest sensible number of altitudes
	if (Anzahl_Hoehen < 0 || Anzahl_Hoehen > no_of_alt - 1)
		Anzahl_Hoehen = no_of_alt - 1;
	Ergebnisvektor
		= make_messung_limb_vector(Dateiname, Limbdaten, Wellenlaengen,
				no_of_alt, no_of_pix, Datum, Center_Lat_Lon, orbit_phase,
				Anzahl_Hoehen, no_of_alt - Anzahl_Hoehen - 1, 1, dark_bg);

	//Teile von Schritt 4 nochmal für die Troposhärische Säule
	//Eigentlich reichen Intensitäten
	for (int j = 0; j < no_of_pix; j++) {
		Troposphaerische_Saeule.m_Intensitaeten.push_back(Limbdaten[1].m_radiance[j]);
	}
	Troposphaerische_Saeule.m_TP_SZA = Limbdaten[1].m_TP_SZA;
	// und für den Mittelwert, der Höhen 10 bis 20
	for (int j = 0; j < no_of_pix; j++) {
		double mean = 0.;
		for (int k = 10; k < 21; k++) {
			mean += Limbdaten[k].m_radiance[j];
		}
		mean /= 11.0;
		mean_10_20.m_Intensitaeten.push_back(mean);
	}

	return Ergebnisvektor;
}
////////////////////////////////////////////////////////////////////////////////
//
// ENDE ReadL1C_Limb_mpl_binary
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Funktionsstart ReadL1C_Limb_meso_thermo_mpl_binary
//
////////////////////////////////////////////////////////////////////////////////
// this is a superfluous function and should be deprecated
// since it might be used (somewhere) we keep it for now
vector<Messung_Limb> ReadL1C_Limb_meso_thermo_mpl_binary(string Dateiname,
		Messung_Limb &niedrigste_Hoehe, Messung_Limb &space)
{
	///////////////////////////////////////////////////////////
	// ähnlich zur üblichen Limbroutine...nur andere TH Reihenfolge und alle
	// Messungen über 70 km werden geladen und die niedrigste bei 53.5 km
	// 0 bis 24   z.B. bei 0 148,7km und bei 25 69.9 km
	///////////////////////////////////////////////////////////

	return ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(Dateiname,
			niedrigste_Hoehe, space, 25);
}
////////////////////////////////////////////////////////////////////////////////
//
// ENDE ReadL1C_Limb_meso_thermo_mpl_binary
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Funktionsstart ReadL1C_Limb_meso_thermo_mpl_binary_reduziert
//
////////////////////////////////////////////////////////////////////////////////
vector<Messung_Limb>
ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(string Dateiname,
		Messung_Limb &niedrigste_Hoehe, Messung_Limb &space, int Anzahl_Hoehen,
		double dark_bg)
{
	// Hier wieder nur Höhen von 70 bis 90 km....(einziger unterschied liegt in
	// der for schleife die nur bis 7 geht)
	///////////////////////////////////////////////////////////
	// ähnlich zur üblichen Limbroutine...nur andere TH Reihenfolge und alle
	// Messungen über 70 km werden geladen und die niedrigste bei 53.5 km
	// 0 bis 24   z.B. bei 0 148,7km und bei 25 69.9 km
	///////////////////////////////////////////////////////////

	//binärdateien sind nicht gepackt(das wär einfach nicht effizient)...
	//ansonsten hier packen und später entpacken
	//..siehe alte versionen
	// erst wird alles geladen, dann analyse und nachberarbeitung durchgeführt

	// 1. Zur Verfügung Stellung der Speicherstrukturen zur Aufnahme der Datei
	// 2. Laden der Datei
	// 3. Nachbearbeitung/Ausschlusskriterien
	//    -> auf seperate Funktion nach laden verschoben worden
	// 4. Erstellung des Übergabevektors
	// //ACHTUNG für Troposphaerische_Saeule nur Intensitäten
	// 5. Speicherfreigabe
	// 6. Rückgabe

	// 1. Zur Verfügung Stellung der Speicherstrukturen zur Aufnahme der Datei
	//int lang_textheader=31;
	string textheader[31];
	int no_of_alt = 0;
	int no_of_pix = 0;
	int Orbitstate[5];
	int Datum[6];
	float Center_Lat_Lon[10];
	float orbit_phase;
	std::vector<float> Wellenlaengen;
	std::vector<Limb_Datensatz> Limbdaten;
	std::vector<Messung_Limb> Ergebnisvektor;
	// 2. Laden der Datei
	//cerr<<" 2. Laden der Datei\n";
	int err = Load_Limb_l_mpl_binary(Dateiname,
						   textheader, no_of_alt, no_of_pix, Orbitstate, Datum,
						   Center_Lat_Lon, orbit_phase, Wellenlaengen,
						   Limbdaten);
	// return empty vector on error
	if (err) return Ergebnisvektor;
	// 3. Nachbearbeitung/Ausschlusskriterien
	// Wird jetzt nach dem Laden durchgeführt
	// 4. Erstellung des Übergabevektors
	// check Anzahl_Hoehen for correctness and change it to
	// the largest sensible number of altitudes
	if (Anzahl_Hoehen < 0 || Anzahl_Hoehen > no_of_alt - 1)
		Anzahl_Hoehen = no_of_alt - 1;
	Ergebnisvektor
		= make_messung_limb_vector(Dateiname, Limbdaten, Wellenlaengen,
				no_of_alt, no_of_pix, Datum, Center_Lat_Lon, orbit_phase,
				Anzahl_Hoehen, Anzahl_Hoehen - 1, -1, dark_bg);

	//Teile von Schritt 4 nochmal für die niedrigste Höhe
	//Eigentlich reichen Intensitäten
	niedrigste_Hoehe = Ergebnisvektor.front();

	// Prepares a one element vector with the last entry in the measurements
	// vector containing the "dark" scan at around 360 km.
	// It then returns the only element as the "space" limb scan as requested.
	space = make_messung_limb_vector(Dateiname, Limbdaten, Wellenlaengen,
				no_of_alt, no_of_pix, Datum, Center_Lat_Lon, orbit_phase,
				1, Anzahl_Hoehen, +1, dark_bg).front();

	return Ergebnisvektor;
}
////////////////////////////////////////////////////////////////////////////////
//
// ENDE ReadL1C_Limb_meso_thermo_mpl_binary_reduziert
//
////////////////////////////////////////////////////////////////////////////////


// helper function to copy Limb_Datensatz *Limbdaten and
// float *Wellenlaengen into a vector<Messung_Limb>
vector<Messung_Nadir> make_messung_nadir_vector(string Dateiname,
		std::vector<Nadir_Datensatz> &Nadirdaten,
		std::vector<float> &Wellenlaenge, int No_of_Messungen, int No_of_Pix)
{
	// 4. Erstellung des Übergabevektors
	vector<Messung_Nadir> mn_vec;

	for (int i = 0; i < No_of_Messungen; i++) {
		Messung_Nadir mn(Dateiname);
		//Herkunftsmerkmale
		mn.m_Messung_ID = Nadirdaten[i].m_Messung_ID;
		//Datum
		mn.m_Jahr = Nadirdaten[i].m_Jahr;
		mn.m_Monat = Nadirdaten[i].m_Monat;
		mn.m_Tag = Nadirdaten[i].m_Tag;
		mn.m_Stunde = Nadirdaten[i].m_Stunde;
		mn.m_Minute = Nadirdaten[i].m_Minute;
		mn.m_Sekunde = Nadirdaten[i].m_Sekunde;
		//Geolokationen
		mn.m_Latitude_Sat = Nadirdaten[i].m_Sat_Lat;
		mn.m_Longitude_Sat = Nadirdaten[i].m_Sat_Lon;
		mn.m_Hoehe_Sat = Nadirdaten[i].m_Hoehe;
		mn.m_Latitude_Ground = Nadirdaten[i].m_geo_nadir_center_lat;
		mn.m_Longitude_Ground = Nadirdaten[i].m_geo_nadir_center_lon;
		mn.m_Erdradius = Nadirdaten[i].m_Sat_Erdradius;
		mn.m_orbit_phase = Nadirdaten[i].m_orbit_phase;

		int idx_282nm = 0;
		//Füllbare Felder
		mn.m_Number_of_Wavelength = No_of_Pix;
		//Felder allokieren
		//Deep Copy der Wellenlängen und Intensitäten
		for (int j = 0; j < No_of_Pix; j++) {
			mn.m_Wellenlaengen.push_back(Wellenlaenge[j]);
			mn.m_Intensitaeten.push_back(Nadirdaten[i].m_radiance[j]);
			mn.m_Intensitaeten_relativer_Fehler.push_back(Nadirdaten[i].m_error[j]);
			mn.m_Intensitaeten_durch_piF.push_back(0.);
			mn.m_Intensitaeten_durch_piF_Gamma.push_back(0.);
			mn.m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand.push_back(0.);
			// save the index for 282.03 nm for later
			if (Wellenlaenge[j] < 282.04) idx_282nm = j;
		}
		// Der Pixel 552 (282,03nm zeigt bei nadir( und nur dort, einen Peak)
		// Das ist die Kanalgrenze zwischen Kanal1a und Kanal1b
		// Interpolieren zwischen 282 und 283 nm
		//mn.m_Intensitaeten[536]=(mn.m_Intensitaeten[535]+mn.m_Intensitaeten[537])/2;
		//mn.m_Intensitaeten_relativer_Fehler[536]=(mn.m_Intensitaeten_relativer_Fehler[535]+mn.m_Intensitaeten_relativer_Fehler[537])/2;
		mn.m_Intensitaeten.at(idx_282nm) = 0.5 *
			(mn.m_Intensitaeten.at(idx_282nm - 1) +
			 mn.m_Intensitaeten.at(idx_282nm + 1));

		mn_vec.push_back(mn);
	}

	return mn_vec;
}

////////////////////////////////////////////////////////////////////////////////
//
// Funktionsstart ReadL1C_Nadir_mpl_binary
//
////////////////////////////////////////////////////////////////////////////////
vector<Messung_Nadir> ReadL1C_Nadir_mpl_binary(string Dateiname, int &Anzahl_Messungen)
{
	// 1. Zur Verfügung Stellung der Speicherstrukturen zur Aufnahme der Datei
	// 2. Laden der Datei
	// 3. Nachbearbeitung/Ausschlusskriterien
	//    -> wird nach dieser Routine als eigenen Funktion implementiert
	// 4. Erstellung des Übergabevektors
	// 5. Speicherfreigabe
	// 6. Rückgabe

	// 1. Zur Verfügung Stellung der Speicherstrukturen zur Aufnahme der Datei
	string textheader[7];
	int No_of_Messungen;  // Das variiert bei nadir doch...standard ist 65
	int No_of_Pix;
	std::vector<int> Kanal_Nr;
	std::vector<float> Wellenlaenge;
	std::vector<Nadir_Datensatz> Nadirdaten;
	std::vector<Messung_Nadir> aus;
	// 2. Laden der Datei
	int err = Load_Nadir_n_mpl_binary(Dateiname,
							textheader, No_of_Messungen, No_of_Pix,
							Kanal_Nr, Wellenlaenge, Nadirdaten);
	// return empty vector on error
	if (err) {
		Anzahl_Messungen = 0;
		return aus;
	}
	// 3. Nachbearbeitung/Ausschlusskriterien
	// Dieser Schritt wird erst nach dem Laden aufgerufen
	// 4. Erstellung des Übergabefelds
	Anzahl_Messungen = No_of_Messungen;
	aus
		= make_messung_nadir_vector(Dateiname, Nadirdaten, Wellenlaenge,
				No_of_Messungen, No_of_Pix);

	// 6. Rückgabe
	return aus;
}
////////////////////////////////////////////////////////////////////////////////
//
// ENDE ReadL1C_Nadir_mpl_binary
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//START int Ausgabe_Zeilendichten_Limb(string Dateiname);
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////
void Ausgabe_Saeulendichten(string Dateiname,
		vector<Ausgewertete_Messung_Limb> &A_Messung_L)
{
	size_t lang = A_Messung_L.size();
	if (lang == 0)
		return;
	//Formatierte Ausgabe
	FILE *outfile;
	//Datei öffnen
	//cerr<<"Datei öffnen\n";
	//cerr<<"Dateiname.c_str(): "<<Dateiname.c_str()<<"\n";
	outfile = fopen(Dateiname.c_str(), "w");
	//Überschrift
	//cerr<<"Überschrift\n";
	fprintf(outfile, "%4s %5s %3s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
			"Jahr", "Monat", "Tag",
			"Lat_Sat[°]", "Lon_Sat[°]", "Lat_TP[°]", "Lon_TP[°]",
			"Hoehe_TP[km]", "Erdradius[km]", "Deklinationswinkel[°]",
			"Sonne_Lon[°]", "Zeilendichte[cm^-2]", "Fehler_Zeilendichte[cm^-2]");
	//cerr<<"Matrix schreiben\n";
	for (size_t i = 0; i < lang; i++) {
		//die letzte Zeile der Datei ist leer, da \n in der Vorletzten steht
		fprintf(outfile, "%4i %5i %3i"
				" %1.5E %1.5E %1.5E %1.5E"
				"  %1.5E  %1.5E            %1.5E %1.5E"
				"        %1.5E      %1.5E\n",
				A_Messung_L[i].m_Jahr, A_Messung_L[i].m_Monat, A_Messung_L[i].m_Tag,
				A_Messung_L[i].m_Latitude_Sat, A_Messung_L[i].m_Longitude_Sat,
				A_Messung_L[i].m_Latitude_TP,
				A_Messung_L[i].m_Longitude_TP,
				A_Messung_L[i].m_Hoehe_TP, A_Messung_L[i].m_Erdradius,
				A_Messung_L[i].m_Deklination,
				A_Messung_L[i].m_Sonnen_Longitude,
				A_Messung_L[i].m_Zeilendichte,
				A_Messung_L[i].m_Fehler_Zeilendichten);
	}
	///////////////////////////////////////////////////////////
	//cerr<<"Datei schließen\n";
	// Datei schließen
	fclose(outfile);
}
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//Ende int Ausgabe_Zeilendichten_Limb(string Dateiname);
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////
/* prints the back-inserted columns to a file */
void Ausgabe_Saeulendichten_back(std::string Dateiname,
		std::vector<Ausgewertete_Messung_Limb> &aml_vec, MPL_Matrix &y)
{
	std::vector<Ausgewertete_Messung_Limb> aml_vec_neu;
	std::vector<Ausgewertete_Messung_Limb>::iterator aml_it;

	/* build a new vector with the back-inserted columns
	 * instead of the original ones */
	for (aml_it = aml_vec.begin(); aml_it != aml_vec.end(); ++aml_it) {
		long i = std::distance(aml_vec.begin(), aml_it);
		aml_vec_neu.push_back(*aml_it);
		aml_vec_neu.at(i).m_Zeilendichte = y(i);
	}

	Ausgabe_Saeulendichten(Dateiname, aml_vec_neu);
}
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//START int Ausgabe_Zeilendichten_Nadir
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////
void Ausgabe_Saeulendichten(string Dateiname,
		vector<Ausgewertete_Messung_Nadir> &A_Messung_N)
{
	//Formatierte Ausgabe
	FILE *outfile;
	//Datei öffnen
	outfile = fopen(Dateiname.c_str(), "w");
	//Überschrift
	fprintf(outfile, "%4s %5s %3s "
			"%11s %11s "
			"%11s %11s "
			"%11s %11s %11s"
			"%11s %11s \n",
			"Jahr", "Monat", "Tag",
			"Lat_Sat", "Lon_Sat",
			"Lat_Ground", "Long_Ground",
			"Erdradius", "Deklination[°]", "Sonne_Lon[°]",
			"Säulendichte[cm^2]", "Fehler_Säulendichte[cm^2]");
	size_t lang = A_Messung_N.size();

	for (size_t i = 0; i < lang; i++) {
		//die letzte Zeile der Datei ist leer, da \n in der Vorletzten steht
		fprintf(outfile, "%4i %3i %5i "
				"%1.5E %1.5E "
				"%1.5E %1.5E "
				"%1.5E  %1.5E   %1.5E "
				"    %1.5E       %1.5E\n",
				A_Messung_N[i].m_Jahr, A_Messung_N[i].m_Monat, A_Messung_N[i].m_Tag,
				A_Messung_N[i].m_Latitude_Sat, A_Messung_N[i].m_Longitude_Sat,
				A_Messung_N[i].m_Latitude_Ground,
				A_Messung_N[i].m_Longitude_Ground,
				A_Messung_N[i].m_Erdradius, A_Messung_N[i].m_Deklination,
				A_Messung_N[i].m_Sonnen_Longitude,
				A_Messung_N[i].m_Zeilendichte,
				A_Messung_N[i].m_Fehler_Zeilendichten);
	}
	///////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile);
}
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//Ende int Ausgabe_Zeilendichten_Nadir
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////

MPL_Matrix Read_Atmodatei(string Dateiname)
{
	std::ifstream infile(Dateiname.c_str());
	if (!(infile.is_open())) {
		std::cout << "Datei " << Dateiname << " kann nicht gefunden werden.\n";
		MPL_Matrix dummy;
		return dummy;
	}
	int Zeilenzahl, Spaltenzahl;
	infile >> Zeilenzahl;
	infile >> Spaltenzahl; // in der Datei stehen Spaltenzahl -1 Spalten drin
	Spaltenzahl += 1;
	MPL_Matrix Out(Zeilenzahl, Spaltenzahl);
	for (int i = 0; i < Zeilenzahl; i++) {
		for (int j = 0; j < Spaltenzahl; j++) {
			infile >> Out(i, j);
		}// ende for j
	}//ende for i
	return Out;
}
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//Ende Read_Atmodatei
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Ausgabe_Dichten
////////////////////////////////////////////////////////////////////////////////
void Ausgabe_Dichten(string Dateiname_out, Retrievalgitter &Grid,
		MPL_Matrix &Dichten, MPL_Matrix &Dichten_tot, MPL_Matrix &apriori,
		MPL_Matrix &S_x, MPL_Matrix &S_x_meas, MPL_Matrix &AKM,
		bool save_sx, bool save_akm)
{
	// Die Ausgabe erfolgt in 3 Dateien mit zusätzlichem Namen
	// _Dichten.txt, _Sx.txt und  _AKM.txt
	// Die erste Datei ist die interessanteste davon mit den Dichten.
	// Sie hat folgende Spalten:
	// GP_Nummer Max_H H Min_H Max_Lat Lat Min_Lat Dichte Standardabweichung
	// In den anderen beiden Matrizen stehen jeweils pro Zeile ein Element der
	// Matrizen....  die zuordnung erfolgt aus den Gitterpunktnummern der
	// ersten Datei


	//Formatierte Ausgabe
	FILE *outfile1;
	string Dateiname1, Dateiname2, Dateiname2_meas, Dateiname3;
	Dateiname1 = Dateiname_out + "_Dichten.txt";
	Dateiname2 = Dateiname_out + "_Sx.nc";
	Dateiname2_meas = Dateiname_out + "_Sx_meas.nc";
	Dateiname3 = Dateiname_out + "_AKM.nc";
	int i;
	double stabw = 0;
	//Datei öffnen
	outfile1 = fopen(Dateiname1.c_str(), "w");
	////////////////////////////////////////////////////////////////////////////
	//Überschrift
	fprintf(outfile1, "%5s "
			"%13s %12s %13s "
			"%14s  %12s %14s "
			"  %12s  "
			"%12s %12s %12s %12s %12s %12s\n",
			"GP_ID",
			"Max_Hoehe[km]", "Hoehe[km]", "Min_Hoehe[km]",
			"Max_Breite[°]", "Breite[°]", "Min_Breite[°]",
			"Laenge[°]",
			"Dichte[cm^-3]", "Fehler Mess[cm^-3]",
			"Fehler tot[cm^-3]", "Gesamtdichte[cm^-3]",
			"apriori[cm^-3]", "AKdiag");
	// Alle Zeilen bis auf die letzte
	for (i = 0; i < Grid.m_Anzahl_Punkte; i++) {
		stabw = sqrt(S_x(i, i));
		double stdabw_meas = std::sqrt(S_x_meas(i, i));
		fprintf(outfile1, "%5i  "
				"%+1.5E %+1.5E  %+1.5E "
				" %+1.5E %+1.5E  %+1.5E "
				" %+1.5E  "
				" %+1.5E       %+1.5E      %+1.5E        %+1.5E   %+1.5E   %+1.5E\n",
				i,
				Grid.m_Gitter[i].m_Max_Hoehe, Grid.m_Gitter[i].m_Hoehe, Grid.m_Gitter[i].m_Min_Hoehe,
				Grid.m_Gitter[i].m_Max_Breite, Grid.m_Gitter[i].m_Breite, Grid.m_Gitter[i].m_Min_Breite,
				Grid.m_Gitter[i].longitude,
				Dichten(i), stdabw_meas, stabw, Dichten_tot(i), apriori(i), AKM(i, i));
	}
	////////////////////////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile1);

	// Die anderen beiden kurz und schmerzlos
	//S_x
	// Zeilenweise ausgeben
	if (save_sx) {
		S_x.save_to_netcdf(Dateiname2);
		S_x_meas.save_to_netcdf(Dateiname2_meas);
	}

	//AKM
	// Zeilenweise ausgeben
	if (save_akm)
		AKM.save_to_netcdf(Dateiname3);
}
////////////////////////////////////////////////////////////////////////////////
// Funktionsende Ausgabe_Dichten
////////////////////////////////////////////////////////////////////////////////
