/*
 * Datei_IO.cpp
 *
 *  Created on: 19.04.2010
 *      Author: martin
 */

#include <cmath>
#include "Messung_Limb.h"
#include "Messung_Nadir.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib> //system
#include "Nachricht_Schreiben.h"
#include "MPL_Matrix.h"
#include "LimbNadir_IO.h"
#include "Datei_IO.h"
#include "Retrievalgitter.h"

using namespace std;

extern int Prioritylevel;

// helper function to copy Limb_Datensatz *Limbdaten and
// float *Wellenlaengen into a vector<Messung_Limb>
vector<Messung_Limb> make_messung_limb_vector(string Dateiname,
		Limb_Datensatz *Limbdaten, float *Wellenlaengen,
		int no_of_pix, float orbit_phase, int Datum[6],
		int no_of_heights, int offset, int direction)
{
	// 4. Erstellung des Übergabevektors
	vector<Messung_Limb> Ergebnisvektor;

	for (int i = 0; i < no_of_heights; i++) {
		Messung_Limb ml;
		ml.m_Dateiname_L1C = Dateiname;
		ml.m_Jahr = Datum[0];
		ml.m_Monat = Datum[1];
		ml.m_Tag = Datum[2];
		ml.m_Stunde = Datum[3];
		ml.m_Minute = Datum[4];
		ml.m_orbit_phase = orbit_phase;
		ml.m_Lattidude_Sat = Limbdaten[offset + direction * i].m_Sub_Sat_Lat; //achtung geodätische Koordinaten
		ml.m_Longitude_Sat = Limbdaten[offset + direction * i].m_Sub_Sat_Lon;
		ml.m_Hoehe_Sat = Limbdaten[offset + direction * i].m_Sat_Hoehe;
		ml.m_Lattidude_TP = Limbdaten[offset + direction * i].m_TP_Lat;
		ml.m_Longitude_TP = Limbdaten[offset + direction * i].m_TP_Lon;
		ml.m_Hoehe_TP = Limbdaten[offset + direction * i].m_Tangentenhoehe;
		ml.m_Erdradius = Limbdaten[offset + direction * i].m_Erdradius;
		ml.m_TP_SZA = Limbdaten[offset + direction * i].m_TP_SZA;
		ml.m_Number_of_Wavelength = no_of_pix;

		for (int j = 0; j < no_of_pix; j++) {
			ml.m_Wellenlaengen.push_back(Wellenlaengen[j]);
			ml.m_Intensitaeten.push_back(Limbdaten[offset + direction * i].m_radiance[j]);
			// - Limbdaten[30].m_radiance[j];(nicht gut bei MgI)
			ml.m_Intensitaeten_relativer_Fehler.push_back(Limbdaten[offset + direction * i].m_error[j]);
			// - Limbdaten[30].m_radiance[j];
			ml.m_Sonne.push_back(0.);
			ml.m_Intensitaeten_durch_piF.push_back(0.);
			ml.m_Intensitaeten_durch_piF_Gamma.push_back(0.);
			ml.m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand.push_back(0.);
		}

		// Der Pixel 552 (282,03nm zeigt bei nadir( und nur dort, einen Peak)
		// ....interpretation dead pixel
		ml.m_Intensitaeten.at(536) = 0.5 * (ml.m_Intensitaeten.at(535) + ml.m_Intensitaeten.at(537));
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
		Messung_Limb &Troposphaerische_Saeule, Messung_Limb &mean_10_20)
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
	float *Wellenlaengen;
	Limb_Datensatz *Limbdaten;
	// 2. Laden der Datei
	//cerr<<" 2. Laden der Datei\n";
	Load_Limb_l_mpl_binary(Dateiname,
						   textheader, no_of_alt, no_of_pix, Orbitstate, Datum,
						   Center_Lat_Lon, orbit_phase, Wellenlaengen,
						   Limbdaten);
	// 3. Nachbearbeitung/Ausschlusskriterien
	// Wird jetzt nach dem Laden durchgeführt
	// 4. Erstellung des Übergabevektors
	vector<Messung_Limb> Ergebnisvektor
		= make_messung_limb_vector(Dateiname, Limbdaten, Wellenlaengen,
				no_of_pix, orbit_phase, Datum, 7, 23, 1);

	//Teile von Schritt 4 nochmal für die Troposhärische Säule
	//Eigentlich reichen Intensitäten
	for (int j = 0; j < no_of_pix; j++) {
		Troposphaerische_Saeule.m_Intensitaeten.push_back(Limbdaten[2].m_radiance[j]);
	}
	Troposphaerische_Saeule.m_TP_SZA = Limbdaten[2].m_TP_SZA;
	// und für den Mittelwert, der Höhen 10 bis 20
	for (int j = 0; j < no_of_pix; j++) {
		double mean = 0.;
		for (int k = 10; k < 21; k++) {
			mean += Limbdaten[k].m_error[j];
		}
		mean /= 11.0;
		mean_10_20.m_Intensitaeten.push_back(mean);
	}

	// 5. Speicherfreigabe
	delete[] Wellenlaengen;
	for (int i = 0; i < no_of_alt; i++) {
		if (Limbdaten[i].m_radiance != 0) {
			delete[] Limbdaten[i].m_radiance;
			Limbdaten[i].m_radiance = 0;
		}
		if (Limbdaten[i].m_error != 0) {
			delete[] Limbdaten[i].m_error;
			Limbdaten[i].m_error = 0;
		}
	}
	delete[] Limbdaten;
	// 6. Rückgabe
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
		Messung_Limb &niedrigste_Hoehe)
{
	///////////////////////////////////////////////////////////
	// ähnlich zur üblichen Limbroutine...nur andere TH Reihenfolge und alle
	// Messungen über 70 km werden geladen und die niedrigste bei 53.5 km
	// 0 bis 24   z.B. bei 0 148,7km und bei 25 69.9 km
	///////////////////////////////////////////////////////////

	return ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(Dateiname,
			niedrigste_Hoehe, 25);
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
		Messung_Limb &niedrigste_Hoehe, int Anzahl_Hoehen)
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
	float *Wellenlaengen;
	Limb_Datensatz *Limbdaten;
	// 2. Laden der Datei
	//cerr<<" 2. Laden der Datei\n";
	Load_Limb_l_mpl_binary(Dateiname,
						   textheader, no_of_alt, no_of_pix, Orbitstate, Datum,
						   Center_Lat_Lon, orbit_phase, Wellenlaengen,
						   Limbdaten);
	// 3. Nachbearbeitung/Ausschlusskriterien
	// Wird jetzt nach dem Laden durchgeführt
	// 4. Erstellung des Übergabevektors
	vector<Messung_Limb> Ergebnisvektor
		= make_messung_limb_vector(Dateiname, Limbdaten, Wellenlaengen,
				no_of_pix, orbit_phase, Datum, Anzahl_Hoehen, 24, -1);

	//Teile von Schritt 4 nochmal für die niedrigste Höhe
	//Eigentlich reichen Intensitäten
	for (int j = 0; j < no_of_pix; j++) {
		niedrigste_Hoehe.m_Intensitaeten.push_back(Limbdaten[29].m_radiance[j]);
	}
	niedrigste_Hoehe.m_TP_SZA = Limbdaten[29].m_TP_SZA;
	// 5. Speicherfreigabe
	delete[] Wellenlaengen;
	for (int i = 0; i < no_of_alt; i++) {
		if (Limbdaten[i].m_radiance != 0) {
			delete[] Limbdaten[i].m_radiance;
			Limbdaten[i].m_radiance = 0;
		}
		if (Limbdaten[i].m_error != 0) {
			delete[] Limbdaten[i].m_error;
			Limbdaten[i].m_error = 0;
		}
	}
	delete[] Limbdaten;
	// 6. Rückgabe
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
		Nadir_Datensatz *Nadirdaten, float *Wellenlaenge,
		int No_of_Messungen, int No_of_Pix, int *Kanal_Nr)
{
	// 4. Erstellung des Übergabevektors
	vector<Messung_Nadir> mn_vec;

	for (int i = 0; i < No_of_Messungen; i++) {
		Messung_Nadir mn;
		//Herkunftsmerkmale
		mn.m_Dateiname_L1C = Dateiname;
		mn.m_Messung_ID = Nadirdaten[i].m_Messung_ID;
		//Datum
		mn.m_Jahr = Nadirdaten[i].m_Jahr;
		mn.m_Monat = Nadirdaten[i].m_Monat;
		mn.m_Tag = Nadirdaten[i].m_Tag;
		mn.m_Stunde = Nadirdaten[i].m_Stunde;
		mn.m_Minute = Nadirdaten[i].m_Minute;
		//Geolokationen
		mn.m_Lattitude_Sat = Nadirdaten[i].m_Sat_Lat;
		mn.m_Longitude_Sat = Nadirdaten[i].m_Sat_Lon;
		mn.m_Hoehe_Sat = Nadirdaten[i].m_Hoehe;
		mn.m_Lattitude_Ground = Nadirdaten[i].m_geo_nadir_center_lat;
		mn.m_Longitude_Ground = Nadirdaten[i].m_geo_nadir_center_lon;
		mn.m_Erdradius = Nadirdaten[i].m_Sat_Erdradius;
		mn.m_orbit_phase = Nadirdaten[i].m_orbit_phase;

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
		}
		// Der Pixel 552 (282,03nm zeigt bei nadir( und nur dort, einen Peak)
		// Das ist die Kanalgrenze zwischen Kanal1a und Kanal1b
		// Interpolieren zwischen 282 und 283 nm
		//mn.m_Intensitaeten[536]=(mn.m_Intensitaeten[535]+mn.m_Intensitaeten[537])/2;
		//mn.m_Intensitaeten_relativer_Fehler[536]=(mn.m_Intensitaeten_relativer_Fehler[535]+mn.m_Intensitaeten_relativer_Fehler[537])/2;
		double f1 = 0.9, f2 = 0.1;
		for (int i = 536; i < 545; i++) {
			mn.m_Intensitaeten[i]
				= f1 * mn.m_Intensitaeten[535] + f2 * mn.m_Intensitaeten[545];
			mn.m_Intensitaeten_relativer_Fehler[536]
				= f1 * mn.m_Intensitaeten_relativer_Fehler[535]
				+ f2 * mn.m_Intensitaeten_relativer_Fehler[545];
			f1 -= 0.1;
			f2 += 0.1;
		}

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
	int *Kanal_Nr = 0;
	float *Wellenlaenge = 0;
	Nadir_Datensatz *Nadirdaten = 0;
	// 2. Laden der Datei
	Load_Nadir_n_mpl_binary(Dateiname,
							textheader, No_of_Messungen, No_of_Pix,
							Kanal_Nr, Wellenlaenge, Nadirdaten);
	// 3. Nachbearbeitung/Ausschlusskriterien
	// Dieser Schritt wird erst nach dem Laden aufgerufen
	// 4. Erstellung des Übergabefelds
	Anzahl_Messungen = No_of_Messungen;
	vector<Messung_Nadir> aus
		= make_messung_nadir_vector(Dateiname, Nadirdaten, Wellenlaenge,
				No_of_Messungen, No_of_Pix, Kanal_Nr);

	// 5. Speicherfreigabe
	delete[] Kanal_Nr;
	delete[] Wellenlaenge;
	for (int i = 0; i < No_of_Messungen; i++) {
		if (Nadirdaten[i].m_radiance != 0) {
			delete[] Nadirdaten[i].m_radiance;
			Nadirdaten[i].m_radiance = 0;
		}
		if (Nadirdaten[i].m_error != 0) {
			delete[] Nadirdaten[i].m_error;
			Nadirdaten[i].m_error = 0;
		}
	}
	delete[] Nadirdaten;
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
int Ausgabe_Saeulendichten(string Dateiname,
		vector<Ausgewertete_Messung_Limb> A_Messung_L)
{
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
	int lang = A_Messung_L.size();
	//cerr<<"Matrix schreiben\n";
	for (int i = 0; i < lang; i++) {
		//die letzte Zeile der Datei ist leer, da \n in der Vorletzten steht
		fprintf(outfile, "%4i %5i %3i"
				" %1.5E %1.5E %1.5E %1.5E"
				"  %1.5E  %1.5E            %1.5E %1.5E"
				"        %1.5E      %1.5E\n",
				A_Messung_L[i].m_Jahr, A_Messung_L[i].m_Monat, A_Messung_L[i].m_Tag,
				A_Messung_L[i].m_Lattidude_Sat, A_Messung_L[i].m_Longitude_Sat,
				A_Messung_L[i].m_Lattidude_TP,
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
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//Ende int Ausgabe_Zeilendichten_Limb(string Dateiname);
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//START int Ausgabe_Zeilendichten_Nadir
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////
int Ausgabe_Saeulendichten(string Dateiname,
		vector<Ausgewertete_Messung_Nadir> A_Messung_N)
{
	//Formatierte Ausgabe
	FILE *outfile;
	//Datei öffnen
	outfile = fopen(Dateiname.c_str(), "w");
	//Überschrift
	fprintf(outfile, "%4s %3s %5s "
			"%11s %11s "
			"%11s %11s "
			"%11s %11s %11s"
			"%11s %11s \n",
			"Jahr", "Tag", "Monat",
			"Lat_Sat", "Lon_Sat",
			"Lat_Ground", "Long_Ground",
			"Erdradius", "Deklination[°]", "Sonne_Lon[°]",
			"Säulendichte[cm^2]", "Fehler_Säulendichte[cm^2]");
	int lang = A_Messung_N.size();

	for (int i = 0; i < lang; i++) {
		//die letzte Zeile der Datei ist leer, da \n in der Vorletzten steht
		fprintf(outfile, "%4i %3i %5i "
				"%1.5E %1.5E "
				"%1.5E %1.5E "
				"%1.5E  %1.5E   %1.5E "
				"    %1.5E       %1.5E\n",
				A_Messung_N[i].m_Jahr, A_Messung_N[i].m_Monat, A_Messung_N[i].m_Tag,
				A_Messung_N[i].m_Lattitude_Sat, A_Messung_N[i].m_Longitude_Sat,
				A_Messung_N[i].m_Lattitude_Ground,
				A_Messung_N[i].m_Longitude_Ground,
				A_Messung_N[i].m_Erdradius, A_Messung_N[i].m_Deklination,
				A_Messung_N[i].m_Sonnen_Longitude,
				A_Messung_N[i].m_Zeilendichte,
				A_Messung_N[i].m_Fehler_Zeilendichten);
	}
	///////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile);
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
//********************************************************************//
//Ende int Ausgabe_Zeilendichten_Nadir
//********************************************************************//
////////////////////////////////////////////////////////////////////////////////

MPL_Matrix Read_Atmodatei(string Dateiname)
{
	ifstream infile;
	infile.open(Dateiname.c_str());
	if (!(infile.is_open())) {
		cout << "Datei " << Dateiname << " kann nicht gefunden werden.\n";
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
int Ausgabe_Dichten(string Dateiname_out, Retrievalgitter Grid,
		MPL_Matrix Dichten, MPL_Matrix S_x, MPL_Matrix AKM)
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
	ofstream outfile2, outfile3;
	string Dateiname1, Dateiname2, Dateiname3;
	Dateiname1 = Dateiname_out + "_Dichten.txt";
	Dateiname2 = Dateiname_out + "_Sx.txt";
	Dateiname3 = Dateiname_out + "_AKM.txt";
	int Anzahl = Grid.m_Anzahl_Breiten * Grid.m_Anzahl_Hoehen;
	int i, j;
	double stabw = 0;
	//Datei öffnen
	outfile1 = fopen(Dateiname1.c_str(), "w");
	////////////////////////////////////////////////////////////////////////////
	//Überschrift
	fprintf(outfile1, "%5s "
			"%13s %12s %13s "
			"%14s  %12s %14s "
			"%12s %12s\n",
			"GP_ID",
			"Max_Hoehe[km]", "Hoehe[km]", "Min_Hoehe[km]",
			"Max_Breite[°]", "Breite[°]", "Min_Breite[°]",
			"Dichte[cm^-3]", " Standardabweichung[cm^-3]");
	// Alle Zeilen bis auf die letzte
	for (i = 0; i < Grid.m_Anzahl_Punkte - 1; i++) {
		stabw = sqrt(S_x(i, i));
		fprintf(outfile1, "%5i  "
				"%+1.5E %+1.5E  %+1.5E "
				" %+1.5E %+1.5E  %+1.5E "
				" %+1.5E               %+1.5E\n",
				i,
				Grid.m_Gitter[i].m_Max_Hoehe, Grid.m_Gitter[i].m_Hoehe, Grid.m_Gitter[i].m_Min_Hoehe,
				Grid.m_Gitter[i].m_Max_Breite, Grid.m_Gitter[i].m_Breite, Grid.m_Gitter[i].m_Min_Breite,
				Dichten(i), stabw);


	}
	//letzte Zeile (ohne \n am Ende)
	i = Anzahl - 1;
	stabw = sqrt(S_x(i, i));
	fprintf(outfile1, "%5i  "
			"%+1.5E %+1.5E  %+1.5E "
			" %+1.5E %+1.5E  %+1.5E "
			" %+1.5E               %+1.5E\n",
			i,
			Grid.m_Gitter[i].m_Max_Hoehe, Grid.m_Gitter[i].m_Hoehe, Grid.m_Gitter[i].m_Min_Hoehe,
			Grid.m_Gitter[i].m_Max_Breite, Grid.m_Gitter[i].m_Breite, Grid.m_Gitter[i].m_Min_Breite,
			Dichten(i), stabw);
	////////////////////////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile1);
	// Die anderen beiden kurz und schmerzlos
	//S_x
	// Zeilenweise ausgeben
	outfile2.open(Dateiname2.c_str());
	for (i = 0; i < Anzahl; i++) {
		for (j = 0; j < Anzahl; j++) {
			outfile2 << S_x(i, j);
			if (j != Anzahl - 1)
				outfile2 << "\t";
		}// for j
		if (i != Anzahl - 1)
			outfile2 << "\n";
	}//for i
	outfile2.close();
	outfile2.clear();
	//AKM
	// Zeilenweise ausgeben
	outfile3.open(Dateiname3.c_str());
	for (i = 0; i < Anzahl; i++) {
		for (j = 0; j < Anzahl; j++) {
			outfile3 << AKM(i, j);
			if (j != Anzahl - 1)
				outfile3 << "\t";
		}// for j
		if (i != Anzahl - 1)
			outfile3 << "\n";
	}//for i
	outfile3.close();
	outfile3.clear();
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// Funktionsende Ausgabe_Dichten
////////////////////////////////////////////////////////////////////////////////
