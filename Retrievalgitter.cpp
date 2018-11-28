/*
 * Retreivalgitter.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 27.05.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include "Retrievalgitter.h"

#include "Ausgewertete_Messung_Limb.h"
#include "Konfiguration.h"
#include "Glaetten.h"

#include <vector>
#include <cmath>
#include "MPL_Vektor.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>


using std::vector;

//Methoden
//////////////////////////////////////////////////////////
//
// Funktionsstart Retrievalgitter_erzeugen
//
/////////////////////////////////////////////////////////
void Retrievalgitter::Retrievalgitter_erzeugen(
		vector<Ausgewertete_Messung_Limb> &AM_Limb, double Epsilon,
		Konfiguration &Konf)
{
	// deg -> rad conversion factor
	const double deg = M_PI / 180.;
	// emprical longitude factor
	const double phi_fac = 11.91;
	// Das Epsilon gibt den Mindestabstand zweier Gitterpunkte in Breitengrad an

	// ACHTUNG falls die Lats beim auf und absteigen des Satelliten zufällig
	// gleich sind, gibt es Fehler....  das ist aber extrem
	// unwahrscheinlich...davon sollte man sich aber auch nochmal überzeugen

	// Zur Zeitaufwandsabschätzung: Das Gitter wird am Ende aus Rund 20 Lats
	// bestehen, und für jede Nord_Sued Bestimmung werden auch weniger als 1000
	// Rechenschritte gemacht
	//Zunächst Latitudes in einen Vektor schreiben, aber nur die, welche noch
	//nicht vorgekommen sind Es wird davon ausgegangen, dass der Datensatz nach
	//Zeit sortiert ist
	vector<double> Lats_Messung;
	vector<double> CLons0, SLons0;
	vector<double>::iterator lat_it;
	vector<Ausgewertete_Messung_Limb>::iterator aml_it;

	std::cout << "# measurement geolocations begin" << std::endl;
	for (aml_it = AM_Limb.begin(); aml_it != AM_Limb.end(); ++aml_it) {
		bool doppelt = false;

		for (lat_it = Lats_Messung.begin(); lat_it != Lats_Messung.end(); ++lat_it) {
			if (abs(aml_it->center_lat - (*lat_it)) < Epsilon) {
				doppelt = true;
				break;
			}
		}

		/* For each "unique" point we save the latitude and the calculated
		 * equatorial longitude lon0 according to: lon = C * tg(lat) + lon0,
		 * which translates to: lon0 = lon - C * tg(lat).
		 * (C is our experimental factor 11.91.)
		 * Simply interpolating to latitude = 0 sometimes does not work because
		 * of removed scans, e.g. from the SAA, or because of lack of data.
		 * To reconstruct the equatorial longitude of the orbit, we then calculate
		 * the mean of all these calculated values. To avoid periodicity
		 * problems around 0/360 degrees, we actually save the sine and cosine
		 * values and use atan2() later. */
		if (!(doppelt)) {
			double phi0 = aml_it->center_lon -
				phi_fac * std::tan(aml_it->center_lat * deg);
			Lats_Messung.push_back(aml_it->center_lat);
			CLons0.push_back(std::cos(phi0 * deg));
			SLons0.push_back(std::sin(phi0 * deg));
			std::cout << aml_it->center_lon << "\t"
				<< Lats_Messung.back() << std::endl;
		}
	}
	std::cout << "# measurement geolocations end" << std::endl;
	// Lats_Messung enthält nun alle Lats nur einmal und in Zeitgeordneter
	// Reihenfolge Der Satellit auf dem Stück zwischen Maximaler und Minimaler
	// Höhe in Nord nach Süd Richtung
	// Durch die Ekliptik der Erde und der zugehörigen Deklination der
	// Sonnenzenitbreite, gibt es im Sommer am Anfang der Messreihe Messpunkte
	// wo der Satellit von Süd nach Nord fliegt, und schon messen kann
	// (Polartag am Nordpol).
	// Im Winter gibt es dementsprechend zusätzliche Süd-Nord Messungen am Ende
	// der Messreihe Wir suchen nun zunächst den index der Maximalen und der
	// Minimalen Lat um das Gitter zu bauen
	int Max_Index = Get_Index_of_Maximum(Lats_Messung);
	int Min_Index = Get_Index_of_Minimum(Lats_Messung, Max_Index);
	// Wir bauen ein äquidistantes Gitter zwischen diesen beiden Punkten auf,
	// dafür bestimmen wir zunächst die Gitterkonstante
	double MaxLat = Lats_Messung[Max_Index];
	double MinLat = Lats_Messung[Min_Index];
	int Breitenzahl = Min_Index - Max_Index + 1;
	/* Reconstruct the equatorial longitude from the mean of the individual
	 * scan values (sine and cosine) saved above.
	 * We can just sum here because we divide the two values, i.e. any
	 * constant factor (e.g. CLons0.size() for the real mean) drops out.
	 * We remove the endpoints to avoid that the large outliers at the edge of
	 * the orbits skew the mean too much. */
	double CLon0s = std::accumulate(CLons0.begin() + 1, CLons0.end() - 1, 0.);
	double SLon0s = std::accumulate(SLons0.begin() + 1, SLons0.end() - 1, 0.);
	double lon0 = std::fmod(
			std::atan2(SLon0s, CLon0s) + 2. * M_PI, 2. * M_PI) / deg;

	/////////////////////////////
	// selbst gesetzt.....
	MaxLat = Konf.m_MaxLat;
	MinLat = Konf.m_MinLat;
	Breitenzahl = Konf.m_NLat; //20
	if (Breitenzahl == 0) Breitenzahl = 1;
	double Gitterkonstante = (MaxLat - MinLat) / (double)(Breitenzahl);
	// adjust grid box centre latitudes
	MinLat += 0.5 * Gitterkonstante;
	MaxLat -= 0.5 * Gitterkonstante;
	std::cout << "# grid: maxlat = " << MaxLat << ", minlat = " << MinLat << std::endl;
	std::cout << "# lon0 = " << lon0 << std::endl;
	// Hoeheneinteilung
	// Bisherige Werte für: (mehr Höhen, mehr Rechenzeit,
	// mehr unterbestimmtheit der Gelichungssysteme
	//  MgI 82 Gitterpunkte bis 150 km in 1km Schritten
	//  MgII 132 Gitterpunkte bis 200 km in 1km Schritten
	//  unbekannte Spezies bei niedrigen Hoehen......teste bis 110
	double h_start = Konf.m_MinAlt;
	double h_end = Konf.m_MaxAlt;
	double dh = Konf.m_dAlt, dhh = 0.5 * dh;
	int Anzahl_Hoehen = int((h_end - h_start) / dh) + 1;   // 82  //42  //132
	int add_alts = Konf.m_Anzahl_zusaetzliche_Hoehengitterpunkte;
	/* Check if the extensions really start above the top grid altitude.
	 * If not, simply skip the extensions entirely. This ensures backward
	 * compatibility with old code versions and config files where the number
	 * of grid extensions is non-zero but starting at the ground. Since we did
	 * not extend the grid at all before, we use this as an indicator to not
	 * extend the grid. */
	if (add_alts > 0 && Konf.m_Grid_ext_low.at(0) <= h_end)
		add_alts = 0;
	vector<double> untere_Hoehe(Anzahl_Hoehen + add_alts);
	vector<double> mittlere_Hoehe(Anzahl_Hoehen + add_alts);
	vector<double> obere_Hoehe(Anzahl_Hoehen + add_alts);
	for (int i = 0; i < Anzahl_Hoehen; i++) {
		// dh km Schritte
		untere_Hoehe[i]   = h_start - dhh + i * dh;
		mittlere_Hoehe[i] = h_start + i * dh;
		obere_Hoehe[i]    = h_start + dhh + i * dh;
	}
	/////////////////////////////
	for (int i = 0; i < add_alts; ++i) {
		untere_Hoehe[Anzahl_Hoehen + i] = Konf.m_Grid_ext_low.at(i);
		obere_Hoehe[Anzahl_Hoehen + i] = Konf.m_Grid_ext_high.at(i);
		mittlere_Hoehe[Anzahl_Hoehen + i]
			= 0.5 * (obere_Hoehe[Anzahl_Hoehen + i]
					+ untere_Hoehe[Anzahl_Hoehen + i]);
	}
	Anzahl_Hoehen += add_alts;

	this->m_Anzahl_Hoehen = Anzahl_Hoehen;
	this->m_Anzahl_Breiten = Breitenzahl;
	std::cout << "m_Anzahl_Breiten: " << m_Anzahl_Breiten << "\n";
	std::cout << "m_Anzahl_Hoehen: " << m_Anzahl_Hoehen << "\n";
	this->m_Anzahl_Punkte = m_Anzahl_Breiten * m_Anzahl_Hoehen;
	// Es muss noch Speicherplatz für das Gitter reserviert werden ////////
	m_Gitter.resize(m_Anzahl_Punkte);
	Gitterpunkt GP;
	//Gitterpunkte des Gitters erzeugen sortiert nach Höhen und Breiten
	//cout<<"Anfang for i\n";
	for (int i = 0; i < Breitenzahl; i++) { //Breiten von Nord nach Süd
		//cout<<"Anfang for j\n";
		for (int j = 0; j < Anzahl_Hoehen; j++) { //Höhen von unten nach oben
			//Index des Punktes
			// so sind alle Höhen in einer Zeile und Breiten in einer Spalte
			GP.m_eigener_Index = i + j * Breitenzahl;
			// z.b. erste Zeile j==0; d.h. Höhe konstant und i wandert
			// (die schleife ist etwas unglücklich,
			//  weil sie Spaltenweise auffüllt, nicht
			//  zeileweise, wie sonst in c++ üblich)
			// Index der Nachbarpunkte...falls Nachbar nicht existiert..-> -1
			if ((i - 1) >= 0) {
				GP.m_Index_Nord_Nachbar = (i - 1) + j * Breitenzahl;   //Nord
			} else {
				GP.m_Index_Nord_Nachbar = -1;
			}
			if ((i + 1) < Breitenzahl) {
				GP.m_Index_Sued_Nachbar = (i + 1) + j * Breitenzahl;   //Sued
			} else {
				GP.m_Index_Sued_Nachbar = -1;
			}
			if ((j - 1) >= 0) {
				GP.m_Index_unterer_Nachbar = i + (j - 1) * Breitenzahl; //unten
			} else {
				GP.m_Index_unterer_Nachbar = -1;
			}
			if ((j + 1) < Anzahl_Hoehen) {
				GP.m_Index_oberer_Nachbar = i + (j + 1) * Breitenzahl;   //oben
			} else {
				GP.m_Index_oberer_Nachbar = -1;
			}
			// und die Diagonalen Punkte
			if (((i - 1) >= 0) && ((j + 1) < Anzahl_Hoehen)) {
				GP.m_Index_oberer_Nord_Nachbar
					= (i - 1) + (j + 1) * Breitenzahl;   //Nord oben
			} else {
				GP.m_Index_oberer_Nord_Nachbar = -1;
			}
			if (((i - 1) >= 0) && ((j - 1) >= 0)) {
				GP.m_Index_unterer_Nord_Nachbar
					= (i - 1) + (j - 1) * Breitenzahl;   //Nord unten
			} else {
				GP.m_Index_unterer_Nord_Nachbar = -1;
			}
			if (((i + 1) < Breitenzahl) && ((j + 1) < Anzahl_Hoehen)) {
				GP.m_Index_oberer_Sued_Nachbar
					= (i + 1) + (j + 1) * Breitenzahl;   //Süd oben
			} else {
				GP.m_Index_oberer_Sued_Nachbar = -1;
			}
			if (((i + 1) < Breitenzahl) && ((j - 1) >= 0)) {
				GP.m_Index_unterer_Sued_Nachbar
					= (i + 1) + (j - 1) * Breitenzahl;   //Süd unten
			} else {
				GP.m_Index_unterer_Sued_Nachbar = -1;
			}
			//  3 Hoehen
			GP.m_Hoehe = mittlere_Hoehe[j];
			GP.m_Max_Hoehe = obere_Hoehe[j];
			GP.m_Min_Hoehe = untere_Hoehe[j];
			//  3 Breiten
			GP.m_Breite = MaxLat - i * Gitterkonstante;
			GP.m_Max_Breite = GP.m_Breite + 0.5 * Gitterkonstante;
			GP.m_Min_Breite = GP.m_Breite - 0.5 * Gitterkonstante;
			// calculate the grid point longitude
			GP.longitude = std::tan(GP.m_Breite * deg) * phi_fac + lon0;
			//SZA initialisieren
			// Die beiden gibts nicht mehr
			// TODO entgültig löschen
			//GP.m_SZA=0;
			//GP.m_Streuwinkel=0;
			// Am Ende Gitterpunkt übergeben
			this->m_Gitter[GP.m_eigener_Index] = GP;
		}// for j ende
	}// for i ende
	//cout<<"Ende Retrievalgitter erzeugen\n";
}
//////////////////////////////////////////////////////////
// ENDE Retrievalgitter_erzeugen
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// Funktionsstart Alle_Durchstosspunkte_Null_setzen
//////////////////////////////////////////////////////////
void Retrievalgitter::Alle_Durchstosspunkte_Null_setzen()
{
	std::vector<Gitterpunkt>::iterator gp;
	for (gp = m_Gitter.begin(); gp != m_Gitter.end(); ++gp) {
		gp->m_vorderer_Durchstosspunkt.Null_Initialisierung();
		gp->m_hinterer_Durchstosspunkt.Null_Initialisierung();
	}
}
//////////////////////////////////////////////////////////
// ENDE Alle_Durchstosspunkte_Null_setzen
//////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////
// Funktionsstart In_Datei_Ausgeben
/////////////////////////////////////////////////////////
void Retrievalgitter::In_Datei_Ausgeben(std::string Dateiname)
{
	//Gebe alle Gitterpunkte aus
	int i;
	int Anzahl = this->m_Anzahl_Punkte;
	FILE *outfile;
	//Datei öffnen
	outfile = fopen(Dateiname.c_str(), "w");
	//Überschrift
	fprintf(outfile, "%11s " \
			"%11s %11s %11s "\
			"%11s %11s "\
			"%11s %11s %11s "\
			"%11s %11s %11s "\
			"%11s %11s %11s "\
			"%11s %11s %11s "\
			"%11s %11s %11s\n", \
			"GP_Index", \
			"Nachbar_oben_Nord", "Nachbar_oben", "Nachbar_oben_Süd", \
			"Nachbar_Nord", "Nachbar_Süd", \
			"Nachbar_unten_Nord", "Nachbar_unten", "Nachbar_unten_Süd", \
			"Max_Hoehe[km]", "Hoehe[km]", "Min_Hoehe[km]", \
			"Max_Breite[°]", "Breite[°]", "Min_Breite[°]", \
			"Durchstoß1_x", "Durchstoß1_y", "Durchstoß1_z", \
			"Durchstoß2_x", "Durchstoß2_y", "Durchstoß2_z");
	// Alle Zeilen bis auf die letzte
	for (i = 0; i < Anzahl - 1; i++) {
		fprintf(outfile, "%4i "
				"%4i %4i %4i "
				"%4i %4i "
				"%4i %4i %4i "
				"%1.5E %1.5E %1.5E "
				"%1.5E %1.5E %1.5E "
				"%1.5E %1.5E %1.5E "
				"%1.5E %1.5E %1.5E\n",
				m_Gitter[i].m_eigener_Index,
				m_Gitter[i].m_Index_oberer_Nord_Nachbar,
				m_Gitter[i].m_Index_oberer_Nachbar,
				m_Gitter[i].m_Index_oberer_Sued_Nachbar,
				m_Gitter[i].m_Index_Nord_Nachbar,
				m_Gitter[i].m_Index_Sued_Nachbar,
				m_Gitter[i].m_Index_unterer_Nord_Nachbar,
				m_Gitter[i].m_Index_unterer_Nachbar,
				m_Gitter[i].m_Index_unterer_Sued_Nachbar,
				m_Gitter[i].m_Max_Hoehe, m_Gitter[i].m_Hoehe,
				m_Gitter[i].m_Min_Hoehe,
				m_Gitter[i].m_Max_Breite, m_Gitter[i].m_Breite,
				m_Gitter[i].m_Min_Breite,
				m_Gitter[i].m_vorderer_Durchstosspunkt(0),
				m_Gitter[i].m_vorderer_Durchstosspunkt(1),
				m_Gitter[i].m_vorderer_Durchstosspunkt(2),
				m_Gitter[i].m_hinterer_Durchstosspunkt(0),
				m_Gitter[i].m_hinterer_Durchstosspunkt(1),
				m_Gitter[i].m_hinterer_Durchstosspunkt(2));
	}
	//letzte Zeile (ohne \n am Ende)
	i = Anzahl - 1;
	fprintf(outfile, "%4i "
			"%4i %4i %4i "
			"%4i %4i "
			"%4i %4i %4i "
			"%1.5E %1.5E %1.5E "
			"%1.5E %1.5E %1.5E "
			"%1.5E %1.5E %1.5E "
			"%1.5E %1.5E %1.5E\n",
			m_Gitter[i].m_eigener_Index,
			m_Gitter[i].m_Index_oberer_Nord_Nachbar,
			m_Gitter[i].m_Index_oberer_Nachbar,
			m_Gitter[i].m_Index_oberer_Sued_Nachbar,
			m_Gitter[i].m_Index_Nord_Nachbar, m_Gitter[i].m_Index_Sued_Nachbar,
			m_Gitter[i].m_Index_unterer_Nord_Nachbar,
			m_Gitter[i].m_Index_unterer_Nachbar,
			m_Gitter[i].m_Index_unterer_Sued_Nachbar,
			m_Gitter[i].m_Max_Hoehe, m_Gitter[i].m_Hoehe,
			m_Gitter[i].m_Min_Hoehe,
			m_Gitter[i].m_Max_Breite, m_Gitter[i].m_Breite,
			m_Gitter[i].m_Min_Breite,
			m_Gitter[i].m_vorderer_Durchstosspunkt(0),
			m_Gitter[i].m_vorderer_Durchstosspunkt(1),
			m_Gitter[i].m_vorderer_Durchstosspunkt(2),
			m_Gitter[i].m_hinterer_Durchstosspunkt(0),
			m_Gitter[i].m_hinterer_Durchstosspunkt(1),
			m_Gitter[i].m_hinterer_Durchstosspunkt(2));
	// Datei schließen
	fclose(outfile);
}
//////////////////////////////////////////////////////////
// ENDE In_Datei_Ausgeben
/////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// Funktionsstart SZA_und_Streuwinkel_NULL_setzen
//////////////////////////////////////////////////////////
/*void Retrievalgitter::SZA_und_Streuwinkel_NULL_setzen()
{
// Bestimmung der Winkel während des Raytracings alle x Schritte sinnvoller

    for(int i=0; i<this->m_Anzahl_Punkte;i++)
    {
        this->m_Gitter[i].m_SZA=0;
        this->m_Gitter[i].m_Streuwinkel=0;
    }
}*/
//////////////////////////////////////////////////////////
// ENDE SZA_und_Streuwinkel_NULL_setzen
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// Funktionsstart SZA_und_Streuwinkel_bestimmen
//////////////////////////////////////////////////////////
/*void Retrievalgitter::SZA_und_Streuwinkel_bestimmen(double Sat_Lat,
 *   double Sat_Lon,double Sat_Hoehe,
    //Richtpunkt (Tangentenpunkt oder Grundpunkt)
     double RP_Lat,double RP_Lon, double RP_Hoehe,
     double Erdradius, double Deklination, int Stunde, int Minute)
{
    // Das birgt gerade an den Umkehrpunkten große Ungenauigkeiten,
    // den Längengrad der Box zu bestimmen
    // Bestimmung der Winkel während des Raytracings alle x Schritte sinnvoller
    const double pi=3.1415926535897932;
    //Sonnentageszeitwinkel bestimmen(Längengrad der Sonne)
    double Sonnentageszeitwinkel; //in Grad
    double d_Stunden;
    d_Stunden=Stunde+(double) Minute/60.0;
        Sonnentageszeitwinkel= 180-360*(d_Stunden/24.0);
    double Sin_Lat_OP;         double Cos_Lat_OP;            // OP wie Ortspunkt
    double Sin_Lon_OP;        double Cos_Lon_OP;
    double Sin_Lat_Sat;        double Cos_Lat_Sat;              //
    double Sin_Lon_Sat;       double Cos_Lon_Sat;
    double Sin_Lat_RP;        double Cos_Lat_RP;            // RP wie Richtpunkt
    double Sin_Lon_RP;       double Cos_Lon_RP;
    double Sin_Lat_Sonne;    double Cos_Lat_Sonne;
    double Sin_Lon_Sonne;   double Cos_Lon_Sonne;
    // Hier weiter
    Sin_Lat_Sat =sin(*pi/180.0);    Cos_Lat_Sat=cos(*pi/180.0);
    Sin_Lon_Sat =sin(*pi/180.0);    Cos_Lon_Sat=cos(*pi/180.0);


    Sin_Lat_Sonne =sin(Deklination*pi/180.0);
    Cos_Lat_Sonne=cos(Deklination*pi/180.0);
    Sin_Lon_Sonne =sin(Sonnentageszeitwinkel*pi/180.0);
    Cos_Lon_Sonne=cos(Sonnentageszeitwinkel*pi/180.0);
    //Ortsvektor des Satelliten,
    //des Richtpunkt und Einheitsvektor der auf Verbindungsachse
    MPL_Vektor Sat, RP, e_Sat_RP;

    //Schleife über alle Gitterpunkte
    for(int i=0; i< this->m_Anzahl_Punkte;i++)
    {
        // erstmal Längengrad des Gitterpunkts bestimmen

        // Sin_Lat_OP      =sin(this->m_Gitter[i].*pi/180.0);
        // Cos_Lat_OP     =cos(*pi/180.0);
        // Sin_Lon_OP     =sin(*pi/180.0);
        // Cos_Lon_OP    =cos(*pi/180.0);
        //SZA bestimmen
        //Das Skalarprodukt des Ortsvektors mit dem Sonnenort
        //ergibt gerade den Kosinus des gesuchten Winkels

        //Streuwinkel bestimmen
    }

}*/
//////////////////////////////////////////////////////////
// ENDE SZA_und_Streuwinkel_bestimmen
//////////////////////////////////////////////////////////


/*  Veraltet...kann nach Test gelöscht werden
 * void Retrievalgitter::Retrievalgitter_erzeugen_und_Messrichtung_herausfinden(
 *   vector<Ausgewertete_Messung_Limb>& AM_Limb)
{
    // ACHTUNG falls die Lats beim auf und absteigen des Satelliten zufällig
    // gleich sind, gibt es Fehler....  das ist aber extrem
    // unwahrscheinlich...davon sollte man sich aber auch nochmal überzeugen

    // Zur Zeitaufwandsabschätzung: Das Gitter wird am Ende aus Rund 20 Lats
    // bestehen, und für jede Nord_Sued Bestimmung werden auch weniger als 1000
    // Rechenschritte gemacht

    //Zunächst Latitudes in einen Vektor schreiben, aber nur die, welche noch
    //nicht vorgekommen sind Es wird davon ausgegangen, dass der Datensatz nach
    //Zeit sortiert ist
    vector<double> Lats_Messung;
    for(unsigned int i=0;i<AM_Limb.size();i++)
    {
        bool doppelt=false;
        for(unsigned int j=0;j<Lats_Messung.size();j++)
        {
            if(AM_Limb[i].m_Lat_TP==Lats_Messung[j])
            {
                doppelt=true;
                break;
            }
        }
        if(!(doppelt))
        {
            Lats_Messung.push_back(AM_Limb[i].m_Lat_TP);
        }
    }// ende for i
    // Lats_Messung enthält nun alle Lats nur einmal und in Zeitgeordneter
    // Reihenfolge Der Satellit auf dem Stück zwischen Maximaler und Minimaler
    // Höhe in Nord nach Süd Richtung
    // Durch die Ekliptik der Erde und der zugehörigen Deklination der
    // Sonnenzenitbreite, gibt es im Sommer am Anfang der Messreihe Messpunkte
    // wo der Satellit von Süd nach Nord fliegt, und schon messen kann
    // (Polartag am Nordpol).
    // Im Winter gibt es dementsprechend zusätzliche Süd-Nord Messungen am Ende
    // der Messreihe Wir suchen nun zunächst den index der Maximalen und der
    // Minimalen Lat um das Gitter zu bauen
    int Max_Index =Get_Index_of_Maximum(Lats_Messung);
    int Min_Index  =Get_Index_of_Minimum(Lats_Messung);
    // Wir bauen nun ein äquidistantes Gitter zwischen diesen beiden Punkten
    // auf, dafür bestimmen wir zunächst die Gitterkonstante
    double MaxLat=Lats_Messung[Max_Index];
    double MinLat=Lats_Messung[Min_Index];
    int Breitenzahl=Min_Index-Max_Index+1;
    double Gitterkonstante=(MaxLat-MinLat)/((double) Breitenzahl);

    //Wir brauchen noch die Höhen, die wir aber kennen
    // und statisch deklarieren
    double mittlere_Hoehe[11] =
      { 74, 77, 80, 83, 86, 89, 92, 106.75, 135, 200, 375};
    double untere_Hoehe[11] =
      { 72.5, 75.5, 78.5, 81.5, 84.5, 87.5, 90.5, 93.5, 120, 150, 250 };
    double obere_Hoehe[11] =
      { 75.5, 78.5, 81.5, 84.5, 87.5, 90.5, 93.5, 120, 150, 250, 500 };
    int Anzahl_Hoehen=11;
    //Gitterpunkte des Gitters erzeugen sortiert nach Höhen und Breiten
    for(int i=0;i<Breitenzahl;i++) //Breiten von Nord nach Süd
    {
        for(int j=0;j<Anzahl_Hoehen;j++) //Höhen von unten nach oben
        {
            Gitterpunkt GP;
            //Index des Punktes
            GP.m_eigener_Index=i+j*Breitenzahl;
            // Index der Nachbarpunkte...falls Nachbar nicht existiert..-> -1
            if((i-1)>=0)
                {GP.m_Index_Nord_Nachbar=(i-1)+j*Breitenzahl;}         //Nord
            else
                {GP.m_Index_Nord_Nachbar=-1;}
            if((i+1)<=Breitenzahl)
                {GP.m_Index_Sued_Nachbar=(i+1)+j*Breitenzahl;}       //Sued
            else
                {GP.m_Index_Sued_Nachbar=-1;}
            if((j-1)>=0)
                {GP.m_Index_unterer_Nachbar=i+(j-1)*Breitenzahl;}      //unten
            else
                {GP.m_Index_unterer_Nachbar=-1;}
            if((j+1)<=Anzahl_Hoehen)
                {GP.m_Index_oberer_Nachbar=i+(j+1)*Breitenzahl;}      //oben
            else
                {GP.m_Index_oberer_Nachbar=-1;}
            // und die Diagonalen Punkte
            if(((i-1)>=0) && ((j+1)<=Anzahl_Hoehen))
                {GP.m_Index_oberer_Nord_Nachbar=(i-1)+(j+1)*Breitenzahl;}         //Nord oben
            else
                {GP.m_Index_oberer_Nord_Nachbar=-1;}
            if(((i-1)>=0) && ((j-1)>=0))
                {GP.m_Index_unterer_Nord_Nachbar=(i-1)+(j-1)*Breitenzahl;}         //Nord unten
            else
                {GP.m_Index_unterer_Nord_Nachbar=-1;}
            if(((i+1)<=Breitenzahl) && ((j+1)<=Anzahl_Hoehen))
                {GP.m_Index_oberer_Sued_Nachbar=(i+1)+(j+1)*Breitenzahl;}         //Süd oben
            else
                {GP.m_Index_oberer_Sued_Nachbar=-1;}
            if(((i+1)<=Breitenzahl) && ((j-1)>=0))
                {GP.m_Index_oberer_Sued_Nachbar=(i+1)+(j-1)*Breitenzahl;}         //Süd unten
            else
                {GP.m_Index_oberer_Sued_Nachbar=-1;}

            //  3 Hoehen
            GP.m_Hoehe=mittlere_Hoehe[j];
            GP.m_Max_Hoehe=obere_Hoehe[j];
            GP.m_Min_Hoehe=untere_Hoehe[j];
            //  3 Breiten
            GP.m_Breite=MaxLat+i*Gitterkonstante;
            GP.m_Max_Breite=GP.m_Breite+0.5*Gitterkonstante;
            GP.m_Min_Breite=GP.m_Breite-0.5*Gitterkonstante;
            //SZA initialisieren
            GP.m_SZA=0;
            // Am Ende Gitterpunkt übergeben
            this->m_Gitter[GP.m_eigener_Index]=GP;
        }// for j ende
    }// for i ende
    // Nun noch die Nord-Sued-Richtung der Messung bestimmen
    for(unsigned int i=0;i<AM_Limb.size();i++)//Alle Messungen durchgehen
    {
        for(unsigned int j=0;j<Lats_Messung.size();j++) // Alle Breiten durchgehn
        {
            if(AM_Limb[i].m_Lat_TP==Lats_Messung[j])
                // Lat finden....(hoffen das die lats nicht zufällig gleich sind,
                // was echt krasser Zufall wär,
                // da auf lat auf 4 Nachkommastellen genau angegeben ist)
            {
                if((j<=Max_Index) || (j>Min_Index))
                {
                    AM_Limb[i].Nord_Sued_Messung=false; // dann haben wir eine Süd-Nordmessung
                }
                else
                {
                    //Punkt liegt zwischen Minimum und Maximum also auf Nord-Süd-Richtung
                    AM_Limb[i].Nord_Sued_Messung=true;
                }
                break; //Schleife beenden
            }
        }//ende for j
    }//ende for i
}// Ende Retrievalgitter_erzeugen_und_Messrichtung_herausfinden */

// globale Hilfsfunktionen
int Get_Index_of_Maximum(vector<double> A)
{
	vector<double>::iterator it;
	it = max_element(A.begin(), A.end());

	return distance(A.begin(), it);
}//Ende Get_Index_of_Maximum
int Get_Index_of_Minimum(vector<double> A, int Startindex)
{
	//so kommt Minimum nach maximum...alles andere ist komisch...
	//TODO so schreiben, das nachricht kommt, falls das absolute minimum einen
	//kleineren index als das den Startindex hat
	vector<double>::iterator it;
	it = min_element(A.begin() + Startindex, A.end());

	return distance(A.begin(), it);
}//Ende Get Index of Minimum
