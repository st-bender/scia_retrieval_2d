/*
 * Sonnenspektrum.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 20.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
#include "Sonnenspektrum.h"
#include <string>
#include<fstream>
#include <sstream>
#include<iostream>
#include <iterator>
#include <algorithm>

#include "Messung_Limb.h"
#include "Glaetten.h"

using std::cout;
using std::string;

//=============================================================================
int Sonnenspektrum::Laden_GOME(string Dateiname, string Fallback_Dateiname)
{
	double dummy;
	// Es werden die korrigierten Sonnenspektren von Mark Weber verwendet;
	// die liegen täglich von 2002 bis 2006 vor
	// 5 Spalten
	// erste WL, 2te Intensität
	// Es muss nicht das volle Spektrum geladen werden,
	// die ersten 800 zeilen decken den 200nm-300nm Bereich ab

	//cout<<"Aufruf von Sonnenspektrum laden\n";
	std::ifstream infile(Dateiname.c_str());
	if (!(infile.is_open())) {
		cout << "Achtung!!!: verwende Fallbackspektrum\n";
		infile.open(Fallback_Dateiname.c_str());
	}
	if (!(infile.is_open())) {
		cout << "Achtung!!!: kein Fallbackspektrum\n";
		cout << "kritischer Fehler!!!\n";
		return 1;
	}
	m_Anzahl_WL = 850;
	for (int i = 0; i < m_Anzahl_WL; i++) {
		// etwas mehr Punkte einlesen,
		// damit nacher auf 826 Punkte des Erdscheinspektrums
		// Interpoliert werden kann
		infile >> m_Wellenlaengen[i]
			   >> m_Intensitaeten[i]
			   >> dummy
			   >> dummy
			   >> dummy;
	}
	// ok, einlesen geht also
	//cout<<m_Wellenlaengen[25]<<"\t"<<m_Intensitaeten[25];
	return 0;
}//Sonnenspektrum::Laden(string Dateiname, string Fallback_Dateiname) ende
//=============================================================================
//=============================================================================
////////////////////////////////////////////////////
// Methodenstart  Laden_Scia
/////////////////////////////////////////////////////
int Sonnenspektrum::Laden_SCIA(string Dateiname, string Fallback_Dateiname)
{
	// Lädt das Sonnenspektrum das mit SCIAMACHY während des gleichen Orbits
	// gemessen wird
	int nh = 5, columns = 2;
	double read_wl, read_int, dummy;
	string s_dummy;
	std::stringstream ss;
	std::ifstream infile(Dateiname.c_str());
	if (!(infile.is_open())) {
		cout << "Achtung!!!: verwende Fallbackspektrum\n";
		infile.open(Fallback_Dateiname.c_str());
	}
	if (!(infile.is_open())) {
		cout << "Achtung!!!: kein Fallbackspektrum\n";
		cout << "kritischer Fehler!!!\n";
		return 1;
	}
	// check first line for number of header lines
	getline(infile, s_dummy);
	if (s_dummy[0] != '#') {
		ss << s_dummy;
		ss >> nh;
	}
	// nh dummyzeilen
	for (int i = 0; i < nh; i++) {
		getline(infile, s_dummy);
		if (s_dummy.find("accuracy") != std::string::npos)
			columns = 3;
	}

	infile >> m_Anzahl_WL;
	if (m_Anzahl_WL == 0) {
		std::cerr << "Sonnenspektrum ist unbrauchbar!" << std::endl;
		return 2;
	}

	getline(infile, s_dummy); // Rest der Zeile
	if (nh > 5) {
		// "D0" files contain more lines with
		// meta data which we don't need.
		for (int i = 0; i < 3; ++i)
			getline(infile, s_dummy);
	}
	for (int i = 0; i < m_Anzahl_WL; i++) {
		if (nh > 5 && columns < 3)
			// "D0" files contain only tow columns.
			infile >> read_wl >> read_int;
		else
			infile >> read_wl >> read_int >> dummy;
		m_Wellenlaengen.push_back(read_wl);
		m_Intensitaeten.push_back(read_int);
	}

	int_ref1 = get_rad_at_wl(228.6);
	int_ref2 = get_rad_at_wl(240.2);
	return 0;
}
////////////////////////////////////////////////////
// ENDE  Laden_Scia
/////////////////////////////////////////////////////

void Sonnenspektrum::moving_average(int window_size)
{
	m_Intensitaeten = my_moving_average(m_Intensitaeten, window_size);
}
void Sonnenspektrum::savitzky_golay(int window_size)
{
	m_Intensitaeten = my_savitzky_golay(m_Intensitaeten, window_size);
}
void Sonnenspektrum::saoref_to_sciamachy()
{
	m_Intensitaeten = my_sciamachy_blur(m_Intensitaeten);
}
double Sonnenspektrum::get_rad_at_wl(double wl)
{
	std::vector<double>::iterator wlit =
		std::lower_bound(m_Wellenlaengen.begin(),
			m_Wellenlaengen.end(), wl);
	std::ptrdiff_t i = std::distance(m_Wellenlaengen.begin(), wlit);
	double rad = interpolate(m_Wellenlaengen, m_Intensitaeten, wl);
	std::cerr << "wls = " << m_Wellenlaengen.at(i - 1)
		<< ", " << m_Wellenlaengen.at(i) << std::endl;
	std::cerr << "int0 = " << m_Intensitaeten.at(i - 1)
		<< ", " << m_Intensitaeten.at(i)
		<< " -> " << rad << std::endl;

	return rad;
}

void Sonnenspektrum::Interpolieren(Messung_Limb &Messung_Erdschein)
{
	double I1, I2;
	//Beide Spektren haben nur einen offset, aber bei fast allen Messungen den
	//selben Der 16+1te Punkt des Sonnenspektrums entspircht in etwa dem 1
	//Punkt des Limbspektrums (Bei GOME)
	//
	//Dies ist ein guter Startwert für den Linken und den Rechten
	//Interpolationspunkt, bei bedarf können dann noch bis zu 5 Schritte in
	//beide Richtungen gemacht werden...falls das mal nicht so geht, so
	//Programmabbruch
	// Damit müssen in fast allen Fällen nur 2 bis 3 Schritte gemacht werden,
	// was deutlich schneller geht, als eine allgemeine Interpolation ala
	// quicksearch (auch wenn die vermutlich bei 800 Messpunkte auch nur ca. 10
	// Schritte braucht)

	// die sciaspektren müssen im Prinzip nicht interpoliert werden, da die
	// gleich sind...aber trotzdem( geht schnell, kein Risiko)

	double kleine_WL, grosse_WL;
	double int_interp, wl_interp;
	long Index_kleine_WL = 0, Index_grosse_WL = 1;
	std::vector<double>::iterator low, me_wl_it;
	// cout<<Messung_Erdschein.m_Wellenlaengen[0]<<"\n";
	//cout<<Messung_Erdschein.m_Intensitaeten[0]<<"\n";
	m_Int_interpoliert.resize(0);
	m_WL_interpoliert.resize(0);
	for (me_wl_it = Messung_Erdschein.m_Wellenlaengen.begin();
			me_wl_it != Messung_Erdschein.m_Wellenlaengen.end();
			++me_wl_it) {
		//Für alle Punkte des neuen Sonnenspektrums
		// Startpunkt
		low = lower_bound(m_Wellenlaengen.begin(), m_Wellenlaengen.end(),
				(*me_wl_it));
		if (low == m_Wellenlaengen.begin()) {
			Index_kleine_WL = 0;
			Index_grosse_WL = 1;
			I1 = 1.;
			I2 = 0.;
		} else {
			if (low == m_Wellenlaengen.end()) --low; // use the last one
			Index_kleine_WL = distance(m_Wellenlaengen.begin(), low) - 1;
			Index_grosse_WL = Index_kleine_WL + 1;
			kleine_WL = m_Wellenlaengen[Index_kleine_WL];
			grosse_WL = m_Wellenlaengen[Index_grosse_WL];
			I2 = ((*me_wl_it) - kleine_WL) / (grosse_WL - kleine_WL);
			I1 = 1.0 - I2;
		}
		// Einfach altes Fenster überschreiben, das wird nicht
		int_interp = I1 * m_Intensitaeten[Index_kleine_WL]
			+ I2 * m_Intensitaeten[Index_grosse_WL];
		wl_interp = (*me_wl_it);
		m_Int_interpoliert.push_back(int_interp);
		m_WL_interpoliert.push_back(wl_interp);
	}
	//cout<<Messung_Erdschein.m_Wellenlaengen[0]<<"\n";
}//Interpolieren(Messung_Limb Messung_Erdschein) ende
//=============================================================================
///////////////////////////////////////////////////
// Methodenstart nicht interpolieren
///////////////////////////////////////////////////
void Sonnenspektrum::nicht_interpolieren()
{
	// vector copies
	m_Int_interpoliert = m_Intensitaeten;
	m_WL_interpoliert = m_Wellenlaengen;
}
///////////////////////////////////////////////////
// ENDE nicht interpolieren
///////////////////////////////////////////////////




int Sonnenspektrum::Speichern(string Dateiname)  //zur Kontrolle
{
	std::ofstream outfile(Dateiname.c_str());
	if (!(outfile.is_open())) {
		cout << Dateiname << " konnte nicht zum Speichern geöffnet werden\n";
		return 1;
	}
	outfile << "WL\t\t" << "I" << "\n";
	for (unsigned int i = 0; i < m_WL_interpoliert.size(); i++) {
		outfile << m_WL_interpoliert[i] << "\t" << m_Int_interpoliert[i] << "\n";
	}
	return 0;
} //Ende Speichern

int Sonnenspektrum::Speichern_was_geladen_wurde(string Dateiname)
{
	std::ofstream outfile(Dateiname.c_str());
	if (!(outfile.is_open())) {
		cout << Dateiname << " konnte nicht zum Speichern geöffnet werden\n";
		return 1;
	}
	outfile << "WL\t\t" << "I\n";
	for (int i = 0; i < m_Anzahl_WL; i++) {
		outfile << this->m_Wellenlaengen[i] << "\t"
				<< this->m_Intensitaeten[i] << "\n";
	}
	return 0;
}//Ende Speichern_was_geladen_wurde
