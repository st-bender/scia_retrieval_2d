/*
 * Liniendaten.cpp
 *
 *  Created on: 13.04.2010
 *      Author: martin
 */

// Methoden für die Klasse Liniendaten

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstdio>

#include "Liniendaten.h"

using namespace std;

Liniendaten::Liniendaten()
{
	//Alles mit 0 initialisieren
	m_Wellenlaenge = 0;
	m_rel_Einstein = 0;
	m_f_Wert = 0;
	m_E1 = 0;
	m_E2 = 0;
	m_Gamma = 0;
}


/*************************************************************************/
void Liniendaten::Einlesen(string Dateiname, double Wellenlaenge)
{
	ifstream infile;
	string zeile;
	//Datei öffnen
	infile.open(Dateiname.c_str());
	if (!infile.is_open()) {
		cerr << "Datei " << Dateiname << " fehlt...Programm stürzt ab\n";
		exit(1);
	}
	// erste 8 Zeilen ignorieren
	for (int i = 0; i < 8; i++) {
		getline(infile, zeile);
	}
	while (!(infile.eof())) {
		double current_WL;
		infile >> current_WL;
		// cout<<"current_WL: "<<current_WL<<"\t"<<"Wellenlaenge:"<<Wellenlaenge
		// <<"\n";
		if (current_WL == Wellenlaenge) {
			this->m_Wellenlaenge = Wellenlaenge;
			infile >> this->m_rel_Einstein;
			infile >> this->m_f_Wert;
			infile >> this->m_E1;
			infile >> this->m_E2;
			//this->Auf_Bildschirm_Ausgeben();
			break;
		} else {
			getline(infile, zeile);
		}
	}//ende while
	//Datei schließen
	infile.close();

}
/* Einlesen ende***************************************************************/
/******************************************************************************/
/* Streuwinkel_ermitteln Ende**************************************************/
/******************************************************************************/
void Liniendaten::Emissivitaet_ermitteln()
{
	//!!!! ACHTUNG HIER GIBTS NOCH DISKUSSIONSBEDARF!!!!!
	//
	// In  David Cleary, " Daytime High Lattitude Rocket Observations of the NO
	// gamma, delta, epsilon bands, 1986
	//
	// Wird angeführt, dass zuerst über die gemessene Intensität Integriert
	// werden muss.  Dann gibt die Gesamtintensität der Linie in 1/m^2s gleich
	// Der wellenlängenspezifische Sonnenintensität in 1/m^2s nm * r_elektron[m]
	// mal lamda^2[nm^2] Das ist zwar Konsistent aber nicht so ganz schlüssig...
	// eigentlich sollte F die selbe Einheit wie I haben und der Rest des
	// Gammafaktors sollte der Wirkungsquerschnitt für Resonanzemission in m^2
	//

	//Konstanten
	const double pi = 3.1415926535898; //const pi=3.1415926535897932384626433;
	const double r_elektron_klassisch = 2.8179402894e-6; // in nm
	//Teilformeln
	//double Phasenfunktion=0.75*m_E1*(cos(pi/180*m_theta)*cos(pi/180*m_theta)+1)+m_E2;
	//Phasenfunktion kommt später dazu, da Streuwinkel noch nicht bekannt
	// diese Korrektur ist aber in etwa ein Faktor 1
	//cout<<"Phasenfunktion: "<<Phasenfunktion<<"\n";
	double Volumenfaktor =
		pi * r_elektron_klassisch * m_f_Wert * (m_Wellenlaenge * 1E-7)
		* (m_Wellenlaenge * 1E-7); // in cm^2nm also wellenlängen in cm r in nm
	//hier wirds mysteriöser 1E7 =100nm ist in etwa wieder
	//....falls das in cm ist ist das eher 1e9...sollte eigentlich
	//die Intervalläge sein //des F Intervalls
	//Der Wert in m^2 ist zumindest der gleiche wie Marcos für cm^2,
	//da warns 1e-7 anderer_Faktor /=1E-7;
	// Emmisivität ausrechnen
	m_Gamma =/*Phasenfunktion* */Volumenfaktor * m_rel_Einstein;
	//-> ist zumindest konsistent
}
/*Emissivitaet_ermitteln ende**************************************************/
/******************************************************************************/
void Liniendaten::Auf_Bildschirm_Ausgeben()
{
	cout << "Wellenlaenge:          " << m_Wellenlaenge << "\n"; //in nm
	cout << "rel Einstein:             " << m_rel_Einstein << "\n";
	cout << "Oszillatorstaerke:     " << m_f_Wert << "\n"; // Oszillatorstärke
	cout << "E1:                         " << m_E1 << "\n";
	cout << "E2:                         " << m_E2 << "\n";
	//cout<<"Streuwinkel in Grad: "<<m_theta<<"\n";
	cout << "Emissivitaet:            " << m_Gamma << "\n";    //Emissivität
}
/*  Auf_Bildschirm_Ausgeben ende***********************************************/
