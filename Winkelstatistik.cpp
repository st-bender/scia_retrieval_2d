/*
 * Winkelstatistik.cpp
 *
 *  Created on: 12.10.2010
 *      Author: martin
 */

#include"Winkelstatistik.h"
#include"MPL_Vektor.h"
#include<string>
#include <cmath>
#include<iostream>
#include<fstream>

using namespace std;
//////////////////////////////////////////////////////
// KONSTRUKTOR
//////////////////////////////////////////////////////
Winkelstatistik::Winkelstatistik()
{
	m_counter_0_001 = 0; //winkel zwischen 0 und 0,001Grad
	m_counter_0_002 = 0;
	m_counter_0_005 = 0;
	m_counter_0_01 = 0;
	m_counter_0_02 = 0;   //Alles bis hier ist OK
	m_counter_0_05 = 0;   //Das geht noch
	m_counter_0_1 = 0;    // Ab hier wird so ziemlich alles falsch
	m_counter_0_2 = 0;
	m_counter_0_5 = 0;
	m_counter_1 = 0;
	m_counter_2 = 0;
	m_counter_5 = 0;
	m_counter_mehr = 0;

}
//////////////////////////////////////////////////////
// ENDE KONSTRUKTOR
//////////////////////////////////////////////////////

//Methoden
//////////////////////////////////////////////////////
// Methodenstart Winkel_berechnen_und_einordnen
//////////////////////////////////////////////////////
int Winkelstatistik::Winkel_berechnen_und_einordnen(
		MPL_Vektor Verbindungsvektor, MPL_Vektor Tangentenpunkt,
		int &Winkel_OK)
{
	const double pi = M_PI;
	// Winkel OK-> 0 ; nicht OK -> 1
	//Rückgabewert der Funktion 0 Winkel Klassifizierbar,
	// 1 nicht (irgendwas ist dann falsch im Programm)
	// Der Tangentenpunkt ist so definiert,
	// dass dort die LOS des Satelliten dort senkrecht zum Zenit steht,
	// Die Abweichungen davon werden hier berechnet
	MPL_Vektor V_Norm(Verbindungsvektor), TP_Norm(Tangentenpunkt);
	V_Norm.Normieren();
	TP_Norm.Normieren();
	double COS_Winkel = TP_Norm * V_Norm;
	double W = 180.0 / pi * acos(COS_Winkel) - 90.0;
	double W_Abs = sqrt(W * W);

	if ((0.0 < W_Abs) && (W_Abs <= 0.001)) {
		m_counter_0_001++;
		Winkel_OK = 0;
		return 0;
	}
	if ((0.001 < W_Abs) && (W_Abs <= 0.002)) {
		m_counter_0_002++;
		Winkel_OK = 0;
		return 0;
	}
	if ((0.002 < W_Abs) && (W_Abs <= 0.005)) {
		m_counter_0_005++;
		Winkel_OK = 0;
		return 0;
	}
	if ((0.005 < W_Abs) && (W_Abs <= 0.01)) {
		m_counter_0_01++;
		Winkel_OK = 0;
		return 0;
	}
	if ((0.01 < W_Abs) && (W_Abs <= 0.02)) {
		m_counter_0_02++;
		Winkel_OK = 0;
		return 0;
	}
	if ((0.02 < W_Abs) && (W_Abs <= 0.05)) {
		m_counter_0_05++;
		Winkel_OK = 1; //Hier setz ich mal die Grenze
		return 0;
	}
	if ((0.05 < W_Abs) && (W_Abs <= 0.1)) {
		m_counter_0_1++;
		Winkel_OK = 1;
		return 0;
	}

	if ((0.1 < W_Abs) && (W_Abs <= 0.2)) {
		m_counter_0_2++;
		Winkel_OK = 1;
		return 0;
	}
	if ((0.2 < W_Abs) && (W_Abs <= 0.5)) {
		m_counter_0_5++;
		Winkel_OK = 1;
		return 0;
	}
	if ((0.5 < W_Abs) && (W_Abs <= 1.0)) {
		m_counter_1++;
		Winkel_OK = 1;
		return 0;
	}
	if ((1.0 < W_Abs) && (W_Abs <= 2.0)) {
		m_counter_2++;
		Winkel_OK = 1;
		return 0;
	}
	if ((2.0 < W_Abs) && (W_Abs <= 5.0)) {
		m_counter_5++;
		Winkel_OK = 1;
		return 0;
	}
	if (5.0 < W_Abs) {
		m_counter_mehr++;
		Winkel_OK = 1;
		return 0;
	}

	return 1; //Merkwürdiges Ende
}
//////////////////////////////////////////////////////
// ENDE Winkel_berechnen_und_einordnen
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
// Methodenstart Statistik_auf_Bildschirm_ausgeben
//////////////////////////////////////////////////////
void Winkelstatistik::Statistik_auf_Bildschirm_ausgeben()
{
	cout << "Winkelstatistik\n";
	cout << "Limbmessungen: " << m_counter_0_001 + m_counter_0_002 +
		 m_counter_0_005 + m_counter_0_01 +
		 m_counter_0_02 + m_counter_0_05 + m_counter_0_1 + m_counter_0_2 +
		 m_counter_0_5 + m_counter_1 +
		 m_counter_2 + m_counter_5 << "\n";
	cout << "0     bis 0.001: " << m_counter_0_001 << "\t" << "Fehler Tangentenpunkt:   60m " << "\n";
	cout << "0.001 bis 0.002: " << m_counter_0_002 << "\t" << "Fehler Tangentenpunkt:  110m " << "\n";
	cout << "0.002 bis 0.005: " << m_counter_0_005 << "\t" << "Fehler Tangentenpunkt:  280m " << "\n";
	cout << "0.005 bis 0.01 : " << m_counter_0_01 << "\t" << "Fehler Tangentenpunkt:  570m " << "\n";
	cout << "0.01  bis 0.02 : " << m_counter_0_02 << "\t" << "Fehler Tangentenpunkt: 1130m " << "\n";
	cout << "0.02  bis 0.05 : " << m_counter_0_05 << "\t" << "Fehler Tangentenpunkt: 2840m " << "\n";
	cout << "0.05  bis 0.1  : " << m_counter_0_1 << "\t" << "Fehler Tangentenpunkt: 5670m " << "\n";
	cout << "0.1   bis 0.2  : " << m_counter_0_2 << "\t" << "Fehler Tangentenpunkt: 11,3km " << "\n";
	cout << "0.2   bis 0.5  : " << m_counter_0_5 << "\t" << "Fehler Tangentenpunkt:  28km " << "\n";
	cout << "0.5   bis 1    : " << m_counter_1 << "\t" << "Fehler Tangentenpunkt:  56km " << "\n";
	cout << "1     bis   2  : " << m_counter_2 << "\t" << "Fehler Tangentenpunkt: 113km " << "\n";
	cout << "2     bis   5  : " << m_counter_5 << "\t" << "Fehler Tangentenpunkt: 284km " << "\n";
	cout << "mehr           : " << m_counter_mehr << "\n";
	cout << "ok             : " << m_counter_0_001 + m_counter_0_002 +
		 m_counter_0_005 + m_counter_0_01 +
		 m_counter_0_02 << "\n";
	cout << "nicht ok       : "
		 << m_counter_0_05 + m_counter_0_1 + m_counter_0_2
		  + m_counter_0_5 + m_counter_1 + m_counter_2 + m_counter_5 << "\n";

	return;
}
//////////////////////////////////////////////////////
// ENDE Statistik_auf_Bildschirm_ausgeben
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
// Methodenstart Statistik_in_Datei_ausgeben
//////////////////////////////////////////////////////
int Winkelstatistik::Statistik_in_Datei_ausgeben(string Dateiname)
{
	ofstream outfile;
	outfile.open(Dateiname.c_str());
	if (!outfile.is_open()) {
		cerr << "Datei " << Dateiname
			 << " kann nicht zum schreiben der Winkelstatistik geöffnet werden."
			 << endl;
		return 1;
	}
	outfile << "Winkelstatistik\n";
	outfile << "Limbmessungen: " << m_counter_0_001 + m_counter_0_002 +
			m_counter_0_005 + m_counter_0_01 +
			m_counter_0_02 + m_counter_0_05 + m_counter_0_1 + m_counter_0_2 +
			m_counter_0_5 + m_counter_1 +
			m_counter_2 + m_counter_5 << "\n";
	outfile << "0____bis_0.001: " << m_counter_0_001 << "\t" << "Fehler Tangentenpunkt:   60m " << "\n";
	outfile << "0.001_bis_0.002: " << m_counter_0_002 << "\t" << "Fehler Tangentenpunkt:  110m " << "\n";
	outfile << "0.002_bis_0.005: " << m_counter_0_005 << "\t" << "Fehler Tangentenpunkt:  280m " << "\n";
	outfile << "0.005_bis_0.01_: " << m_counter_0_01 << "\t" << "Fehler Tangentenpunkt:  570m " << "\n";
	outfile << "0.01__bis_0.02_: " << m_counter_0_02 << "\t" << "Fehler Tangentenpunkt: 1130m " << "\n";
	outfile << "0.02__bis_0.05_: " << m_counter_0_05 << "\t" << "Fehler Tangentenpunkt: 2840m " << "\n";
	outfile << "0.05__bis_0.1__: " << m_counter_0_1 << "\t" << "Fehler Tangentenpunkt: 5670m " << "\n";
	outfile << "0.1___bis_0.2__: " << m_counter_0_2 << "\t" << "Fehler Tangentenpunkt: 11,3km " << "\n";
	outfile << "0.2___bis_0.5__: " << m_counter_0_5 << "\t" << "Fehler Tangentenpunkt:  28km " << "\n";
	outfile << "0.5___bis_1____: " << m_counter_1 << "\t" << "Fehler Tangentenpunkt:  56km " << "\n";
	outfile << "1_____bis_2____: " << m_counter_2 << "\t" << "Fehler Tangentenpunkt: 113km " << "\n";
	outfile << "2_____bis_5____: " << m_counter_5 << "\t" << "Fehler Tangentenpunkt: 284km " << "\n";
	outfile << "mehr__________: " << m_counter_mehr << "\n";
	outfile << "ok__________: " << m_counter_0_001 + m_counter_0_002 +
			m_counter_0_005 + m_counter_0_01 +
			m_counter_0_02 << "\n";
	outfile << "nicht_ok____: "
			<< m_counter_0_05 + m_counter_0_1 + m_counter_0_2
			 + m_counter_0_5 + m_counter_1 + m_counter_2 + m_counter_5 << "\n";
	outfile.close();
	outfile.clear();

	return 0;
}
//////////////////////////////////////////////////////
// ENDE Statistik_in_Datei_ausgeben
//////////////////////////////////////////////////////
