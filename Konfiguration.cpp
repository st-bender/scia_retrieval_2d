/*
 * Konfiguration.cpp
 *
 *  Created on: 06.04.2010
 *      Author: martin
 */
/*********************************************************************
 Diese Datei enthält die Methoden der Klasse Konfiguration
 ********************************************************************/

//eingebundene Header Dateien
#include "Konfiguration.h" // eigener header
#include <fstream>         // ifstream
#include <string>          // string
#include <cstdio>         // cout
#include <iostream>        // cout
#include <cstdlib>        // für atoi

//namespace
using namespace std;

//Protos für Hilfsfunktionen
vector<double> String_mit_double_Zahlen_in_Vektor_schreiben(string Zeile);
vector<int>     String_mit_int_Zahlen_in_Vektor_schreiben(string Zeile);

//Destructor
Konfiguration::~Konfiguration()
{
	if (m_Retrieval_Kovarianzen != NULL) {
		delete[] m_Retrieval_Kovarianzen;
		m_Retrieval_Kovarianzen = NULL;

	}

}

// Konfiguration_einlesen /////////////////////////////////////////////
int Konfiguration::Konfiguration_einlesen()
{
	/***************************************************
	In dieser Funktion werden alle Parameter aus Scia2d.conf
	eingelesen.
	***************************************************/

	//Datei Öffnen
	ifstream infile;
	//cout<<"Datei einlesen\n";
	infile.open("SCIA2D.conf");
	if (!infile.is_open()) {
		cerr << "Konfigurationsdatei SCIA2D.conf fehlt....Programm stürzt ab\n";
		exit(1);
	}
	//cout<<"Datei geöffnet\n";
	//Datei zeilenweise Einlesen und reagieren
	// am besten SCIA2D.conf ausdrucken und verfolgen
	string Zeile;
	while (!(infile.eof())) {
		getline(infile, Zeile);
		// Einträge abarbeiten
		// Directory Structure /////
		if (Zeile == "Number of emitters") {
			//cout<<"Number of emitters\n";
			getline(infile, Zeile);
			this->m_Anzahl_der_Emitter = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "Solar spectra") {
			//cout<<"solar\n";
			getline(infile, Zeile);
			this->m_Pfad_Solar_Spektrum = Zeile;
			continue;
		}
		if (Zeile == "Fallback solar spectrum") {
			//cout<<"fallbacksolar\n";
			getline(infile, Zeile);
			this->m_Pfad_Solar_Fallback_Spektrum = Zeile;
			continue;
		}
		if (Zeile == "Line parameters") {
			//cout<<"Linepars\n";
			getline(infile, Zeile);
			this->m_Pfad_Linienparameter_Metalle = Zeile;
			continue;
		}
		if (Zeile == "Atmospheric structure") {
			//cout<<"Atmos\n";
			getline(infile, Zeile);
			this->m_Pfad_Dichten_der_Atmosphaerengase = Zeile;
			continue;
		}
		if (Zeile == "Cross sections") {
			//cout<<"Xsect\n";
			getline(infile, Zeile);
			this->m_Pfad_Wirkungsquerschnitte_der_Atmosphaerengase = Zeile;
			continue;
		}
		if (Zeile == "Absorption wavelengths") {
			//cout<<"ABS WL\n";
			getline(infile, Zeile);
			// Zeile bei whitespace zeichen teilen
			m_AbsorbtionsWL_der_Atmosphaerengase = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
			continue;
		}
		// Input DATA /////
		if (Zeile == "Scan ID data path") {
			//cout<<"Scan ID\n";
			getline(infile, Zeile);
			this->m_Pfad_Datei_mit_Dateinamen_fuer_Messungen_eines_Orbits = Zeile;
			continue;
		}
		if (Zeile == "Correction factors") {
			//cout<<"Correction factors\n";
			getline(infile, Zeile);
			this->m_Pfad_Korrekturfaktoren = Zeile;
			continue;
		}
		if (Zeile == "MinAlt and MaxAlt") {
			//cout<<"Minalt maxalt\n";
			getline(infile, Zeile);
			vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
			this->m_MinAlt = dummy[0];
			this->m_MaxAlt = dummy[1];
			continue;
		}
		if (Zeile == "Altitude grid extensions") {
			//cout<<"Altitude grid extensions\n";
			getline(infile, Zeile);
			int Zeilenzahl = atoi(Zeile.c_str());
			this->m_Anzahl_zusaetzliche_Hoehengitterpunkte = Zeilenzahl;
			for (int i = 0; i < Zeilenzahl; i++) {
				getline(infile, Zeile);
				vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
				this->m_Grid_ext_low.push_back(dummy[0]);
				this->m_Grid_ext_high.push_back(dummy[1]);
			}
			continue;
		}
		if (Zeile == "TOA") {
			//cout<<"TOA\n";
			getline(infile, Zeile);
			this->m_TOA = atof(Zeile.c_str());
			continue;
		}
		// Selection rules /////
		if (Zeile == "Nadir retrieval") {
			//cout<<"Nadir Retrieval\n";
			getline(infile, Zeile);
			this->m_Nadir_only = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "Night") {
			//cout<<"Night\n";
			getline(infile, Zeile);
			this->m_Nachtmessung = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "Geolocation") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_Geolocation = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "Large SZA") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_Large_SZA = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "NLC") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_NLC = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "Maximal SZA") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_Maximaler_SZA = atof(Zeile.c_str());
			continue;
		}
		if (Zeile == "Geolocation boundaries") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			//cout<<Zeile<<"\n";
			m_Geolocation_Grenzen = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
			continue;
		}
		// Baseline fit parameters /////
		if (Zeile == "Baseline windows") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			int Zeilenzahl = atoi(Zeile.c_str());
			this->m_Anzahl_Baseline_Intervalle = Zeilenzahl;
			for (int i = 0; i < Zeilenzahl; i++) {
				getline(infile, Zeile);
				vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
				this->m_Baselinefenster_WL_low.push_back(dummy[0]);
				this->m_Baselinefenster_WL_high.push_back(dummy[1]);
			}
			continue;
		}
		// Retrieval species fit parameters /////
		if (Zeile == "Retrieval windows") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			int Zeilenzahl = atoi(Zeile.c_str());
			this->m_Anzahl_Retrieval_Intervalle = Zeilenzahl;
			//cout<<this->m_Anzahl_Retrieval_Intervalle<<"\n\n";
			for (int i = 0; i < Zeilenzahl; i++) {
				getline(infile, Zeile);
				vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
				this->m_Retrievalfenster_WL_low.push_back(dummy[0]);
				this->m_Retrievalfenster_WL_high.push_back(dummy[1]);
			}
			continue;
		}
		if (Zeile == "Assignment of wl windows") {
			getline(infile, Zeile);
			//cout<<Zeile<<"\n\n";
			m_Assignment_of_WL_Windows = String_mit_int_Zahlen_in_Vektor_schreiben(Zeile);
			//cout<<Zeile<<"blabla\n\n";
			continue;
		}
		//Regularization /////
		if (Zeile == "Retrieval covariances") {
			//cout<<Zeile<<"\n";

			this->m_Retrieval_Kovarianzen = new double[this->m_Anzahl_der_Emitter * 3];

			for (int i = 0; i < 3; i++) {
				getline(infile, Zeile);
				//D in E umwandeln
				for (uint j = 0; j < Zeile.size(); j++) {
					if (Zeile[j] == 'D')
						Zeile[j] = 'E';
				}
				//cout<<Zeile<<"\n";
				vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
				for (int j = 0; j < m_Anzahl_der_Emitter; j++) {
					m_Retrieval_Kovarianzen[j + i * m_Anzahl_der_Emitter] = dummy[j];
				}
			}
			//	continue;
		}
		// Miscellaneous
		if (Zeile == "Error thresholds") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_Fehlergrenzen = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
			continue;
		}
		if (Zeile == "Spectral FWHM") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_FWHM = atof(Zeile.c_str());
			continue;
		}
		if (Zeile == "Do correction of radiances") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_Do_Corrections_of_Radiances = atoi(Zeile.c_str());
			continue;
		}
		if (Zeile == "OS type") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			this->m_Betriebssystem = Zeile;
			continue;
		}
		if (Zeile == "Maximal number of LM steps and scaling factor") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
			this->m_Max_Zahl_Levenberg_Schritte = dummy[0];
			this->m_Levenberg_Schrittweite = (int) dummy[1];
			continue;
		}
		if (Zeile == "Maximal number of iteration steps and convergence threshold") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			vector<double> dummy = String_mit_double_Zahlen_in_Vektor_schreiben(Zeile);
			this->m_Max_Zahl_Iterationen = (int) dummy[0];
			this->m_Convergence_Treshold = dummy[1];
			continue;
		}
	} //ende while !eof

	//Datei Schließen
	//cout<<"Datei schließen\n";
	infile.close();
	//cout<<"Datei geschlossen\n";

	return 0;
}
// Ende Konfiguration_einlesen ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int Konfiguration::Konfiguration_anzeigen()
{
	// Directory structure /////
	cout << "Anzahl der Emitter: " << this->m_Anzahl_der_Emitter << "\n";
	cout << "Solarpfad: " << this->m_Pfad_Solar_Spektrum << "\n";
	cout << "Ersatz-Solarpfad: " << this->m_Pfad_Solar_Fallback_Spektrum << "\n";
	cout << "Pfad Linienparameter: " << this->m_Pfad_Linienparameter_Metalle << "\n";
	cout << "Pfad Atmosphärendichten: " << this->m_Pfad_Dichten_der_Atmosphaerengase << "\n";
	cout << "Pfad AtmosphärenWQ: " << this->m_Pfad_Wirkungsquerschnitte_der_Atmosphaerengase << "\n";
	for (uint i = 0; i < this->m_AbsorbtionsWL_der_Atmosphaerengase.size(); i++) {
		cout << "Absorptionswellenlänge " << i + 1 << ": " << this->m_AbsorbtionsWL_der_Atmosphaerengase[i] << "\n";
	}
	cout << "Pfad Orbitliste: " << this->m_Pfad_Datei_mit_Dateinamen_fuer_Messungen_eines_Orbits << "\n";
	cout << "Pfad Korrekturfaktoren: " << this->m_Pfad_Korrekturfaktoren << "\n";
	cout << "MinAlt: " << this->m_MinAlt << "\n";
	cout << "MaxAlt: " << this->m_MaxAlt << "\n";
	cout << "Zahl der Gittererweiterungen: " << this->m_Anzahl_zusaetzliche_Hoehengitterpunkte << "\n";
	for (int i = 0; i < this->m_Anzahl_zusaetzliche_Hoehengitterpunkte; i++) {
		cout << "Erweiterung: " << i + 1 << " von " << this->m_Grid_ext_low[i] << " bis " << this->m_Grid_ext_high[i] << "\n";
	}
	cout << "TOA: " << this->m_TOA << "\n";
	cout << "Nadir Only: " << this->m_Nadir_only << "\n";
	cout << "Nachtmessung: " << this->m_Nachtmessung << "\n";
	cout << "Geolocation: " << this->m_Geolocation << "\n";
	cout << "Large SZA: " << this->m_Large_SZA << "\n";
	cout << "Mesosphärische Wolken: " << this->m_NLC << "\n";
	cout << "Max SZA: " << this->m_Maximaler_SZA << "\n";
	cout << "Geolocation Boundaries: ";
	for (int i = 0; i < 4; i++) {
		cout << this->m_Geolocation_Grenzen[i] << "\t";
	}
	cout << "\n";
	cout << "Anzahl Baseline Fenster:" << this->m_Anzahl_Baseline_Intervalle << "\n";
	for (int i = 0; i < m_Anzahl_Baseline_Intervalle; i++) {
		cout << "BL Intervall "
			 << i
			 << ": "
			 << this->m_Baselinefenster_WL_low[i]
			 << " " << this->m_Baselinefenster_WL_high[i]
			 << "\n";
	}
	cout << "Anzahl Retrieval Fenster:" << this->m_Anzahl_Retrieval_Intervalle << "\n";
	for (int i = 0; i < this->m_Anzahl_Retrieval_Intervalle; i++) {
		cout << "Retrieval Intervall " << i << ": " << this->m_Retrievalfenster_WL_low[i] << " " << this->m_Retrievalfenster_WL_high[i] << "\n";
	}
	cout << "Assignment of WL:";
	//cout<<this->m_Assignment_of_WL_Windows.size();
	for (uint i = 0; i < this->m_Assignment_of_WL_Windows.size(); i++) {
		cout << m_Assignment_of_WL_Windows[i] << "\t";
	}
	cout << "\n";
	cout << "S_apriori: ";
	for (int i = 0; i < this->m_Anzahl_der_Emitter; i++) {
		cout << m_Retrieval_Kovarianzen[i + 0 * m_Anzahl_der_Emitter] << "\t";

	}
	cout << "\n";
	cout << "S_lat: ";
	for (int i = 0; i < this->m_Anzahl_der_Emitter; i++) {
		cout << m_Retrieval_Kovarianzen[i + 1 * m_Anzahl_der_Emitter] << "\t";
	}
	cout << "\n";
	cout << "S_lon: ";
	for (int i = 0; i < this->m_Anzahl_der_Emitter; i++) {
		cout << m_Retrieval_Kovarianzen[i + 2 * m_Anzahl_der_Emitter] << "\t";
	}
	cout << "\n";
	cout << "error tresholds: ";
	for (uint i = 0; i < this->m_Fehlergrenzen.size(); i++) {
		cout << this->m_Fehlergrenzen[i] << "\t";
	}
	cout << "\n";
	cout << "FWHM: " << this->m_FWHM << "\n";
	cout << "Do Corrections of Radiances? " << this->m_Do_Corrections_of_Radiances << "\n";
	cout << "OS: " << this->m_Betriebssystem << "\n";
	cout << "Max LM Steps: " << this->m_Max_Zahl_Levenberg_Schritte << "\n";
	cout << "LM Schrittweite: " << this->m_Levenberg_Schrittweite << "\n";
	cout << "Max Iterations: " << this->m_Max_Zahl_Iterationen << "\n";
	cout << "Convergence treshold: " << this->m_Convergence_Treshold << "\n";
	cout << "\n";
	return 0;
}//ende Konfiguration_anzeigen

// HILFSFUNTIONEN /////////////////////////////////////////////////////
vector<double> String_mit_double_Zahlen_in_Vektor_schreiben(string Zeile)
{
	vector<double> Zahlen;
	int IndexWortanfang = 0;
	int IndexWortende = 0;
	for (uint i = 0; i < Zeile.size(); i++) {
		if (Zeile[i] == ' ') {
			IndexWortende = i;
			string Wort = Zeile.substr(IndexWortanfang, IndexWortende - IndexWortanfang + 1);
			if (IndexWortende > IndexWortanfang + 1) {
				double Zahl = atof(Wort.c_str());
				Zahlen.push_back(Zahl);
				//	cout<<IndexWortanfang<<" "<<IndexWortende<<" "<<Zahl<<"\n";
			}
			IndexWortanfang = IndexWortende;
		}
	}
	//letzte Zahl
	if (Zeile[Zeile.size() - 1] != ' ') {
		IndexWortende = Zeile.size() - 1;
		string Wort = Zeile.substr(IndexWortanfang, IndexWortende - IndexWortanfang + 1);
		if (IndexWortende > IndexWortanfang + 1) {
			double Zahl = atof(Wort.c_str());
			Zahlen.push_back(Zahl);
			//	cout<<IndexWortanfang<<" "<<IndexWortende<<" "<<Zahl<<"\n";
		}
		//IndexWortanfang=IndexWortende;  // wird nicht gebraucht...letzte Zeile
	}
	return Zahlen;
}//ende String_mit_double_Zahlen_in_Vektor_schreiben

////////////////////////////////////////////////

vector<int> String_mit_int_Zahlen_in_Vektor_schreiben(string Zeile)
{
	vector<int> Zahlen;
	int IndexWortanfang = 0;
	int IndexWortende = 0;
	//cout<<"blalbla\n";
	for (uint i = 0; i < Zeile.size(); i++) {
		if (Zeile[i] == ' ') {
			IndexWortende = i;
			string Wort = Zeile.substr(IndexWortanfang, IndexWortende - IndexWortanfang + 1);
			//cout<<Wort<<"\n";
			int Buchstaben = 0; //nichtleerzeichen im String zählen
			for (uint j = 0; j < Wort.size(); j++) {
				if (Wort[j] != ' ')
					Buchstaben += 1;
			}
			if (Buchstaben != 0) {
				int Zahl = atoi(Wort.c_str());
				//	cout<<Zahl<<"\n\n";
				Zahlen.push_back(Zahl);
			}
			IndexWortanfang = IndexWortende;
		}
	}
	//letzte Zahl
	if (Zeile[Zeile.size() - 1] != ' ') {
		IndexWortende = Zeile.size() - 1;
		string Wort = Zeile.substr(IndexWortanfang, IndexWortende - IndexWortanfang + 1);
		int Buchstaben = 0;
		for (uint j = 0; j < Wort.size(); j++) {
			if (Wort[j] != ' ')
				Buchstaben += 1;
		}
		if (Buchstaben != 0) {
			int Zahl = atoi(Wort.c_str());
			Zahlen.push_back(Zahl);
		}

	}
	return Zahlen;
}//ende String_mit_int_Zahlen_in_Vektor_schreiben
// HILFSFUNTIONEN ende /////////////////////////////////////////////////////


