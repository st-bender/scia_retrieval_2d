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
#include <sstream>
#include <cstdlib>        // für atoi

//namespace
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;

//Protos für Hilfsfunktionen
template <class T> vector<T> string_to_vector(string zeile);

// default constructor to initialise the config default values
Konfiguration::Konfiguration() :
	m_Pfad_Solar_Correction_Factors("DATA/sol_corrfac_sao-scia.dat"),
	m_Pfad_Linienparameter_Metalle("DATA/LinePars.dat"),
	m_Pfad_Dichten_der_Atmosphaerengase("DATA/ACE_air_o3_50_200.dat"),
	m_Pfad_Wirkungsquerschnitte_der_Atmosphaerengase("DATA/XSECT_air_o3.dat"),
	m_Pfad_NO_parameters("DATA/Luqueetal.dat"),
	m_Pfad_Ap_index("DATA/spidr_ap_2000-2012.dat"),
	m_Pfad_Kp_index("DATA/spidr_kp_2000-2012.dat"),
	m_Pfad_f107_index("DATA/spidr_f107_2000-2012.dat"),
	m_MinLat(-90.0), m_MaxLat(90.0), m_NLat(72),
	m_TOA(200.0), m_BOA(50.0),
	skip_SAA(true), SAA_cutoff(8.8e10),
	atmo_Temp(200.), NO_apriori(false),
	NO_apriori_bottom(40.0), NO_apriori_top(160.0),
	NO_apriori_scale(1.0), NO_apriori_smoothness(4.0),
	retrieval_algo(1)
{
}

// Konfiguration_einlesen /////////////////////////////////////////////
int Konfiguration::Konfiguration_einlesen(std::string file)
{
	/***************************************************
	In dieser Funktion werden alle Parameter aus Scia2d.conf
	eingelesen.
	***************************************************/

	//Datei Öffnen
	std::ifstream infile(file.c_str());
	//cout<<"Datei einlesen\n";
	if (!infile.is_open()) {
		std::cerr << "Konfigurationsdatei " << file.c_str()
			 << " fehlt....Programm stürzt ab\n";
		exit(1);
	}
	//cout<<"Datei geöffnet\n";
	//Datei zeilenweise Einlesen und reagieren
	// am besten SCIA2D.conf ausdrucken und verfolgen
	string Zeile;
	while (!(infile.eof())) {
		stringstream ss;
		int Zeilenzahl;
		getline(infile, Zeile);
		// Einträge abarbeiten
		// Directory Structure /////
		if (Zeile == "Number of emitters") {
			//cout<<"Number of emitters\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> m_Anzahl_der_Emitter;
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
		if (Zeile == "solar correction factors") {
			getline(infile, Zeile);
			this->m_Pfad_Solar_Correction_Factors = Zeile;
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
		if (Zeile == "NO parameters") {
			getline(infile, Zeile);
			this->m_Pfad_NO_parameters = Zeile;
			continue;
		}
		if (Zeile == "Ap index file") {
			getline(infile, Zeile);
			this->m_Pfad_Ap_index = Zeile;
			continue;
		}
		if (Zeile == "Kp index file") {
			getline(infile, Zeile);
			this->m_Pfad_Kp_index = Zeile;
			continue;
		}
		if (Zeile == "f10.7 index file") {
			getline(infile, Zeile);
			this->m_Pfad_f107_index = Zeile;
			continue;
		}
		if (Zeile == "Absorption wavelengths") {
			//cout<<"ABS WL\n";
			getline(infile, Zeile);
			// Zeile bei whitespace zeichen teilen
			m_AbsorptionsWL_der_Atmosphaerengase = string_to_vector<double>(Zeile);
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
			vector<double> dummy = string_to_vector<double>(Zeile);
			this->m_MinAlt = dummy[0];
			this->m_MaxAlt = dummy[1];
			continue;
		}
		if (Zeile == "MinLat and MaxLat") {
			getline(infile, Zeile);
			vector<double> dummy = string_to_vector<double>(Zeile);
			this->m_MinLat = dummy[0];
			this->m_MaxLat = dummy[1];
			continue;
		}
		if (Zeile == "NLat") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_NLat;
			continue;
		}
		if (Zeile == "Altitude grid extensions") {
			//cout<<"Altitude grid extensions\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> Zeilenzahl;
			this->m_Anzahl_zusaetzliche_Hoehengitterpunkte = Zeilenzahl;
			for (int i = 0; i < Zeilenzahl; i++) {
				getline(infile, Zeile);
				vector<double> dummy = string_to_vector<double>(Zeile);
				this->m_Grid_ext_low.push_back(dummy[0]);
				this->m_Grid_ext_high.push_back(dummy[1]);
			}
			continue;
		}
		if (Zeile == "TOA") {
			//cout<<"TOA\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> m_TOA;
			continue;
		}
		if (Zeile == "BOA") {
			//cout<<"TOA\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> m_BOA;
			continue;
		}
		// Selection rules /////
		if (Zeile == "Nadir retrieval") {
			//cout<<"Nadir Retrieval\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_Nadir_only;
			continue;
		}
		if (Zeile == "Night") {
			//cout<<"Night\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_Nachtmessung;
			continue;
		}
		if (Zeile == "Geolocation") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_Geolocation;
			continue;
		}
		if (Zeile == "Large SZA") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_Large_SZA;
			continue;
		}
		if (Zeile == "NLC") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_NLC;
			continue;
		}
		if (Zeile == "Skip SAA") {
			int saa = 0;
			getline(infile, Zeile);
			ss << Zeile;
			ss >> saa;
			if (saa == 0) skip_SAA = false;
			else skip_SAA = true;
			continue;
		}
		if (Zeile == "SAA cut-off value") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> SAA_cutoff;
			continue;
		}
		if (Zeile == "Maximal SZA") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_Maximaler_SZA;
			continue;
		}
		if (Zeile == "Geolocation boundaries") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			//cout<<Zeile<<"\n";
			m_Geolocation_Grenzen = string_to_vector<double>(Zeile);
			continue;
		}
		// Baseline fit parameters /////
		if (Zeile == "Baseline windows") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> Zeilenzahl;
			this->m_Anzahl_Baseline_Intervalle = Zeilenzahl;
			for (int i = 0; i < Zeilenzahl; i++) {
				getline(infile, Zeile);
				vector<double> dummy = string_to_vector<double>(Zeile);
				this->m_Baselinefenster_WL_low.push_back(dummy[0]);
				this->m_Baselinefenster_WL_high.push_back(dummy[1]);
			}
			continue;
		}
		// Retrieval species fit parameters /////
		if (Zeile == "Retrieval windows") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> Zeilenzahl;
			this->m_Anzahl_Retrieval_Intervalle = Zeilenzahl;
			//cout<<this->m_Anzahl_Retrieval_Intervalle<<"\n\n";
			for (int i = 0; i < Zeilenzahl; i++) {
				getline(infile, Zeile);
				vector<double> dummy = string_to_vector<double>(Zeile);
				this->m_Retrievalfenster_WL_low.push_back(dummy[0]);
				this->m_Retrievalfenster_WL_high.push_back(dummy[1]);
			}
			continue;
		}
		if (Zeile == "Assignment of wl windows") {
			getline(infile, Zeile);
			//cout<<Zeile<<"\n\n";
			m_Assignment_of_WL_Windows = string_to_vector<int>(Zeile);
			//cout<<Zeile<<"blabla\n\n";
			continue;
		}
		//Regularization /////
		if (Zeile == "Retrieval covariances") {
			//cout<<Zeile<<"\n";

			this->m_Retrieval_Kovarianzen.resize(this->m_Anzahl_der_Emitter * 3);

			for (int i = 0; i < 3; i++) {
				getline(infile, Zeile);
				//D in E umwandeln
				string::size_type pos = 0;
				while ((pos = Zeile.find("D", pos)) != string::npos) {
					Zeile.replace(pos, 1, "E");
					pos++;
				}

				//cout<<Zeile<<"\n";
				vector<double> dummy = string_to_vector<double>(Zeile);
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
			this->m_Fehlergrenzen = string_to_vector<double>(Zeile);
			continue;
		}
		if (Zeile == "Spectral FWHM") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_FWHM;
			continue;
		}
		if (Zeile == "Do correction of radiances") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			ss << Zeile;
			ss >> this->m_Do_Corrections_of_Radiances;
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
			vector<double> dummy = string_to_vector<double>(Zeile);
			this->m_Max_Zahl_Levenberg_Schritte = static_cast<int>(dummy[0]);
			this->m_Levenberg_Schrittweite = dummy[1];
			continue;
		}
		if (Zeile == "Maximal number of iteration steps and convergence threshold") {
			//cout<<Zeile<<"\n";
			getline(infile, Zeile);
			vector<double> dummy = string_to_vector<double>(Zeile);
			this->m_Max_Zahl_Iterationen = (int) dummy[0];
			this->m_Convergence_Treshold = dummy[1];
			cout << "max no. of iterations = " << m_Max_Zahl_Iterationen << endl;
			cout << "convergence threshold = " << m_Convergence_Treshold << endl;
			continue;
		}
		if (Zeile == "atmosphere temperature") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> atmo_Temp;
			continue;
		}
		if (Zeile == "number of NO transitions") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> no_NO_transitions;
			continue;
		}
		if (Zeile == "NO transitions") {
			for (unsigned i = 0; i < no_NO_transitions; i++) {
				getline(infile, Zeile);
				std::vector<int> NO_data = string_to_vector<int>(Zeile);
				NO_v_u.push_back(NO_data.at(0));
				NO_v_l.push_back(NO_data.at(1));
				NO_v_l_abs.push_back(NO_data.at(2));
			}
			continue;
		}
		if (Zeile == "NO apriori") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> NO_apriori;
			continue;
		}
		if (Zeile == "NO apriori altitude range") {
			getline(infile, Zeile);
			vector<double> dummy = string_to_vector<double>(Zeile);
			this->NO_apriori_bottom = dummy[0];
			this->NO_apriori_top = dummy[1];
			continue;
		}
		if (Zeile == "NO apriori scale") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> NO_apriori_scale;
			continue;
		}
		if (Zeile == "NO apriori smoothness") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> NO_apriori_smoothness;
			continue;
		}
		if (Zeile == "Retrieval algorithm") {
			getline(infile, Zeile);
			ss << Zeile;
			ss >> retrieval_algo;
			continue;
		}
		if (Zeile == "MLT") {
			int mlt = 0;
			getline(infile, Zeile);
			ss << Zeile;
			ss >> mlt;
			if (mlt == 0) MLT = false;
			else MLT = true;
			continue;
		}
	} //ende while !eof

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
	for (uint i = 0; i < this->m_AbsorptionsWL_der_Atmosphaerengase.size(); i++) {
		cout << "Absorptionswellenlänge " << i + 1 << ": " << this->m_AbsorptionsWL_der_Atmosphaerengase[i] << "\n";
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
	cout << "Skip SAA: ";
	if (skip_SAA) cout << "yes";
	else cout << "no";
	cout << endl;
	cout << "Max SZA: " << this->m_Maximaler_SZA << "\n";
	cout << "SAA cut-off: " << SAA_cutoff << "\n";
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
	cout << "S_alt: ";
	for (int i = 0; i < this->m_Anzahl_der_Emitter; i++) {
		cout << m_Retrieval_Kovarianzen[i + 2 * m_Anzahl_der_Emitter] << "\t";
	}
	cout << "\n";
	cout << "error thresholds: ";
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
	cout << "atmosphere temperature: " << this->atmo_Temp << "\n";
	cout << "number of NO transitions: " << this->no_NO_transitions << "\n";
	cout << "NO transitions:\n";
	for (unsigned i = 0; i < no_NO_transitions; i++) {
		cout << "v_u = " << NO_v_u.at(i) << ", v_l = " << NO_v_l.at(i)
			 << ", v_l_abs = " << NO_v_l_abs.at(i) << endl;
	}
	cout << "NO apriori: ";
	switch (NO_apriori) {
	case 1:
		cout << "NOEM";
		break;
	case 2:
		cout << "Regression";
		break;
	case 0:
	default:
		cout << "null";
	}
	cout << endl;
	cout << "NO apriori bottom: " << this->NO_apriori_bottom << "\n";
	cout << "NO apriori top: " << this->NO_apriori_top << "\n";
	cout << "NO apriori scale: " << this->NO_apriori_scale << "\n";
	cout << "NO apriori smoothness: " << this->NO_apriori_smoothness << "\n";
	cout << "Retrieval algorithm: ";
	switch (retrieval_algo) {
	case 0:
		cout << "old";
		break;
	case 1:
	default:
		cout << "new";
	}
	cout << "\n";
	cout << "MLT: ";
	if (MLT) cout << "true";
	else cout << "false";
	cout << endl;
	return 0;
}//ende Konfiguration_anzeigen

// HILFSFUNTION /////////////////////////////////////////////////////
template <class T> vector<T> string_to_vector(string zeile)
{
	vector<T> zahlen;
	T zahl;
	std::istringstream iss(zeile);

	while (iss >> zahl)
		zahlen.push_back(zahl);

	return zahlen;
}
// HILFSFUNTION ende /////////////////////////////////////////////////////


