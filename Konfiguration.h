#ifndef KONFIGURATION_HH_
#define KONFIGURATION_HH_

#include <string>
#include <vector>

using namespace std;


class Konfiguration
{
	//alles public machen
public:
	// Methoden
	int Konfiguration_einlesen();  // die ist in Scia.conf drin
	int Konfiguration_anzeigen();  // gucken obs geklappt hat
	// Member_Variablen
	// Directory Structures //////////////////////////////////////////
	int m_Anzahl_der_Emitter;
	string m_Pfad_Solar_Spektrum;
	string m_Pfad_Solar_Fallback_Spektrum;
	string m_Pfad_Linienparameter_Metalle;
	string m_Pfad_Dichten_der_Atmosphaerengase;
	string m_Pfad_Wirkungsquerschnitte_der_Atmosphaerengase;
	vector<double> m_AbsorbtionsWL_der_Atmosphaerengase;
	// Input_Data ////////////////////////////////////////////////////
	string m_Pfad_Datei_mit_Dateinamen_fuer_Messungen_eines_Orbits;
	string m_Pfad_Korrekturfaktoren;   // wird  nicht genutzt
	// Altitude Grid /////////////////////////////////////////////////
	double m_MinAlt;  // höhen für retrieval 71
	double m_MaxAlt;  // 93
	int m_Anzahl_zusaetzliche_Hoehengitterpunkte;
	vector<double> m_Grid_ext_low;
	vector<double> m_Grid_ext_high;
	double m_TOA; //top of Atmosphere, ab hier keine Absorption mehr
	// Switches für Selection Rules //////////////////////////////////
	int m_Nadir_only;  // switch 0 ja 1 nein
	int m_Nachtmessung;// switch 0 ja 1 nein
	int m_Geolocation; // switch 0 ja 1 nein
	int m_Large_SZA;   // switch 0 ja 1 nein //SZA Sonnenwinkel zum Zenit
	int m_NLC;         // switch 0 ja 1 nein
	double m_Maximaler_SZA;         // falls Auswahlkriterium Large_SZA aktiv ist dies die Grenze
	vector<double> m_Geolocation_Grenzen;//LonMin,LonMax,LatMin,LatMax
	// Baseline Fitparameters ////////////////////////////////////////
	int m_Anzahl_Baseline_Intervalle;
	vector<double> m_Baselinefenster_WL_low;
	vector<double> m_Baselinefenster_WL_high;
	// Columndensity Fitparameters ///////////////////////////////////
	int m_Anzahl_Retrieval_Intervalle;
	vector<double> m_Retrievalfenster_WL_low;
	vector<double> m_Retrievalfenster_WL_high;
	vector<int> m_Assignment_of_WL_Windows; /// zuordnung der linien zu den spezies
	// Regularisierungswichtungsfaktoren /////////////////////////////
	vector<double> m_Retrieval_Kovarianzen;
	// Sonstiges /////////////////////////////////////////////////////
	vector<double> m_Fehlergrenzen;   //sollte man auch dynamisch machen nach anzahl der spezies
	double m_FWHM;
	string m_Betriebssystem;
	int m_Do_Corrections_of_Radiances; // switch 0 ja 1 nein
	int m_Max_Zahl_Levenberg_Schritte;
	double m_Levenberg_Schrittweite;
	int m_Max_Zahl_Iterationen;
	double m_Convergence_Treshold;
};


#endif /* KONFIGURATION_HH_ */
