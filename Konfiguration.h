#ifndef KONFIGURATION_HH_
#define KONFIGURATION_HH_

#include <string>
#include <vector>

class Konfiguration
{
	//alles public machen
public:
	// default constructor
	Konfiguration();
	// Methoden
	void Konfiguration_einlesen(std::string file);  // die ist in Scia.conf drin
	void Konfiguration_anzeigen();  // gucken obs geklappt hat
	// Member_Variablen
	// Directory Structures //////////////////////////////////////////
	int m_Anzahl_der_Emitter;
	std::string m_Pfad_Solar_Spektrum;
	std::string m_Pfad_Solar_Fallback_Spektrum;
	std::string m_Pfad_Solar_Correction_Factors;
	std::string m_Pfad_Linienparameter_Metalle;
	std::string m_Pfad_Dichten_der_Atmosphaerengase;
	std::string m_Pfad_Wirkungsquerschnitte_der_Atmosphaerengase;
	std::string m_Pfad_NO_parameters;
	std::string m_Pfad_Ap_index;
	std::string m_Pfad_Kp_index;
	std::string m_Pfad_f107_index;
	std::string m_Pfad_f107a_index;
	std::string m_Pfad_f107_adj_index; // F10.7 adjusted flux (1 AU) for NOEM
	std::vector<double> m_AbsorptionsWL_der_Atmosphaerengase;
	// Input_Data ////////////////////////////////////////////////////
	std::string m_Pfad_Datei_mit_Dateinamen_fuer_Messungen_eines_Orbits;
	std::string m_Pfad_Korrekturfaktoren;   // wird  nicht genutzt
	// Altitude Grid /////////////////////////////////////////////////
	double m_MinAlt;  // höhen für retrieval 71
	double m_MaxAlt;  // 93
	double m_MinLat;
	double m_MaxLat;
	int m_NLat;
	int m_Anzahl_zusaetzliche_Hoehengitterpunkte;
	std::vector<double> m_Grid_ext_low;
	std::vector<double> m_Grid_ext_high;
	double m_TOA; //top of Atmosphere, ab hier keine Absorption mehr
	double m_BOA; //bottom of Atmosphere, darunte nur noch Absorption
	double m_min_TP, m_max_TP; // tangent points altitude range to use
	// Switches für Selection Rules //////////////////////////////////
	unsigned short m_Nadir_only;  // switch 0 ja 1 nein
	unsigned short m_Nachtmessung;// switch 0 ja 1 nein
	unsigned short m_Geolocation; // switch 0 ja 1 nein
	unsigned short m_Large_SZA;   // switch 0 ja 1 nein //SZA Sonnenwinkel zum Zenit
	unsigned short m_NLC;         // switch 0 ja 1 nein
	bool skip_SAA;
	double m_Maximaler_SZA;         // falls Auswahlkriterium Large_SZA aktiv ist dies die Grenze
	double SAA_cutoff;
	std::vector<double> m_Geolocation_Grenzen;//LonMin,LonMax,LatMin,LatMax
	// Baseline Fitparameters ////////////////////////////////////////
	int m_Anzahl_Baseline_Intervalle;
	std::vector<double> m_Baselinefenster_WL_low;
	std::vector<double> m_Baselinefenster_WL_high;
	// Columndensity Fitparameters ///////////////////////////////////
	int m_Anzahl_Retrieval_Intervalle;
	std::vector<double> m_Retrievalfenster_WL_low;
	std::vector<double> m_Retrievalfenster_WL_high;
	std::vector<int> m_Assignment_of_WL_Windows; /// zuordnung der linien zu den spezies
	// Regularisierungswichtungsfaktoren /////////////////////////////
	std::vector<double> m_Retrieval_Kovarianzen;
	// Sonstiges /////////////////////////////////////////////////////
	std::vector<double> m_Fehlergrenzen;   //sollte man auch dynamisch machen nach anzahl der spezies
	double m_FWHM;
	std::string m_Betriebssystem;
	unsigned short m_Do_Corrections_of_Radiances; // switch 0 ja 1 nein
	int m_Max_Zahl_Levenberg_Schritte;
	double m_Levenberg_Schrittweite;
	int m_Max_Zahl_Iterationen;
	double m_Convergence_Treshold;
	double atmo_Temp;
	unsigned no_NO_transitions;
	std::vector<int> NO_v_u, NO_v_l, NO_v_l_abs;
	bool NO_pol_correction; // polarisation correction of the NO lines
	unsigned short NO_rayleigh_fit_method;
	std::pair<double, double> NO_rayleigh_fit_window;
	unsigned short NO_apriori;
	double NO_apriori_bottom, NO_apriori_top;
	/* NO apriori scaling factor and transition width in km
	 * for a smooth transition from 0 to 1 (or to NO_apriori_scale). */
	double NO_apriori_scale, NO_apriori_smoothness;
	unsigned short retrieval_algo;
	bool MLT;
};


#endif /* KONFIGURATION_HH_ */
