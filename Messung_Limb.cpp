/*****************************************
 * Messung_Limb.cpp
 *
 *  Created on: 20.04.2010
 *      Author: martin
 *****************************************/


//Macro SAVEDELETE
#define SAVEDELETE(array)     \
if(array!=0)                                 \
{                                                \
    delete[] array;                        \
    array=0;                                 \
}                                                \
 

#include "Messung_Limb.h"
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <algorithm>

#include <fstream>  //für Ausgabe
#include <iostream>//für Ausgabe
#include <cstdlib>  //für Ausgabe
#include <cstdio>   //Filekram
#include "Ausdrucke.h"
#include "Fit_Polynom.h"
#include "Glaetten.h"

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
} //Lapackroutine zur Lösung eines linearen GLeichungssystems


using namespace std;

Messung_Limb::Messung_Limb()
{
	//initialisierung
	m_Zeilendichte = 0;
	m_Fehler_Zeilendichten = 0;
	m_Dateiname_L1C = "dummy";
	m_Number_of_Wavelength = 0;
	m_Jahr = 0;
	m_Monat = 0;
	m_Tag = 0;
	m_Stunde = 0;
	m_Minute = 0;
	m_Deklinationswinkel = 0;
	m_Sonnen_Longitude = 0;
	m_Lattidude_TP = 0;
	m_Longitude_TP = 0;
	m_Lattidude_Sat = 0;
	m_Longitude_Sat = 0;
	m_Hoehe_TP = 0;
	m_Hoehe_Sat = 0;
	m_Erdradius = 0;
	//statische Felder werden erstmal nicht 0 gesetzt
}
//========================================
//
//copyconstructor
//
//========================================
Messung_Limb::Messung_Limb(const Messung_Limb &rhs)
{
	*this = rhs;
}//copyconstructor ende

//========================================

//========================================
//
// Assignmentoperator Overload
//
//========================================
Messung_Limb &Messung_Limb::operator =(const Messung_Limb &rhs)
{
	//TODO das nochmal anpassen
	// Prevent self assignment. We say two Strings
	// are equal if their memory addresses are equal.
	if (this == &rhs)
		return *this;
	// Ergebnisse
	m_Zeilendichte = rhs.m_Zeilendichte;
	m_Fehler_Zeilendichten = rhs.m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	m_Deklinationswinkel = rhs.m_Deklinationswinkel;
	m_Sonnen_Longitude = rhs.m_Sonnen_Longitude;
	// Dateiname
	m_Dateiname_L1C = rhs.m_Dateiname_L1C;
	// Datum
	m_Jahr = rhs.m_Jahr;
	m_Monat = rhs.m_Monat;
	m_Tag = rhs.m_Tag;
	m_Stunde = rhs.m_Stunde;
	m_Minute = rhs.m_Minute;
	// Geolocation
	m_Lattidude_Sat = rhs.m_Lattidude_Sat;
	m_Longitude_Sat = rhs.m_Longitude_Sat;
	m_Hoehe_Sat = rhs.m_Hoehe_Sat;
	m_Lattidude_TP = rhs.m_Lattidude_TP;
	m_Longitude_TP = rhs.m_Longitude_TP;
	m_Hoehe_TP = rhs.m_Hoehe_TP;
	m_Erdradius = rhs.m_Erdradius;
	m_TP_SZA = rhs.m_TP_SZA;
	m_Number_of_Wavelength = rhs.m_Number_of_Wavelength;
	// copy vectors
	m_Wellenlaengen = rhs.m_Wellenlaengen;
	m_Intensitaeten = rhs.m_Intensitaeten;
	m_Sonne = rhs.m_Sonne;
	m_Intensitaeten_durch_piF = rhs.m_Intensitaeten_durch_piF;
	m_Intensitaeten_durch_piF_Gamma = rhs.m_Intensitaeten_durch_piF_Gamma;
	m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand
		= rhs.m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand;

	// Return a reference to *this object.
	return *this;
}// Assignmentoperator Overload ende
//========================================
//========================================
//Methoden
//========================================
//========================================
int Messung_Limb::Zeilendichte_Bestimmen(Speziesfenster &Spezfenst, int Index,
		string Arbeitsverzeichnis, string mache_Fit_Plots)
{
	//kurz:
	//Diese Fitroutine ermittelt die Fläche von I/(piF*Gamma) bei der
	//Übergangswellenlänge Da nur wenige Punkte zu einer Linie beitragen,
	//entspricht der Linie der Auflösungsfunktion
	//Diese ist ein Hyperboloide 4ter Ordnung

	//lang:
	// Die Umgebung jeder Linie wird in 3 Bereiche unterteilt. 2 sogenannte
	// Basisfenster links und rechts von der Linie und das Peakfenster selbst,
	// welches in der Umsetzung auch Bereiche aus den beiden Basisfenstern
	// beinhalten darf.  Bis Hierhin haben wir I/(piFGamma) berechnet. Nun
	// wollen wir die "normierte" Intensität für einen Peak berechnen Da der
	// Peak gerade im Kanal zumeist auf einem relativ großen Untergrundsignal
	// sitzt, wird dieses zunächst abgezogen.  Dafür wird im Bereich um die
	// Linie herum, wo keine signifikanten anderen Linien liegen(Basisfenster)
	// ein linearer Fit der Grundlinie durchgeführt, welcher vom der
	// "normierten" Intensität abgezogen wird.  Danach wird im Bereich des
	// Peakfensters die Fläche des Peaks bestimmt.  Da die Linie zumeist nur
	// aus 3 Punkten besteht,  hat das Profil im wesentlichen die Form der
	// Auflösungsfunktion, welche eine hyperboloide 4ten Grades ist. Diese hat
	// im wesentlichen 3 Parameter: Wellenlänge des Peaks, Breite und Fläche.
	// Die robusteste Variante ist es, die Wellenlänge und die Halbwertsbreite
	// der Linie vorzugeben und nur die Fläche anzufitten.  Diese Variante kann
	// durch einen linearen Parameterfit erreicht werden und ist zum einen
	// schnell und was hier viel wichtiger ist robust. Andere nichtlineare
	// Verfahren sind entweder nicht robust(Gauss Newton, insbesondere bei dem
	// schlechten Signal/Noise Verhältnis) oder deutlich langsamer (simulated
	// annealing bzw. ist das aufwändiger dies noch schnell zu implementieren...
	// das wäre Schritt 2)
	// Die ermittelte Fläche ist dann unsere gesuchte Säulenzeilendichte.

	// I/(piFGamma)=integral(AMF n ds) mit AMF = s exp(-tau) ...aber zu der
	// Formel später nochmal zurück Das spätere Retrieval ermittelt dann die
	// Dichte n aus der rechten Seite

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
	int Index_Peakfenster_min = Get_Index(Spezfenst.m_Peakfenster_WLmin[Index]);
	int Index_Peakfenster_max = Get_Index(Spezfenst.m_Peakfenster_WLmax[Index]);
	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int N_Basis = Bas_l + Bas_r;
	int N_Peak = Index_Peakfenster_max - Index_Peakfenster_min + 1;
	// Speicher anfordern
	vector<double> Basisfenster_WL(N_Basis);
	vector<double> Basisfenster_Intensitaet(N_Basis);
	vector<double> Peakfenster_WL(N_Peak);
	vector<double> Peakfenster_Intensitaet(N_Peak);
	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i] = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] = this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_rechts_min + i];
	}
	//Peakfenster WL und I auffüllen
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_WL[i] = m_Wellenlaengen[Index_Peakfenster_min + i];
		Peakfenster_Intensitaet[i] =
			m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Peakfenster_min + i];
	}
	// linearen Fit des Basisfensters durchführen
	// Proto: Fit_Linear(double* x,double* y, double& a0, double& a1,int
	// Anfangsindex, int Endindex)
	double a0, a1;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, 0,
			N_Basis - 1);
	// lineare Funktion von Intensitäten des Peakfenster abziehen
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_Intensitaet[i] -= a0 + a1 * Peakfenster_WL[i];
	}
	// Hyperboloiden an Peakfenster anfitten
	// Proto: Fit_Peak_hyperbolic(double* x,double* y,double x0, double FWHM,
	// double& A, int Anfangsindex, int Endindex)
	double Flaeche;
	Fit_Peak_hyperbolic(Peakfenster_WL, Peakfenster_Intensitaet,
						Spezfenst.m_Wellenlaengen[Index],
						Spezfenst.m_FWHM, Flaeche, 0, N_Peak - 1);
	// Hier Wellenlängen in nm verwendet..das hebt sich mit dem Gitterabstand
	// raus
	//Fehler des Fits bestimmen... da Peakfenster_Intensitaet nicht mehr
	//gebraucht wird die Basislinie für die Fehlerberechnung wieder aufaddiert
	m_Zeilendichte = Flaeche;
	//cout<<m_Zeilendichte<<"\n";
	//cout<<Spezfenst.m_Liniendaten[Index].m_Gamma<<"\n";
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_Intensitaet[i] += a0 + a1 * Peakfenster_WL[i];
	}
	// Funktion double Messung_Limb::Evaluate_Error_primitive(double* x,
	// double* y, double a0,double a1, double A, double FWHM, double x0,
	// int Anfangsindex, int Endindex)
	m_Fehler_Zeilendichten = Evaluate_Error_primitive(Peakfenster_WL,
			Peakfenster_Intensitaet, a0, a1, m_Zeilendichte, Spezfenst.m_FWHM,
			Spezfenst.m_Wellenlaengen[Index], 0, N_Peak - 1);

	////////////////////////////////////////////////////////////////////////////
	// Hier kann man zur Testzwecken noch einen Plot machen  ///////////////////
	if (mache_Fit_Plots == "ja") {
		//TODO das als Funktion implementieren
		vector<double> Funktion(N_Peak);
		const double pi = M_PI;
		double FWHM = Spezfenst.m_FWHM;
		double cnorm = 4.0 * pi * M_SQRT2 / (FWHM * FWHM * FWHM);

		for (int i = 0; i < N_Peak; i++) {
			double Basis = a0 + a1 * Peakfenster_WL[i];
			double Peak = m_Zeilendichte /
				(cnorm * (pow(0.5 * FWHM, 4)
				 + pow(Peakfenster_WL[i] - Spezfenst.m_Wellenlaengen[Index], 4)));
			//cout<<Peak<<"\n";
			//cout<<m_Zeilendichte<<"\n";
			Funktion[i] = Peak + Basis;
		}

		string s1, s_OrbNum, s2;
		int h = this->m_Hoehe_TP;
		char buf[256];
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...
		// falls / im Namen ist das schlecht
		string Datnam = m_Dateiname_L1C.substr(m_Dateiname_L1C.size() - 39, 39);

		//TODO Pfad anpassen
		sprintf(buf, "mkdir %s/Plots 2>/dev/null", Arbeitsverzeichnis.c_str());

		string Befehl = buf;
		system(Befehl.c_str());
		sprintf(buf, "%s_%s_%i_%ikm.ps", Datnam.c_str(),
				Spezfenst.m_Spezies_Name.c_str(), Index, h);
		string new_datnam = buf;
		sprintf(buf, "%s/Plots/%s", Arbeitsverzeichnis.c_str(),
				new_datnam.c_str());
		s1 = buf;
		//s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern, um
		//später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal wie die Datei heißt, die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Limbdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		sprintf(buf, "Orbit %5s Limb TP: Lat: %G Grad  Lon: %G Grad Hoehe: %G km",
				s_OrbNum.c_str(), m_Lattidude_TP, m_Longitude_TP, m_Hoehe_TP);
		s2 = buf;
		//cout<<s1<<"\n";
		//int Plot_2xy(string Dateiname,string title, string xlabel,
		//string ylabel,double* x1,double*y1, double* x2,double* y2,
		//int Startindex,int Endindex);
		//Plot_2xy(s1.c_str(),s1.substr(s1.size()-50,50).c_str(),
		//"$\\lambda$ in nm","$\\frac{I}{\\piF\\gamma}$",Peakfenster_WL,
		//Peakfenster_Intensitaet,Peakfenster_WL,Funktion,0,N_Peak-1);
		//-> Fit geht
		Plot_2xy(Arbeitsverzeichnis.c_str(), s1.c_str(), s2.c_str(),
				 "Wellenlaenge in nm",
				 "Schraege Saeule bei Peakposition in cm^{-2}/nm",
				 Peakfenster_WL, Peakfenster_Intensitaet, Peakfenster_WL,
				 Funktion, 0, N_Peak - 1,
				 m_Zeilendichte, m_Fehler_Zeilendichten);
	}
	// Ende Plot ///////////////////////////////////////////////////////////////

	return 0;
}//int Zeilendichte_Bestimmen() ende
//========================================
//========================================

////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Saeulendichte_Bestimmen_MgI285nm
////////////////////////////////////////////////////////////////////////////////
int Messung_Limb::Saeulendichte_Bestimmen_MgI285nm(Speziesfenster &Spezfenst,
		int Index, string Arbeitsverzeichnis, string mache_Fit_Plots,
		double *mean_10_20)
{
	//kurz:
	//alternative Fitroutine für die Bestimmung  der MgI 285.21275 nm Linie
	//lang:
	// Das Signal zu Rausch Verhältnis für die Limbmessungen scheint ziemlich
	// mies zu sein, sodass die Befürchtung besteht, dass man da Rauschen als
	// Messwerte interpretiert
	// Ein Erster Rettungsversuch ist es, statt den Qoutienten I/piF direkt aus
	// den Messwerten zu bilden, die breiten Peaks im Limb und Sonnenspektrum
	// auszunutzen
	//
	// I/(Fgamma * oder / Integrations wegstück in nm) ist schon vorhanden
	// in dem entsprechendem Fenster wird die Basislinie diese Quotienten der
	// Messwerte gebildet Limb und Sonnenspektrum werden einzeln als Polynome
	// angefittet
	// Da das Minimum eines nicht geraden Polynoms nicht bei der vorgegebenen
	// Wellenlänge liegt wird das Limbspektrum so geshifted, dass die Minima
	// übereinander liegen
	// Der Quotient beider Polynome wird gebildet
	// Der konstante Faktor aus gamma und infinitesimalen Wegstück(ca. 0.11
	// also WL(2)-WL(1)) wird anmultipliziert, sodass Quotient aus Messwerten
	// und Quotient aus Polynomen die gleiche Größenordnung haben
	// Die Basislinie wird nun von beiden abgezogen
	// Wenn überhaupt, sollte man erst jetzt den Quotienten so shiften, dass
	// das Minimum bei der Wellenlänge des Übergangs liegt.  Abschließend
	// werden für jede Messung mehrere Plots erstellt

	double *Basisfenster_WL;
	double *Basisfenster_Intensitaet;
	double *Vollfenster_WL;
	double *Vollfenster_Limb;
	double *Vollfenster_Sonne;

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
	int Index_Uebergangs_Wellenlaenge = Get_Index(Spezfenst.m_Wellenlaengen[Index]) - Index_Basisfenster_links_min;

	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int N_Basis = Bas_l + Bas_r;
	//int N_Peak=Index_Peakfenster_max-Index_Peakfenster_min+1;
	int N_Vollfenster = Index_Basisfenster_rechts_max - Index_Basisfenster_links_min + 1;
	// Speicher anfordern
	Basisfenster_WL = new double[N_Basis];
	Basisfenster_Intensitaet = new double[N_Basis];
	Vollfenster_WL = new double[N_Vollfenster];
	Vollfenster_Limb = new double[N_Vollfenster];
	Vollfenster_Sonne = new double[N_Vollfenster];

	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i] = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] = this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_rechts_min + i];
	}

	//Vollfenster Limb und Sonne auffüllen
	for (int i = 0; i < N_Vollfenster; i++) {
		Vollfenster_WL[i] = m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Vollfenster_Limb[i] = m_Intensitaeten[Index_Basisfenster_links_min + i];
		Vollfenster_Sonne[i] = m_Sonne[Index_Basisfenster_links_min + i];
	}

	double *Vollfenster_Limb_mittlere_atmo; //zwischen 40 und 60km
	Vollfenster_Limb_mittlere_atmo = new double[N_Vollfenster];
	for (int i = 0; i < N_Vollfenster; i++) {
		Vollfenster_Limb_mittlere_atmo[i] = mean_10_20[Index_Basisfenster_links_min + i];
	}
	//Minima vergleichen um Linie herum (5 nachbarpunkte)
	// und shift durchführen
	// Suche mit Glattem Fenster, Verschiebung im unglatten
	int Suchbereich = 5;
	int min_Index_Limb = Index_Uebergangs_Wellenlaenge - Suchbereich;
	int min_Index_Sun = Index_Uebergangs_Wellenlaenge - Suchbereich;
	for (int i = Index_Uebergangs_Wellenlaenge - Suchbereich;
			i <= Index_Uebergangs_Wellenlaenge + Suchbereich; i++) {
		if (Vollfenster_Limb_mittlere_atmo[i] < Vollfenster_Limb_mittlere_atmo[min_Index_Limb]) {
			min_Index_Limb = i;
		}
		if (Vollfenster_Sonne[i] < Vollfenster_Sonne[min_Index_Sun]) {
			min_Index_Sun = i;
		}
	}
	int Verschiebung = min_Index_Sun - min_Index_Limb;
	if ((Verschiebung >= 1) || (Verschiebung <= -1)) {
		cout << "Verschiebung: " << Verschiebung << "\n";
		cout << "Minima auseinander\n";
	}
	if ((Verschiebung >= 3) || (Verschiebung <= -3)) {
		cout << "Verschiebung: " << Verschiebung << "\n";
		cout << "Minima weit auseinander\n";
	}
	double *vor_Verschiebung;
	vor_Verschiebung = new double[N_Vollfenster];
	for (int i = 0; i < N_Vollfenster; i++) {
		vor_Verschiebung[i] = Vollfenster_Limb[i];
	}
	for (int i = 0; i < N_Vollfenster; i++) {
		if (!(((i + Verschiebung) < 0)
			|| (i + Verschiebung >= N_Vollfenster))) {
			Vollfenster_Limb[i + Verschiebung] = vor_Verschiebung[i];
		}
	}
	delete[] vor_Verschiebung;
	delete[] Vollfenster_Limb_mittlere_atmo;
	// Spektrum glätten
	int smooth_Nachbarn = 0; //2;
	int smooth_Iterationen = 0; //4;
	smooth_data(N_Vollfenster, Vollfenster_Limb, smooth_Nachbarn, smooth_Iterationen);
	smooth_data(N_Vollfenster, Vollfenster_Sonne, smooth_Nachbarn, smooth_Iterationen);

	// linearen Fit des Basisfensters durchführen
	// Proto:
	// Fit_Linear(double* x,double* y, double& a0, double& a1,
	//   int Anfangsindex, int Endindex)
	double a0, a1;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, 0,
			N_Basis - 1);

	int Polynomgrad = 4;
	//Get_Index(Spezfenst.m_Wellenlaengen[Index]) darf hier nicht benutzt werden
	int Fitpunktezahl_links = 4;
	int Polyfit_Startindex = Index_Uebergangs_Wellenlaenge - Fitpunktezahl_links;
	int Polyfit_Endindex = Index_Uebergangs_Wellenlaenge + Fitpunktezahl_links;
	double *Limbfit_Parameter;
	double *Sonnefit_Parameter;
	Limbfit_Parameter = new double[Polynomgrad + 1];
	Sonnefit_Parameter = new double[Polynomgrad + 1];
	// Limb-Spektrum mit einem Polynom anfitten
	// später evlt mehrmals mit shifts in beide Richtungen
	// evtl runs-test machen
	Fit_Polynom(Vollfenster_WL, Vollfenster_Limb, Polyfit_Startindex,
				Polyfit_Endindex, Spezfenst.m_Wellenlaengen[Index], Polynomgrad,
				Limbfit_Parameter);
	// Fit sieht ok aus
	// Sonnen-Spektrum mit gleichem Polynom anfitten
	// später evlt mehrmals mit shifts in beide Richtungen
	Fit_Polynom(Vollfenster_WL, Vollfenster_Sonne, Polyfit_Startindex,
				Polyfit_Endindex, Spezfenst.m_Wellenlaengen[Index], Polynomgrad,
				Sonnefit_Parameter);
	// Beide Polynome im gesamten Fenster diskretisieren...
	// dabei 5 mal so viele Punkte nutzen, wie die ursprünglichen Messwerte
	double *Vollfenster_fein_WL;
	double *Vollfenster_fein_Limb;
	double *Vollfenster_fein_Sonne;
	int N_Vollfenster_fein = 1 + (N_Vollfenster - 1) * 5;
	Vollfenster_fein_WL = new double[N_Vollfenster_fein];
	Vollfenster_fein_Limb = new double[N_Vollfenster_fein];
	Vollfenster_fein_Sonne = new double[N_Vollfenster_fein];
	for (int i = 0; i < (N_Vollfenster - 1); i++) {
		Vollfenster_fein_WL[5 * i] = Vollfenster_WL[i];
		Vollfenster_fein_WL[5 * i + 1] = Vollfenster_WL[i] + 0.2 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
		Vollfenster_fein_WL[5 * i + 2] = Vollfenster_WL[i] + 0.4 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
		Vollfenster_fein_WL[5 * i + 3] = Vollfenster_WL[i] + 0.6 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
		Vollfenster_fein_WL[5 * i + 4] = Vollfenster_WL[i] + 0.8 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
	}
	//letzter Punkt, nicht verfünffachen
	Vollfenster_fein_WL[N_Vollfenster_fein - 1] = Vollfenster_WL[N_Vollfenster - 1];
	//Limb und Sonnenspektrum für diskrete Punkte ausrechnen
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		Vollfenster_fein_Limb[i] = 0;
		Vollfenster_fein_Sonne[i] = 0;
	}
	for (int i = 5 * (Polyfit_Startindex - 3); i <= 5 * (Polyfit_Endindex + 3); i++) {
		double h = Vollfenster_fein_WL[i] - Spezfenst.m_Wellenlaengen[Index]; //x-x0
		double Faktor = 1.0;
		Vollfenster_fein_Limb[i] = 0;
		Vollfenster_fein_Sonne[i] = 0;
		for (int j = 0; j <= Polynomgrad; j++) {
			Vollfenster_fein_Limb[i] += Limbfit_Parameter[j] * Faktor;
			Vollfenster_fein_Sonne[i] += Sonnefit_Parameter[j] * Faktor;
			Faktor *= h;
		} //ende for j
	}// ende for i
	//Minimimum beider Polynome bestimmen und Limbspektrum so shiften, dass das
	//Minimum auch auf dem Minimum des Sonnenspektrums liegt
	double Limb_WL_min, Limb_y_min;
	int Limb_Indexmin;
	double Sonne_WL_min, Sonne_y_min;
	int Sonne_Indexmin;
	x_zu_Minimum_von_y_in_Intervall(Vollfenster_fein_WL, Vollfenster_fein_Limb,
									Polyfit_Startindex * 5, Polyfit_Endindex * 5,
									Limb_WL_min, Limb_y_min, Limb_Indexmin);
	x_zu_Minimum_von_y_in_Intervall(Vollfenster_fein_WL, Vollfenster_fein_Sonne,
									Polyfit_Startindex * 5, Polyfit_Endindex * 5,
									Sonne_WL_min, Sonne_y_min, Sonne_Indexmin);
	int shift = Sonne_Indexmin - Limb_Indexmin; //z.b. Sonne 43 Limb 50 shift=-7
	// shiften des Limbspektrums
	// Die Randpunkte sind jetzt falsch, aber die Linie sollte nicht am Rand
	// liegen
	double *vor_shift;
//    cout<<"Sonne_WL_min: "<<Sonne_WL_min<<"\n";
//    cout<<"Limb_WL_min: "<<Limb_WL_min<<"\n";
//    cout<<"Limb_Indexmin: "<<Limb_Indexmin<<"\n";
//    cout<<"Sonne_Indexmin: "<<Sonne_Indexmin<<"\n";
//    cout<<"shift: "<<shift<<"\n";

	vor_shift = new double[N_Vollfenster_fein];
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		vor_shift[i] = Vollfenster_fein_Limb[i];
	}
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		if (!(((i + shift) < 0) || (i + shift >= N_Vollfenster_fein))) {
			Vollfenster_fein_Limb[i + shift] = vor_shift[i];
		}
	}
	delete[] vor_shift;
	//Die beiden gefitteten Spektren dividieren Limb/Sonne
	double *Fit_Quotient;
	double *Messwerte_Quotient;
	Messwerte_Quotient = new double[N_Vollfenster];
	Fit_Quotient = new double[N_Vollfenster_fein];
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] = Vollfenster_Limb[i] / Vollfenster_Sonne[i];
	}
	int Ind = Get_Index(Spezfenst.m_Wellenlaengen[Index]);
	double Delta_WL = (m_Wellenlaengen[Ind + 1] - m_Wellenlaengen[Ind]);
	double Umrechnung = 1 / (Delta_WL * Spezfenst.m_Liniendaten[Index].m_Gamma);
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] *= Umrechnung;
	}
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		Fit_Quotient[i] = Umrechnung * Vollfenster_fein_Limb[i]
						  / Vollfenster_fein_Sonne[i];
		//cout<<Vollfenster_fein_WL[i]<<"\t"<<Fit_Quotient[i]<<"\n";

	}

	//Basislinie abziehen
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] -= a0 + a1 * Vollfenster_WL[i];
	}
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		Fit_Quotient[i]             -= a0 + a1 * Vollfenster_fein_WL[i];
	}
	// Falls an der Stelle des Peaks der Wert jetzt negativ ist, alle positiven
	// Werte abschneiden, sonst negative
	if (Fit_Quotient[Sonne_Indexmin] < 0) {
		for (int i = 0; i < N_Vollfenster_fein; i++) {
			if (Fit_Quotient[i] > 0) {
				Fit_Quotient[i] = 0;
			}
		}
	} else {
		for (int i = 0; i < N_Vollfenster_fein; i++) {
			if (Fit_Quotient[i] < 0) {
				Fit_Quotient[i] = 0;
			}
		}
	}
	//Alles Abschneiden, was nicht zum peak beiträgt
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		if ((Vollfenster_fein_WL[i] < 284.8) || (Vollfenster_fein_WL[i] > 285.6)) {
			Fit_Quotient[i] = 0;
		}
	}

	// Noch ebend eine Fitfunktion in den Quotienten legen
	double Flaeche;
	// Hier Wellenlängen in nm verwendet..
	// as hebt sich mit dem Gitterabstand raus
	Fit_Peak_hyperbolic(Vollfenster_fein_WL, Fit_Quotient, Sonne_WL_min,
						Spezfenst.m_FWHM, Flaeche, 0,
						N_Vollfenster_fein - 1);

	//Plotroutine aufrufen
	if (mache_Fit_Plots == "ja") {
		// Dateinamenschnickschnack
		string s1, s_OrbNum, s2;
		int h = this->m_Hoehe_TP;
		char buf[256];
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...
		//falls / im Namen ist das schlecht
		string Datnam = m_Dateiname_L1C.substr(m_Dateiname_L1C.size() - 39, 39);
		//TODO Pfad anpassen
		sprintf(buf, "mkdir %s/Plots 2>/dev/null", Arbeitsverzeichnis.c_str());
		string Befehl = buf;
		system(Befehl.c_str());
		sprintf(buf, "%s_%s_%i_%ikm.ps",
				Datnam.c_str(), Spezfenst.m_Spezies_Name.c_str(), Index, h);
		string new_datnam = buf;
		sprintf(buf, "%s/Plots/%s",
				Arbeitsverzeichnis.c_str(), new_datnam.c_str());
		s1 = buf;
		//s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern, um
		//später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal, wie die Datei heißt die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Limbdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		sprintf(buf, "Orbit %5s Limb TP: Lat: %G {/Symbol \\260}  Lon: %G {/Symbol \\260} Hoehe: %G km",
				s_OrbNum.c_str(), m_Lattidude_TP, m_Longitude_TP, m_Hoehe_TP);
		s2 = buf;
		//cout<<s1<<"\n";
		Plot_Slantcoloumns_polyfit_MgI(Arbeitsverzeichnis.c_str(), s1.c_str(),
									   s2.c_str(),
									   Vollfenster_WL,
									   0.35 * (N_Vollfenster - 1),
									   0.8 * (N_Vollfenster - 1),
									   Vollfenster_fein_WL,
									   0.35 * (N_Vollfenster_fein - 1),
									   0.8 * (N_Vollfenster_fein - 1),
									   Vollfenster_Limb, Vollfenster_fein_Limb,
									   Vollfenster_Sonne, Vollfenster_fein_Sonne,
									   Messwerte_Quotient, Fit_Quotient);
	}
	//Speicher freimachen
	delete[] Basisfenster_WL;
	delete[] Basisfenster_Intensitaet;
	delete[] Vollfenster_WL;
	delete[] Vollfenster_Limb;
	delete[] Vollfenster_Sonne;
	delete[] Vollfenster_fein_WL;
	delete[] Vollfenster_fein_Limb;
	delete[] Vollfenster_fein_Sonne;
	delete[] Limbfit_Parameter;
	delete[] Sonnefit_Parameter;

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Saeulendichte_Bestimmen_MgI285nm
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plots_der_Spektren_erzeugen
////////////////////////////////////////////////////////////////////////////////
int Messung_Limb::Plots_der_Spektren_erzeugen(Speziesfenster &Spezfenst,
		int Index, string Arbeitsverzeichnis, string mache_Fit_Plots,
		double *mean_10_20)
{
	// Plot der Spektren Sonne und Limb und Quotient mit und ohne Rauschen
	double *Basisfenster_WL;
	double *Basisfenster_Intensitaet;
	double *Vollfenster_WL;
	double *Vollfenster_Limb;
	double *Vollfenster_Sonne;

	double *Vollfenster_Limb_abs_error;

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
//    int Index_Uebergangs_Wellenlaenge=Get_Index(Spezfenst.m_Wellenlaengen[Index])-Index_Basisfenster_links_min;
	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int N_Basis = Bas_l + Bas_r;
	//int N_Peak=Index_Peakfenster_max-Index_Peakfenster_min+1;
	int N_Vollfenster = Index_Basisfenster_rechts_max - Index_Basisfenster_links_min + 1;
	// Speicher anfordern
	Basisfenster_WL = new double[N_Basis];
	Basisfenster_Intensitaet = new double[N_Basis];
	Vollfenster_WL = new double[N_Vollfenster];
	Vollfenster_Limb = new double[N_Vollfenster];
	Vollfenster_Sonne = new double[N_Vollfenster];

	Vollfenster_Limb_abs_error = new double[N_Vollfenster];

	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i] = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] =
			this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_rechts_min + i];
	}

	//Vollfenster Limb und Sonne auffüllen
	for (int i = 0; i < N_Vollfenster; i++) {
		Vollfenster_WL[i] = m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Vollfenster_Limb[i] = m_Intensitaeten[Index_Basisfenster_links_min + i];
		Vollfenster_Sonne[i] = m_Sonne[Index_Basisfenster_links_min + i];
		Vollfenster_Limb_abs_error[i] = Vollfenster_Limb[i]
			* m_Intensitaeten_relativer_Fehler[Index_Basisfenster_links_min + i];
	}

	// TODO braucht man das hier überhaupt noch für irgendwas...wenn nicht weg
	// damit..das verwirrt nur
	//double* Vollfenster_Limb_mittlere_atmo; //zwischen 40 und 60km
	//Vollfenster_Limb_mittlere_atmo=new double[N_Vollfenster];
	//for(int i=0;i<N_Vollfenster;i++)
	//{
	//  Vollfenster_Limb_mittlere_atmo[i]
	//    = mean_10_20[Index_Basisfenster_links_min+i];
	//}

	// Spektrum glätten bei 0 wird nichts geglättet
	// ...könnte man auch auskommentieren
	int smooth_Nachbarn = 0; //1;//2;
	int smooth_Iterationen = 0; //6;//4;
	smooth_data(N_Vollfenster, Vollfenster_Limb, smooth_Nachbarn, smooth_Iterationen);
	smooth_data(N_Vollfenster, Vollfenster_Sonne, smooth_Nachbarn, smooth_Iterationen);

	// linearen Fit des Basisfensters durchführen
	// TODO Basisfenster neu berechnen
	// Proto:
	// Fit_Linear(double* x,double* y, double& a0, double& a1,int Anfangsindex,
	// int Endindex)
	double a0, a1;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, 0,
			N_Basis - 1);
	//Die beiden gefitteten Spektren dividieren Limb/Sonne
	double *Messwerte_Quotient;
	double *Messwerte_Quotient_error;  // Es wird angenommen, der Fehler liegt nur in Limb vor
	Messwerte_Quotient = new double[N_Vollfenster];
	Messwerte_Quotient_error = new double[N_Vollfenster];
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] = Vollfenster_Limb[i] / Vollfenster_Sonne[i];
		// Die relativen Fehler von MW_Q und Vf_L sind gleich...Limb ist linear,
		// also auch Fehler linear skalieren
		Messwerte_Quotient_error[i] = Vollfenster_Limb_abs_error[i]
									  / Vollfenster_Sonne[i];

	}
	int Ind = Get_Index(Spezfenst.m_Wellenlaengen[Index]);
	double Delta_WL = (m_Wellenlaengen[Ind + 1] - m_Wellenlaengen[Ind]);
	double Umrechnung = 1 / (Delta_WL * Spezfenst.m_Liniendaten[Index].m_Gamma);
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] *= Umrechnung;
		Messwerte_Quotient_error[i] *= Umrechnung;
	}

	//Basislinie abziehen
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] -= a0 + a1 * Vollfenster_WL[i];
		// Im besten Fall verändert sich der absolute Fehler nicht...aber der
		// relative, wenn die Baseline im Vergleich zum peak hoch liegt
	}
	// Fehler aus Residuuen abschätzen
	/////////////////////(und nicht den gegebenen Fehler nehmen
	// Zunächst nochmal den Mittelwert bilden
	double *Messwerte_Quotient_stabw;
	Messwerte_Quotient_stabw = new double[N_Vollfenster];
	double Mean = 0;
	for (int i = 0; i < N_Vollfenster; i++)  {
		Mean += Messwerte_Quotient[i];
	}
	Mean /= N_Vollfenster;
	double standardabweichung = 0;
	for (int i = 0; i < N_Vollfenster; i++) {
		standardabweichung += (Messwerte_Quotient[i] - Mean)
							  * (Messwerte_Quotient[i] - Mean);
	}
	standardabweichung /= N_Vollfenster - 1;
	standardabweichung = sqrt(standardabweichung);
	for (int i = 0; i < N_Vollfenster; i++)    {
		Messwerte_Quotient_stabw[i] = standardabweichung;
	}
	// Ende Fehler aus Residuum abschätzen //////////////

	//Plotroutine aufrufen
	if (mache_Fit_Plots == "ja") {
		// Dateinamenschnickschnack
		string s1, s_OrbNum, s2;
		int h = this->m_Hoehe_TP;
		char buf[256];
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt
		// ...falls / im Namen ist das schlecht
		string Datnam = m_Dateiname_L1C.substr(m_Dateiname_L1C.size() - 39, 39);
		//TODO Pfad anpassen
		sprintf(buf, "mkdir %s/Plots 2>/dev/null", Arbeitsverzeichnis.c_str());
		string Befehl = buf;
		system(Befehl.c_str());
		sprintf(buf, "%s_%s_%i_%ikm.ps",
				Datnam.c_str(), Spezfenst.m_Spezies_Name.c_str(), Index, h);
		string new_datnam = buf;
		sprintf(buf, "%s/Plots/%s",
				Arbeitsverzeichnis.c_str(), new_datnam.c_str());
		s1 = buf;
		//s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern,
		//um später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal, wie die Datei heißt die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Limbdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		sprintf(buf, "Orbit %5s Limb TP: Lat: %G {/Symbol \\260}  Lon: %G {/Symbol \\260} Hoehe: %G km",
				s_OrbNum.c_str(), m_Lattidude_TP, m_Longitude_TP, m_Hoehe_TP);
		s2 = buf;
		//cout<<s1<<"\n";
		Plot_Spektren_und_Quotient(Arbeitsverzeichnis.c_str(),
								   s1.c_str(), s2.c_str(),
								   Vollfenster_WL, 0,
								   N_Vollfenster - 1,
								   Vollfenster_Limb, Vollfenster_Limb_abs_error,
								   Vollfenster_Sonne, Messwerte_Quotient,
								   Messwerte_Quotient_error);
		/*Plot_Quotient_mit_Fehler(Arbeitsverzeichnis.c_str(),s1.c_str(), s2.c_str(),
		                   Vollfenster_WL,0 ,N_Vollfenster-1,
		                   Messwerte_Quotient,Messwerte_Quotient_error,
		                   Messwerte_Quotient_stabw);*/
	}
	///////////////////////////////////////////////
	// Fit der Säulendichte  /////////////
	double Flaeche;
	// Hier Wellenlängen in nm verwendet..
	// das hebt sich mit dem Gitterabstand raus
	Fit_Peak_hyperbolic(Vollfenster_WL, Messwerte_Quotient,
						Spezfenst.m_Wellenlaengen[Index],
						Spezfenst.m_FWHM, Flaeche, 0,
						N_Vollfenster - 1);
	//Fehler des Fits bestimmen... da Peakfenster_Intensitaet nicht mehr
	//gebraucht wird die Basislinie für die Fehlerberechnung wieder aufaddiert
	m_Zeilendichte = Flaeche;
	// Funktion double Messung_Limb::Evaluate_Error_primitive(double* x,
	// double* y, double a0,double a1, double A, double FWHM, double x0,
	// int Anfangsindex, int Endindex)
	m_Fehler_Zeilendichten =
		Evaluate_Error_primitive(Vollfenster_WL,
								 Messwerte_Quotient, a0, a1, m_Zeilendichte,
								 Spezfenst.m_FWHM,
								 Spezfenst.m_Wellenlaengen[Index], 0,
								 N_Vollfenster - 1);
	//////////////////////////////////////////////
	/////////////////////////////////////////////


	//Speicher freimachen
	delete[] Basisfenster_WL;
	delete[] Basisfenster_Intensitaet;
	delete[] Vollfenster_WL;
	delete[] Vollfenster_Limb;
	delete[] Vollfenster_Sonne;
	delete[] Vollfenster_Limb_abs_error;
	delete[] Messwerte_Quotient_error;
	delete[] Messwerte_Quotient_stabw;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Plots_der_Spektren_erzeugen
////////////////////////////////////////////////////////////////////////////////
int Messung_Limb::moving_average(int window_size)
{
	return my_moving_average(m_Intensitaeten, window_size);
}
int Messung_Limb::savitzky_golay(int window_size)
{
	return my_savitzky_golay(m_Intensitaeten, window_size);
}
int Messung_Limb::Intensitaeten_normieren(vector<double> &Sonnen_Intensitaet)
{
	//Teiler wurde vorher interpoliert
	//todo prüfen
	// der Teiler ist das interpolierte Sonnenspektrum
	for (int i = 0; i < m_Number_of_Wavelength; i++) {
		this->m_Intensitaeten_durch_piF[i]
			= this->m_Intensitaeten[i] / Sonnen_Intensitaet[i];
	}
	m_Sonne = Sonnen_Intensitaet;

	return 0;
}
//========================================
//========================================
int Messung_Limb::Intensitaeten_durch_piF_Gamma_berechnen(Speziesfenster Spezfenst, int Index)
{

	//Auf dem ganzen Fenster...Verschwendung !!!!!...
	for (int i = 0; i < m_Number_of_Wavelength; i++) { //langsam, optimierbar
		this->m_Intensitaeten_durch_piF_Gamma[i]
			= this->m_Intensitaeten_durch_piF[i]
			  / Spezfenst.m_Liniendaten[Index].m_Gamma;
	}
	return 0;
}
int Messung_Limb::Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen(Speziesfenster Spezfenst, int Index)
{

	//Auf dem ganzen Fenster...Verschwendung !!!!!

	// Wir berechnen den Gitterabstand nur einmal
	// Am besten gleich bei der Wellenlänge des Übergangs....
	// eigentlich reicht 0,11nm, falls es mal schneller gehn soll
	// die Gitterabstände sind aber über große Bereiche doch schon nicht linear
	int Ind = Get_Index(Spezfenst.m_Wellenlaengen[Index]);
	double Delta_WL = (m_Wellenlaengen[Ind + 1] - m_Wellenlaengen[Ind]);
	// Nun alles damit multiplizieren....wie gesagt..das ist etwas langsam,
	// da es sich um nen konstanten Faktor handelt
	for (int i = 0; i < m_Number_of_Wavelength; i++) { //langsam, optimierbar
		//m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[i]=m_Intensitaeten
		//_durch_piF_Gamma[i]*Delta_WL;
		// Delta_Wl ist in nm gegeben...
		// dann muss beim Peakfit nicht in nm umgerechnet werden
		// wenn integriert wird
		m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[i]
			= m_Intensitaeten_durch_piF_Gamma[i] / Delta_WL;
		// glaub man muss dividieren
	}
	return 0;
}
//========================================
//========================================
int Messung_Limb::Deklinationswinkel_bestimmen()
{
	const double pi = M_PI;
	// Formel nach der englischen Wikipedia
	//theta=-23,45*cos(360° *(N+10)/365);
	// dieser Winkel ändert sich nicht sehr stark von Tag zu Tag
	int Monatstage[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	// reicht auch auf Tagesgenauigkeit
	double Tage = 0;
	for (int i = 0; i < (this->m_Monat - 1); i++) {
		Tage += Monatstage[i];
	}
	Tage += this->m_Tag;
	Tage += (double) this->m_Stunde / 24.0;
	//cout<<Tage<<"\n";
	//double bla=cos(360.0/365.0*(Tage+10.0)*pi/180.0);
	//cout<<bla<<"\n";
	this->m_Deklinationswinkel = -23.45 * cos(360.0 / 365.0 * (Tage + 10.0)
								* pi / 180.0);
	return 0;
}// int        Deklinationswinkel_bestimmen() ende

//========================================
//========================================
// Funktionsstart  Sonnen_Longitude_bestimmen
int Messung_Limb::Sonnen_Longitude_bestimmen()
{
	// 12 Uhr Mittags (UTC) ist die Sonne bei Phi=0
	// (glaub Greenwich, oder zumindest in etwa) im Zenit
	double Stunden = 0.0;
	Stunden += this->m_Stunde;
	Stunden += (double) this->m_Minute / 60.0;

	this->m_Sonnen_Longitude = 180.0 - 360.0 * (Stunden / 24.0);

	return 0;
}
//ENDE Sonnen_Longitude_bestimmen
//========================================
//========================================
Ausgewertete_Messung_Limb Messung_Limb::Ergebnis_Zusammenfassen()
{
	Ausgewertete_Messung_Limb aus;
	//Ergebnisse
	aus.m_Zeilendichte = this->m_Zeilendichte;
	aus.m_Fehler_Zeilendichten = this->m_Fehler_Zeilendichten;
	//Zwischenergebnisse
	aus.m_Deklination = this->m_Deklinationswinkel;
	aus.m_Sonnen_Longitude = this->m_Sonnen_Longitude;
	aus.m_Wellenlaenge = 0;
	// Nullinitialisierung...
	// die Wellenlänge des Übergangs steckt im Speziesfenster
	//Datum
	aus.m_Jahr = this->m_Jahr;
	aus.m_Monat = this->m_Monat;
	aus.m_Tag = this->m_Tag;
	aus.m_Stunde = this->m_Stunde;
	aus.m_Minute = this->m_Minute;
	// Geolocation
	aus.m_Lattidude_Sat = this->m_Lattidude_Sat;
	aus.m_Longitude_Sat = this->m_Longitude_Sat;
	aus.m_Hoehe_Sat = this->m_Hoehe_Sat;
	aus.m_Lattidude_TP = this->m_Lattidude_TP;
	aus.m_Longitude_TP = this->m_Longitude_TP;
	aus.m_Hoehe_TP = this->m_Hoehe_TP;
	aus.m_Erdradius = this->m_Erdradius;
	return aus;
}//Ausgewertete_Messung_Limb Ergebnis_Zusammenfassen() ende
//========================================

//========================================
//Methoden ende
//========================================

////////////////////////////////////////////////////////////////////////////////
//
//Hilfsfunktionen
//
////////////////////////////////////////////////////////////////////////////////
int Messung_Limb::Get_Index(double WL)
{
	// Die Wellenlängen sind monoton steigend sortiert
	// Quicksearch algorithmus kann angewendet werden
	// Startelement-> 825/2 =412
	// Wellenlängen sind etwa 0,1 nm auseinander...also gebe Toleranz von 0,08nm
	// Maximale Schrittzahl ist ceil(log_2 N) also 10 Schritte
	bool gefunden = false;
	int unterer_Fensterindex = 0;
	int oberer_Fensterindex = 825;
	int startindex = 412;
	int aktueller_Index = startindex;

	if (abs(this->m_Wellenlaengen[oberer_Fensterindex] - WL) < 0.08) {
		return oberer_Fensterindex;
	}
	while (!gefunden) {
		//eventuell durch forschleife ersetzen,
		//weil Programm sich hier festhängen könnte
		// Wellenlänge gefunden
		if (abs(this->m_Wellenlaengen[aktueller_Index] - WL) < 0.08) {
			return aktueller_Index;// Schleifenabbruch und sofortige Rückgabe
		}
		if (m_Wellenlaengen[aktueller_Index] < WL) {
			unterer_Fensterindex = aktueller_Index;
		} else {
			oberer_Fensterindex = aktueller_Index;
		}//if m_WL<WL
		aktueller_Index = (oberer_Fensterindex + unterer_Fensterindex) / 2;
		//achtung hier wird abgerundet
	}//while
	return -1;
	// soweit sollte es eigentlich nicht kommen...
	// aber damit die Warnung verschwindet geben wir mal was zurück
}//ende int Messung_Limb::Get_Index(double WL)

int Messung_Limb::sb_Get_Index(double WL)
{
	vector<double>::iterator low;
	low = lower_bound(m_Wellenlaengen.begin(), m_Wellenlaengen.end(), WL);

	// catch edge cases
	if (low == m_Wellenlaengen.begin()) return 0;
	if (low == m_Wellenlaengen.end()) --low;

	return distance(m_Wellenlaengen.begin(), low) - 1;
}

int Messung_Limb::sb_Get_closest_index(double WL)
{
	int i = sb_Get_Index(WL);
	return (WL - m_Wellenlaengen[i]) < (m_Wellenlaengen[i + 1] - WL) ? i : i + 1;
}

void Messung_Limb::Fit_Linear(double *x, double *y, double &a0, double &a1,
		int Anfangsindex, int Endindex)
{
	//fit der Funktion y=a0+a1x;
	//Bestimmung von a und b im Intervall zwischen Anfangs und endindex
	a0 = 0;
	a1 = 0;
	int i;
	// benötigt werden die Mittelwerte von x,y,x*y,und x^2 =====================
	double xsum = 0.;
	double ysum = 0.;
	double xysum = 0.;
	double xxsum = 0.;
	for (i = Anfangsindex; i <= Endindex; i++) {
		xsum += x[i];
		ysum += y[i];
		xysum += x[i] * y[i];
		xxsum += x[i] * x[i];
	}
	//Mittelwerte
	double N = Endindex - Anfangsindex + 1.;
	double x_m = xsum / N;
	double y_m = ysum / N;
	double xy_m = xysum / N;
	double xx_m = xxsum / N;
	//==========================================================================
	// Parameter b
	a1 = (xy_m - y_m * x_m) / (xx_m - x_m * x_m);
	//Parameter a
	a0 = y_m - a1 * x_m;
}
void Messung_Limb::Fit_Linear(vector<double> &x, vector<double> &y, double &a0, double &a1,
		int Anfangsindex, int Endindex)
{
	//fit der Funktion y=a0+a1x;
	//Bestimmung von a und b im Intervall zwischen Anfangs und endindex
	a0 = 0.;
	a1 = 0.;
	int i;
	// benötigt werden die Mittelwerte von x,y,x*y,und x^2 =====================
	double xsum = 0.;
	double ysum = 0.;
	double xysum = 0.;
	double xxsum = 0.;
	for (i = Anfangsindex; i <= Endindex; i++) {
		xsum += x[i];
		ysum += y[i];
		xysum += x[i] * y[i];
		xxsum += x[i] * x[i];
	}
	//Mittelwerte
	double N = Endindex - Anfangsindex + 1.;
	double x_m = xsum / N;
	double y_m = ysum / N;
	double xy_m = xysum / N;
	double xx_m = xxsum / N;
	//==========================================================================
	// Parameter b
	a1 = (xy_m - y_m * x_m) / (xx_m - x_m * x_m);
	//Parameter a
	a0 = y_m - a1 * x_m;
}//Ende Fit_linear

void Messung_Limb::Fit_Polynom_4ten_Grades(double *x, double *y, double x0,
		double *Par_a0, double *Par_a1, double *Par_a2, double *Par_a3,
		double *Par_a4, int Anfangsindex, int Endindex)
{
	// Das geht auch analytisch, aber das ist eine laaaaaaaaaaaaaaaaange Formel,
	// deren Ableitung zwar trivial, aber Fehleranfällig ist(so vieeeeel zu
	// schreiben) Das Diagonalisieren des Linearen Gleichungssystems könnte man
	// auslagern,
	//sodass nur das Rückeinsetzen benutzt werden muss...  bei einer 5*6 Matrix
	//ist ist das aber vermutlich noch zu verschmerzen...wir werden sehn

	//Für jede Messung müssen 30 konstanten bestimmt werden, von denen aber
	//einige doppelt vorkommen

	double a0 = 0;
	double b0 = 0;
	double c0 = 0;
	double d0 = 0;
	double e0 = 0;
	double f0 = 0;
	/*b0*/              /*c0*/             /*d0*/             /*e0*/
	double e1 = 0;
	double f1 = 0;
	/*c0*/             /*d0*/             /*e0*/              /*e1*/
	double e2 = 0;
	double f2 = 0;
	/*d0*/             /*e0*/             /*e1*/              /*e2*/
	double e3 = 0;
	double f3 = 0;
	/*e0*/             /*e1*/            /*e2*/               /*e3*/
	double e4 = 0;
	double f4 = 0;
	// Nur 14 Parameter müssen echt bestimmt werden  //sieh 5x5 Matrix oben
	double h; //h wie hilfsvariable (soll nur ein Buchstabe sein)
	for (int i = Anfangsindex; i <= Endindex; i++) {
		h = x[i] - x0;
		a0 += 1;
		b0 += h;
		c0 += h * h;
		d0 += h * h * h;
		e0 += h * h * h * h;
		e1 += h * h * h * h * h;
		e2 += h * h * h * h * h * h;
		e3 += h * h * h * h * h * h * h;
		e4 += h * h * h * h * h * h * h * h;

		f0 += y[i];
		f1 += y[i] * h;
		f2 += y[i] * h * h;
		f3 += y[i] * h * h * h;
		f4 += y[i] * h * h * h * h;
	}
	// Zur Lösung des Gleichungssystems wird eine Lapack routine benutzt da
	// diese in Fortran geschrieben sind, muss die Transponierte Matrix
	// übergeben werden da die RHS nur aus einem Vektor besteht, gibt es keine
	// Verwirrung mit transponierten Matrix und Vektor aufbauen
	double LHS[25];
	double RHS[5];
	// LHS Matrix Spaltenweise eingeben
	LHS[0] = a0;
	LHS[5] = b0;
	LHS[10] = c0;
	LHS[15] = d0;
	LHS[20] = e0;
	LHS[1] = b0;
	LHS[6] = c0;
	LHS[11] = d0;
	LHS[16] = e0;
	LHS[21] = e1;
	LHS[2] = c0;
	LHS[7] = d0;
	LHS[12] = e0;
	LHS[17] = e1;
	LHS[22] = e2;
	LHS[3] = d0;
	LHS[8] = e0;
	LHS[13] = e1;
	LHS[18] = e2;
	LHS[23] = e3;
	LHS[4] = e0;
	LHS[9] = e1;
	LHS[14] = e2;
	LHS[19] = e3;
	LHS[24] = e4;

	RHS[0] = f0;
	RHS[1] = f1;
	RHS[2] = f2;
	RHS[3] = f3;
	RHS[4] = f4;
	// Restliche Vorbereitungen für Lapackroutine
	int N = 5;  //<-- Feldgröße Speed propto N^3 , LHS ist quadratisch, N ist Anzahl der Gitterpunkte
	int *IPIV;  //array mit der Pivotisierungsmatrix sollte so groß wie N sein, alle Elemente 0
	IPIV = new int[N];
	// ------ RHS oben definiert
	int NRHS = 1; //Spalten von RHS 1 nehmen, um keine c/Fortran Verwirrungen zu provozieren
	int LDA = N;
	int LDB = N;
	int INFO;
	//int Anzahl=N*N;  Anzahl sollte die Integer grenzen nicht überschreiten,
	//aber danbn sollte der Aufbau von LHS schon stören
	// AUFRUF   A ist LHS.transponiert und B ist RHS
	dgesv_(&N, &NRHS, LHS, &LDA, IPIV, RHS, &LDB, &INFO);
	// ENDE LU ZERLEGUNG
	delete[] IPIV;
	IPIV = 0;
	// Ergebnis steckt nun in RHS
	*Par_a0 = RHS[0];
	*Par_a1 = RHS[1];
	*Par_a2 = RHS[2];
	*Par_a3 = RHS[3];
	*Par_a4 = RHS[4];
	return;
}

void Messung_Limb::Fit_Peak_hyperbolic(double *x, double *y, double x0,
		double FWHM, double &A, int Anfangsindex, int Endindex)
{
	//Folgende Funktion ist fürs Integral über alles ordentlich auf 1 normiert
	//Spaltfunktion
	//cnorm           = 4.*PI*sqrt(2.) / FWHM**3      ! Normierung stimmt MPL
	//SlitFuncSPEC    = 1./( ( (.5*FWHM)**4 + X**4 ) * cnorm )
	//Im folgenden nenne ich SlitFuncSPEC==g
	//
	//Für einen Linearen Parameterfit der Funktion A*g gilt:
	// cih^2=sum(y-Ag)^2=sum(y^2-2Agy+A^2g^2)
	//dchi^2/dA=-2sum(gy)+2 A sum(g^2) das soll 0 sein
	//-> A=sum(gy)/sum(g^2)
	//
	// FWHM muss gegeben werden und wir werten die Funktion um den Mittelwert
	// x0 aus also statt X-> X-x0

	// In dieser Funktion wird die Fläche A der Spaltfunktion bestimmt, da die
	// Funktionwerte y=I/(piFGamma) sind so ist A dann die Säulendichte

	const double pi = M_PI;
	// Zahl der Messwertpaare
	// double, damit später keine Probleme beim weiterrechnen
	double sum_gy = 0.;
	double sum_gg = 0.;
	double g;
	const double cnorm = 4.0 * pi * M_SQRT2 / (FWHM * FWHM * FWHM); //lambda m
	// (0.5 * FWHM)^4
	const double fwhm2to4 = 0.0625 * FWHM * FWHM * FWHM * FWHM;

	for (int i = Anfangsindex; i <= Endindex; i++) {
		//g berechnen
		g = 1. / (cnorm * (fwhm2to4 + pow(x0 - x[i], 4)));
		//eine Rechnung...nicht  Zeitkritisch
		// sum_gy erhöhen
		sum_gy += g * y[i];
		// sum_gg erhöhen
		sum_gg += g * g;
	}
	A = sum_gy / sum_gg;
}
void Messung_Limb::Fit_Peak_hyperbolic(vector<double> &x, vector<double> &y,
		double x0, double FWHM, double &A, int Anfangsindex, int Endindex)
{
	//Folgende Funktion ist fürs Integral über alles ordentlich auf 1 normiert
	//Spaltfunktion
	//cnorm           = 4.*PI*sqrt(2.) / FWHM**3      ! Normierung stimmt MPL
	//SlitFuncSPEC    = 1./( ( (.5*FWHM)**4 + X**4 ) * cnorm )
	//Im folgenden nenne ich SlitFuncSPEC==g
	//
	//Für einen Linearen Parameterfit der Funktion A*g gilt:
	// cih^2=sum(y-Ag)^2=sum(y^2-2Agy+A^2g^2)
	//dchi^2/dA=-2sum(gy)+2 A sum(g^2) das soll 0 sein
	//-> A=sum(gy)/sum(g^2)
	//
	// FWHM muss gegeben werden und wir werten die Funktion um den Mittelwert
	// x0 aus also statt X-> X-x0

	// In dieser Funktion wird die Fläche A der Spaltfunktion bestimmt, da die
	// Funktionwerte y=I/(piFGamma) sind so ist A dann die Säulendichte

	const double pi = M_PI;
	// Zahl der Messwertpaare
	// double, damit später keine Probleme beim weiterrechnen
	double sum_gy = 0.;
	double sum_gg = 0.;
	double g;
	const double cnorm = 4.0 * pi * M_SQRT2 / (FWHM * FWHM * FWHM); //lambda m
	// (0.5 * FWHM)^4
	const double fwhm2to4 = 0.0625 * FWHM * FWHM * FWHM * FWHM;

	for (int i = Anfangsindex; i <= Endindex; i++) {
		//g berechnen
		g = 1. / (cnorm * (fwhm2to4 + pow(x0 - x[i], 4)));
		//eine Rechnung...nicht  Zeitkritisch
		// sum_gy erhöhen
		sum_gy += g * y[i];
		// sum_gg erhöhen
		sum_gg += g * g;
	}
	A = sum_gy / sum_gg;
}

double Messung_Limb::Evaluate_Error_primitive(double *x, double *y, double a0,
		double a1, double A, double FWHM, double x0, int Anfangsindex,
		int Endindex)
{
	/***************************************************************************
	Wie der Name schon sagt, ist dies eine eher einfache Berechnung des Fehlers.
	Summe der Quadratischen Abweichungen-> Chi^2 hmm nicht gut... aber als
	Wichtungsfaktor noch akzeptabel
	 **************************************************************************/
	double Error = 0;
	const double pi = M_PI;
	double cnorm = 4.0 * pi * M_SQRT2 / (FWHM * FWHM * FWHM);
	//double y_quadrat=0;
	for (int i = Anfangsindex; i < Endindex + 1; i++) {
		//Funktionswert Bestimmen
		double Basis = a0 + a1 * x[i];
		double Peak = A / (cnorm * (pow(0.5 * FWHM, 4) + pow(x0 - x[i], 4)));
		double Funktionswert = Peak + Basis;
		// Quadratische Abweichung des Funktionswerts zum Messwert Bestimmen
		// und aufaddieren
		Error += (Funktionswert - y[i]) * (Funktionswert - y[i]);
		//y_quadrat+=y[i]*y[i];
	}
	Error /= (Endindex - Anfangsindex + 1);
	Error = sqrt(Error);
	//Das ist nach Numerical Recipes der Fehlerbalken der Messpunkte Es ist
	//vermutlich anschaulicher diesen Fehler noch durch den Mittelwert zu
	//teilen;
	//y_quadrat/=(Endindex-Anfangsindex+1);
	//double quot=Error/sqrt(y_quadrat);
	//cout<<quot<<"\n";
	return Error;
}
double Messung_Limb::Evaluate_Error_primitive(vector<double> &x,
		vector<double> &y, double a0, double a1, double A, double FWHM,
		double x0, int Anfangsindex, int Endindex)
{
	/***************************************************************************
	Wie der Name schon sagt, ist dies eine eher einfache Berechnung des Fehlers.
	Summe der Quadratischen Abweichungen-> Chi^2 hmm nicht gut... aber als
	Wichtungsfaktor noch akzeptabel
	 **************************************************************************/
	double Error = 0;
	const double pi = M_PI;
	double cnorm = 4.0 * pi * M_SQRT2 / (FWHM * FWHM * FWHM);
	//double y_quadrat=0;
	for (int i = Anfangsindex; i < Endindex + 1; i++) {
		//Funktionswert Bestimmen
		double Basis = a0 + a1 * x[i];
		double Peak = A / (cnorm * (pow(0.5 * FWHM, 4) + pow(x0 - x[i], 4)));
		double Funktionswert = Peak + Basis;
		// Quadratische Abweichung des Funktionswerts zum Messwert Bestimmen
		// und aufaddieren
		Error += (Funktionswert - y[i]) * (Funktionswert - y[i]);
		//y_quadrat+=y[i]*y[i];
	}
	Error /= (Endindex - Anfangsindex + 1);
	Error = sqrt(Error);
	//Das ist nach Numerical Recipes der Fehlerbalken der Messpunkte Es ist
	//vermutlich anschaulicher diesen Fehler noch durch den Mittelwert zu
	//teilen;
	//y_quadrat/=(Endindex-Anfangsindex+1);
	//double quot=Error/sqrt(y_quadrat);
	//cout<<quot<<"\n";
	return Error;
}
////////////////////////////////////////////////////////////////////////////////
//
//Hilfsfunktionen ENDE
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen
//
////////////////////////////////////////////////////////////////////////////////
int Messung_Limb::Ausgabe_in_Datei(string Dateiname)
{
	//TODO hier kann sich später auch noch was verändern
	/************************************************************
	Diese Funktion gibt den Inhalt des Objekts in eine Datei aus.
	Damit kann überprüft werden:
	a) ob das einlesen aller parameter ordentlich geklappt hat
	b) ob die Unterfunktionen das richtige errechnet haben
	c) Die Felder können geplottet werden, zusammen mit den Fitfunktionen->
	   Überprüfung, ob Fit sinnvoll (z.b. mit Matlab oder Gnuplot)
	************************************************************/
	//Formatierte Ausgabe
	FILE *outfile;
	//Datei öffnen
	outfile = fopen(Dateiname.c_str(), "w");
	//checken, ob Datei auch offen fehlt...aber ok Funktion wird eh beim
	//debuggen eingesetzt...da kriegt man das schon raus..hoffentlich Datei
	//schreiben
	///////////////////////////////////////////////////////////
	//zunächst die randdaten in den header
	//fprintf(outfile,"blblblbllblblb");
	fprintf(outfile, "%20s %1.5E\n", "m_Zeilendichte[cm^-2]: ", m_Zeilendichte);
	fprintf(outfile, "%20s %1.5E\n", "Sigma für alle Messwerte[cm^-2]: ",
			m_Fehler_Zeilendichten);
	fprintf(outfile, "%20s %12s\n", "m_Dateiname_L1C: ",
			m_Dateiname_L1C.c_str()); //hoffentlich passt das auch
	fprintf(outfile, "%20s %5i\n", "m_Number_of_Wavelength: ",
			m_Number_of_Wavelength);
	fprintf(outfile, "%20s %5i\n", "m_Jahr: ", m_Jahr);
	fprintf(outfile, "%20s %5i\n", "m_Monat: ", m_Monat);
	fprintf(outfile, "%20s %5i\n", "m_Tag: ", m_Tag);
	fprintf(outfile, "%20s %5i\n", "m_Stunde: ", m_Stunde);
	fprintf(outfile, "%20s %5i\n", "m_Minute: ", m_Minute);
	fprintf(outfile, "%20s %1.5E\n", "m_Deklinationswinkel: ",
			m_Deklinationswinkel);
	fprintf(outfile, "%20s %1.5E\n", "m_Sonnen_Longitude: ",
			m_Sonnen_Longitude);
	//Geolokationen
	fprintf(outfile, "%20s %1.5E\n", "m_Lattidude_Sat: ", m_Lattidude_Sat);
	fprintf(outfile, "%20s %1.5E\n", "m_Longitude_Sat: ", m_Longitude_Sat);
	fprintf(outfile, "%20s %1.5E\n", "m_Hoehe_Sat: ", m_Hoehe_Sat);
	fprintf(outfile, "%20s %1.5E\n", "m_Lattidude_TP: ", m_Lattidude_TP);
	fprintf(outfile, "%20s %1.5E\n", "m_Longitude_TP: ", m_Longitude_TP);
	fprintf(outfile, "%20s %1.5E\n", "m_Hoehe_TP: ", m_Hoehe_TP);
	fprintf(outfile, "%20s %1.5E\n", "m_Erdradius: ", m_Erdradius);
	//Überschrift
	fprintf(outfile, "%12s %12s %12s %12s %12s\n",
			"m_Wellenlaengen", "m_Intensitaeten", "m_Intensitaeten_durch_piF",
			"m_Intensitaeten_durch_piF_Gamma",
			"m_Intensitaeten_durch_piF_Gamma*ds");
	//for(int i=0;i<m_Number_of_Wavelength;i++)
	for (int i = 0; i < m_Number_of_Wavelength; i++) {
		//die letzte Zeile der Datei ist leer, da \n in der Vorletzten steht
		fprintf(outfile, "%1.5E %1.5E %1.5E %1.5E %1.5E\n",
				m_Wellenlaengen[i], m_Intensitaeten[i],
				m_Intensitaeten_durch_piF[i],
				m_Intensitaeten_durch_piF_Gamma[i],
				m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[i]);
	}
	///////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile);
	return 0;
}//Ausgabe_in_Datei ENDE
////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen ENDE
//
////////////////////////////////////////////////////////////////////////////////

