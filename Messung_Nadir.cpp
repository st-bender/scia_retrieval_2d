/*
 * Messung_Nadir.cpp
 *
 *  Created on: 30.04.2010
 *      Author: martin
 */
//Macro SAVEDELETE
#define SAVEDELETE(array)     \
if(array!=0)                                 \
{                                                \
    delete[] array;                        \
    array=0;                                 \
}                                                \
//Das Makro ist nur Save, falls die arrays alle mit 0 initialisiert werden

#include <stdio.h>
#include "Messung_Nadir.h"
#include <cmath>
#include "stdlib.h"
#include "Ausgewertete_Messung_Nadir.h"
#include <stdio.h>
#include <iostream>
#include "Ausdrucke.h"

using namespace std;
//////////////////////////////////////////////////
//constructor
//////////////////////////////////////////////////
Messung_Nadir::Messung_Nadir()
{
	//Pointer als NULL-Pointer initialisieren
	m_Wellenlaengen = 0;
	m_Intensitaeten = 0;
	m_Intensitaeten_relativer_Fehler = 0;
	m_Intensitaeten_durch_piF = 0;
	m_Intensitaeten_durch_piF_Gamma = 0;
}
//////////////////////////////////////////////////
//constructor ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
// copyconstructor
//////////////////////////////////////////////////
Messung_Nadir::Messung_Nadir(const Messung_Nadir &rhs)
{
	//Pointer als NULL-Pointer initialisieren
	m_Wellenlaengen = 0;
	m_Intensitaeten = 0;
	m_Intensitaeten_relativer_Fehler = 0;
	m_Intensitaeten_durch_piF = 0;
	m_Intensitaeten_durch_piF_Gamma = 0;
	//Restproblem bei Assignement Operator gelöst
	*this = rhs;
}
//////////////////////////////////////////////////
// copyconstructor ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Destructor
//////////////////////////////////////////////////
Messung_Nadir::~Messung_Nadir()
{
	save_delete_all_memory();
}
//////////////////////////////////////////////////
//Destructor ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Assignmentoperator Overload
//////////////////////////////////////////////////
Messung_Nadir &Messung_Nadir::operator =(const Messung_Nadir &rhs)
{
	// Prevent self assignment. a=a;  We say two Strings
	// are equal if their memory addresses are equal.
	if (this == &rhs)
		return *this;
	//cout<<"skalare kopieren\n";
	//Ergebnisse
	m_Zeilendichte                 = rhs.m_Zeilendichte;
	m_Fehler_Zeilendichten     = rhs.m_Fehler_Zeilendichten;
	//Zwischenergebnisse
	m_Deklinationswinkel        = rhs.m_Deklinationswinkel;
	m_Sonnen_Longitude        = rhs.m_Sonnen_Longitude;
	// Herkunftsmerkmale
	m_Dateiname_L1C            = rhs.m_Dateiname_L1C;
	m_Messung_ID                 = rhs.m_Messung_ID;
	// Dartum
	m_Jahr                            = rhs.m_Jahr;
	m_Monat                         = rhs.m_Monat;
	m_Tag                             = rhs.m_Tag;
	m_Stunde                        = rhs.m_Stunde;
	m_Minute                         = rhs.m_Minute;
	//Geolocations
	m_Lattitude_Sat               = rhs.m_Lattitude_Sat;
	m_Longitude_Sat             = rhs.m_Longitude_Sat;
	m_Hoehe_Sat                  = rhs.m_Hoehe_Sat;
	m_Lattitude_Ground         = rhs.m_Lattitude_Ground;
	m_Longitude_Ground        = rhs.m_Longitude_Ground;
	m_Erdradius                    = rhs.m_Erdradius;
	m_orbit_phase                 = rhs.m_orbit_phase;
	//für Füllbare Felder Wichtig
	m_Number_of_Wavelength = rhs.m_Number_of_Wavelength;
	////////////////////////////////////////////////////////////////////////
	//Füllbare Felder   Hier DEEPCOPY machen
	// TODO später nochmal prüfen, eigentlich wird hier 2*zuviel kopiert...evtl eine eigene copy routine für jedes Feld
	//cout<<"speicher freimachen\n";
	// Erstmal die Felder löschen
	save_delete_all_memory();
	// Felder neu anlegen
	//cout<<"speicher neu anmelden\n";
	m_Wellenlaengen                          = new double[this->m_Number_of_Wavelength];
	m_Intensitaeten                             = new double[this->m_Number_of_Wavelength];
	m_Intensitaeten_relativer_Fehler      = new double[this->m_Number_of_Wavelength];
	m_Intensitaeten_durch_piF             = new double[this->m_Number_of_Wavelength];
	m_Intensitaeten_durch_piF_Gamma = new double[this->m_Number_of_Wavelength];
	// Nun Füllbare Felder auffüllen -> Zeitintensiv
	//cout<<"Felder Fuellen\n";
	for (int i = 0; i < m_Number_of_Wavelength; i++) {
		m_Wellenlaengen[i]                           = rhs.m_Wellenlaengen[i];
		m_Intensitaeten[i]                              = rhs.m_Intensitaeten[i];
		m_Intensitaeten_relativer_Fehler[i]       = rhs.m_Intensitaeten_relativer_Fehler[i];
		m_Intensitaeten_durch_piF[i]              = rhs.m_Intensitaeten_durch_piF[i];
		m_Intensitaeten_durch_piF_Gamma[i] = rhs.m_Intensitaeten_durch_piF_Gamma[i];
	}
	////////////////////////////////////////////////////////////////////////
	//Füllbare Felder ENDE
	////////////////////////////////////////////////////////////////////////
	return *this;
}
//////////////////////////////////////////////////
// Assignmentoperator Overload ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Methoden
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////
//save_delete_all_memory()
//////////////////////////////////////////////////////////////////////////
void Messung_Nadir::save_delete_all_memory()// Speicher löschen
{
	// Die 4 Felder Sicher löschen
	SAVEDELETE(this->m_Wellenlaengen);
	SAVEDELETE(this->m_Intensitaeten);
	SAVEDELETE(this->m_Intensitaeten_relativer_Fehler);
	SAVEDELETE(this->m_Intensitaeten_durch_piF);
	SAVEDELETE(this->m_Intensitaeten_durch_piF_Gamma);
}
//////////////////////////////////////////////////
//save_delete_all_memory() ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Zeilendichte_Bestimmen
//////////////////////////////////////////////////
int Messung_Nadir::Zeilendichte_Bestimmen(Speziesfenster &Spezfenst, int Index, string Arbeitsverzeichnis, string mache_Fit_Plots, int MessungsNr)
{
	//kurz:
	//Diese Fitroutine ermittelt die Fläche von I/(piF*Gamma) bei der Übergangswellenlänge
	//Da nur wenige Punkte zu einer Linie beitragen, entspricht der Linie der Auflösungsfunktion
	//Diese ist ein Hyperboloide 4ter Ordnung

	//lang:
	// Die Umgebung jeder Linie wird in 3 Bereiche unterteilt. 2 sogenannte Basisfenster links und rechts von der Linie
	// und das Peakfenster selbst, welches in der Umsetzung auch Bereiche aus den beiden Basisfenstern beinhalten darf.
	// Bis Hierhin haben wir I/(piFGamma) berechnet. Nun wollen wir die "normierte" Intensität für einen Peak berechnen
	// Da der Peak gerade im Kanal zumeist auf einem relativ großen Untergrundsignal sitzt, wird dieses zunächst abgezogen.
	// Dafür wird im Bereich um die Linie herum, wo keine signifikanten anderen Linien liegen(Basisfenster) ein linearer Fit
	// der Grundlinie durchgeführt, welcher vom der "normierten" Intensität abgezogen wird.
	// Danach wird im Bereich des Peakfensters die Fläche des Peaks bestimmt.
	// Da die Linie zumeist nur aus 3 Punkten besteht,  hat das Profil im wesentlichen die Form der Auflösungsfunktion, welche
	// eine hyperboloide 4ten Grades ist. Diese hat im wesentlichen 3 Parameter: Wellenlänge des Peaks, Breite und Fläche.
	// Die robusteste Variante ist es, die Wellenlänge und die Halbwertsbreite der Linie vorzugeben und nur die Fläche anzufitten.
	// Diese Variante kann durch einen linearen Parameterfit erreicht werden und ist zum einen schnell und was hier viel wichtiger ist
	// robust. Andere nichtlineare Verfahren sind entweder nicht robust(Gauss Newton, insbesondere bei dem
	// schlechten Signal/Noise Verhältnis) oder deutlich langsamer (simulated annealing bzw. ist das aufwändiger dies noch schnel
	//l zu implementieren... das wäre Schritt 2)
	// Die ermittelte Fläche ist dann unsere gesuchte Säulenzeilendichte.

	// I/(piFGamma)=integral(AMF n ds) mit AMF = s exp(-tau) ...aber zu der Formel später nochmal zurück
	// Das spätere Retrieval ermittelt dann die Dichte n aus der rechten Seite

	double *Basisfenster_WL;
	double *Basisfenster_Intensitaet;
	double *Peakfenster_WL;
	double *Peakfenster_Intensitaet;
	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min   = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max  = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
	int Index_Peakfenster_min             = Get_Index(Spezfenst.m_Peakfenster_WLmin[Index]);
	int Index_Peakfenster_max            = Get_Index(Spezfenst.m_Peakfenster_WLmax[Index]);
	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int Speicherbedarf_Basis = Bas_l + Bas_r;
	int Speicherbedarf_Peak = Index_Peakfenster_max - Index_Peakfenster_min + 1;
	// Speicher anfordern
	Basisfenster_WL          = new double[Speicherbedarf_Basis];
	Basisfenster_Intensitaet = new double[Speicherbedarf_Basis];
	Peakfenster_WL           = new double[Speicherbedarf_Peak];
	Peakfenster_Intensitaet = new double[Speicherbedarf_Peak];
	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i]          = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] = this->m_Intensitaeten_durch_piF_Gamma[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i]          = this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] = this->m_Intensitaeten_durch_piF_Gamma[Index_Basisfenster_rechts_min + i];
	}
	//Peakfenster WL und I auffüllen
	for (int i = 0; i < Speicherbedarf_Peak; i++) {
		Peakfenster_WL[i]          = m_Wellenlaengen[Index_Peakfenster_min + i];
		Peakfenster_Intensitaet[i] = m_Intensitaeten_durch_piF_Gamma[Index_Peakfenster_min + i];
	}
	// linearen Fit des Basisfensters durchführen
	// Proto: Fit_Linear(double* x,double* y, double& a0, double& a1,int Anfangsindex, int Endindex)
	double a0, a1;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0,  a1, 0, Speicherbedarf_Basis - 1);
	// lineare Funktion von Intensitäten des Peakfenster abziehen
	for (int i = 0; i < Speicherbedarf_Peak; i++) {
		Peakfenster_Intensitaet[i] -= a0 + a1 * Peakfenster_WL[i];
	}
	// Hyperboloiden an Peakfenster anfitten
	// Proto: Fit_Peak_hyperbolic(double* x,double* y,double x0, double FWHM, double& A, int Anfangsindex, int Endindex)
	double Flaeche;
	Fit_Peak_hyperbolic(Peakfenster_WL, Peakfenster_Intensitaet, Spezfenst.m_Wellenlaengen[Index],
						Spezfenst.m_FWHM, Flaeche, 0, Speicherbedarf_Peak - 1);
	//Fehler des Fits bestimmen... da Peakfenster_Intensitaet nicht mehr gebraucht wird die Basislinie für die Fehlerberechnung wieder
	//aufaddiert
	m_Zeilendichte = Flaeche;
	for (int i = 0; i < Speicherbedarf_Peak; i++) {
		Peakfenster_Intensitaet[i] += a0 + a1 * Peakfenster_WL[i];
	}
	// Funktion double Messung_Limb::Evaluate_Error_primitive(double* x,double* y, double a0,double a1, double A, double FWHM, double x0, int Anfangsindex, int Endindex)
	m_Fehler_Zeilendichten = Evaluate_Error_primitive(Peakfenster_WL, Peakfenster_Intensitaet, a0, a1, m_Zeilendichte, Spezfenst.m_FWHM, Spezfenst.m_Wellenlaengen[Index], 0, Speicherbedarf_Peak - 1);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Hier kann man zur Testzwecken noch einen Plot machen  /////////////////////////////////////////
	if (mache_Fit_Plots == "ja") {
		//TODO das als Funktion implementieren
		double *Funktion;
		Funktion           = new double[Speicherbedarf_Peak];
		const double pi = 3.14159265;
		double FWHM  = Spezfenst.m_FWHM;
		double cnorm   =  4.0 * pi * sqrt(2.0) / (FWHM * FWHM * FWHM);

		for (int i = 0; i < Speicherbedarf_Peak; i++) {
			double Basis = a0 + a1 * Peakfenster_WL[i];
			double Peak = m_Zeilendichte / (cnorm * (pow(0.5 * FWHM, 4) + pow(Peakfenster_WL[i] - Spezfenst.m_Wellenlaengen[Index], 4)));
			//cout<<Peak<<"\n";
			//cout<<m_Zeilendichte<<"\n";
			Funktion[i] = Peak + Basis;
		}

		string s1, s_OrbNum, s2;
		char buf[256];
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...falls / im Namen ist das schlecht
		string Datnam = m_Dateiname_L1C.substr(m_Dateiname_L1C.size() - 39, 39);

		//TODO Pfad anpassen
		sprintf(buf, "mkdir %s/Plots 2>/dev/null", Arbeitsverzeichnis.c_str());
		string Befehl = buf;
		system(Befehl.c_str());
		sprintf(buf, "%s_%s_%i_%i.ps", Datnam.c_str(), Spezfenst.m_Spezies_Name.c_str(), MessungsNr, Index);
		string new_datnam = buf;
		sprintf(buf, "%s/Plots/%s", Arbeitsverzeichnis.c_str(), new_datnam.c_str());
		s1 = buf; //s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern, um später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal, wie die Datei heißt... die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Nadirdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		sprintf(buf, "Orbit %5s Nadir GP: Lat: %G G(rad)  Lon: %G G Sat: Lat: %G G Lon: %G",
				s_OrbNum.c_str(), m_Lattitude_Ground, m_Longitude_Ground, m_Lattitude_Sat, m_Longitude_Sat);
		s2 = buf;
		//cout<<s1<<"\n";
		//int Plot_2xy(string Dateiname,string title, string xlabel, string ylabel,double* x1,double*y1, double* x2,double* y2,int Startindex,int Endindex);
		//Plot_2xy(s1.c_str(),s1.substr(s1.size()-50,50).c_str(),"$\\lambda$ in nm","$\\frac{I}{\\piF\\gamma}$",Peakfenster_WL,Peakfenster_Intensitaet,Peakfenster_WL,Funktion,0,Speicherbedarf_Peak-1); //-> Fit geht
		Plot_2xy(Arbeitsverzeichnis.c_str(), s1.c_str(), s2.c_str(), "Wellenlaenge in nm", "Schraege Saeule bei Peakposition in cm^{-2}/nm",
				 Peakfenster_WL, Peakfenster_Intensitaet, Peakfenster_WL, Funktion, 0, Speicherbedarf_Peak - 1,
				 m_Zeilendichte, m_Fehler_Zeilendichten);
		SAVEDELETE(Funktion);
	}
	// Ende Plot //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//dynamische Felder der Funktion löschen
	SAVEDELETE(Basisfenster_WL);
	SAVEDELETE(Basisfenster_Intensitaet);
	SAVEDELETE(Peakfenster_WL);
	SAVEDELETE(Peakfenster_Intensitaet);
	return 0;
}//int Zeilendichte_Bestimmen() ende
//////////////////////////////////////////////////
//Zeilendichte_Bestimmen ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Intensitaeten_normieren
//////////////////////////////////////////////////
int Messung_Nadir::Intensitaeten_normieren(double *Teiler) // Da das Limbspektrum mehr Wellenlängen hat, ist das kein Problem
{
	for (int i = 0; i < this->m_Number_of_Wavelength; i++) {
		this->m_Intensitaeten_durch_piF[i] = this->m_Intensitaeten[i] / Teiler[i];
	}
	return 0;
}
//////////////////////////////////////////////////
//Intensitaeten_normieren ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Intensitaeten_durch_piF_Gamma_berechnen
//////////////////////////////////////////////////
int Messung_Nadir::Intensitaeten_durch_piF_Gamma_berechnen(Speziesfenster Spezfenst, int Index)
{
	for (int i = 0; i < this->m_Number_of_Wavelength; i++) { //langsam, optimierbar
		this->m_Intensitaeten_durch_piF_Gamma[i] = this->m_Intensitaeten_durch_piF[i] / Spezfenst.m_Liniendaten[Index].m_Gamma * 0.11; // Faktor 0.11 siehe Diskussion über Emissivitäten LimbLimb
	}
	//cout<<"m_Intensitaeten_durch_piF[0] "<<m_Intensitaeten_durch_piF[0]<<"\tgamma "<<Spezfenst.m_Liniendaten[Index].m_Gamma<<"\n";
	return 0;
}
//////////////////////////////////////////////////
//Intensitaeten_durch_piF_Gamma_berechnen ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Deklinationswinkel_bestimmen
//////////////////////////////////////////////////
int Messung_Nadir::Deklinationswinkel_bestimmen() // auch gleich wie bei Limb
{
	const double pi = 3.14159265;
	// Formel nach der englischen Wikipedia
	//theta=-23,45*cos(360° *(N+10)/365);
	// dieser Winkel ändert sich nicht sehr stark von Tag zu Tag
	int Monatstage[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // reicht auch auf Tagesgenauigkeit
	double Tage = 0;
	for (int i = 0; i < (this->m_Monat - 1); i++) {
		Tage += Monatstage[i];
	}
	Tage += this->m_Tag;
	Tage += (double) this->m_Stunde / 24.0;
	//double bla=cos(360.0/365.0*(Tage+10.0)*pi/180.0);
	//cout<<bla<<"\n";
	this->m_Deklinationswinkel = -23.45 * cos((double)360 / 365 * (Tage + 10) * pi / 180);
	return 0;
}// int        Deklinationswinkel_bestimmen() ende
//////////////////////////////////////////////////
//Deklinationswinkel_bestimmen ENDE
//////////////////////////////////////////////////
//========================================
// Funktionsstart  Sonnen_Longitude_bestimmen
int Messung_Nadir::Sonnen_Longitude_bestimmen()
{
	// 12 Uhr Mittags (UTC) ist die Sonne bei Phi=0(glaub Greenwich, oder zumindest in etwa) im Zenit
	double Stunden = 0.0;
	Stunden += this->m_Stunde;
	Stunden += (double) this->m_Minute / ((double) 60.0);

	this->m_Sonnen_Longitude = 180.0 - 360.0 * (Stunden / 24.0);

	return 0;
}
//ENDE Sonnen_Longitude_bestimmen
//========================================
//////////////////////////////////////////////////
//Ergebnis_Zusammenfassen
//////////////////////////////////////////////////
Ausgewertete_Messung_Nadir  Messung_Nadir::Ergebnis_Zusammenfassen()
{
	Ausgewertete_Messung_Nadir aus;
	aus.m_Wellenlaenge = 0; // Nullinitialisierung...die Wellenlänge des Übergangs steckt im Speziesfenster
	//Ergebnisse
	aus.m_Zeilendichte = this->m_Zeilendichte;
	aus.m_Fehler_Zeilendichten = this->m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	aus.m_Deklination = this->m_Deklinationswinkel;
	aus.m_Sonnen_Longitude = this->m_Sonnen_Longitude;
	// Datum
	aus.m_Jahr = this->m_Jahr;
	aus.m_Monat = this->m_Monat;
	aus.m_Tag = this->m_Tag;
	aus.m_Stunde = this->m_Stunde;
	aus.m_Minute = this->m_Minute;
	// Geolokationen
	aus.m_Lattitude_Sat = this->m_Lattitude_Sat;
	aus.m_Longitude_Sat = this->m_Longitude_Sat;
	aus.m_Hoehe_Sat = this->m_Hoehe_Sat;
	aus.m_Lattitude_Ground = this->m_Lattitude_Ground;
	aus.m_Longitude_Ground = this->m_Longitude_Ground;
	aus.m_Erdradius = this->m_Erdradius;
	return aus;
}//Ausgewertete_Messung_Nadir  Ergebnis_Zusammenfassen() ENDE
//////////////////////////////////////////////////
//Ergebnis_Zusammenfassen ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Methoden ENDE
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Hilfsfunktionen  //nur intern aufrufen!!!
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////
//Get_Index
//////////////////////////////////////////////////
int Messung_Nadir::Get_Index(double WL)
{
	// Die Wellenlängen sind monoton steigend sortiert
	// Quicksearch algorithmus kann angewendet werden
	// Startelement-> 826/2 =413
	// Wellenlängen sind etwa 0,1 nm auseinander...also gebe Toleranz von 0,08nm
	// Maximale Schrittzahl ist ceil(log_2 N) also 10 Schritte
	bool gefunden = false;
	int unterer_Fensterindex = 0;
	int oberer_Fensterindex = this->m_Number_of_Wavelength - 1;
	int startindex = oberer_Fensterindex / 2;
	int aktueller_Index = startindex;

	if (abs(this->m_Wellenlaengen[oberer_Fensterindex] - WL) < 0.08) {
		return oberer_Fensterindex;
	}
	while (!gefunden) { //eventuell durch forschleife ersetzen, weil Programm sich hier festhängen könnte
		// Wellenlänge gefunden
		if (abs(this->m_Wellenlaengen[aktueller_Index] - WL) < 0.08) {
			return aktueller_Index;// Schleifenabbruch und sofortige Rückgabe
		}
		if (m_Wellenlaengen[aktueller_Index] < WL) {
			unterer_Fensterindex = aktueller_Index;
		} else {
			oberer_Fensterindex = aktueller_Index;
		}//if m_WL<WL
		aktueller_Index = (oberer_Fensterindex + unterer_Fensterindex) / 2; //achtung hier wird abgerundet
	}//while
	return -1; // soweit sollte es eigentlich nicht kommen...aber damit die Warnung verschwindet geben wir mal was zurück
}//ende int Messung_Nadir::Get_Index(double WL)

//////////////////////////////////////////////////
//Get_Index ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Fit_Linear
//////////////////////////////////////////////////
void Messung_Nadir::Fit_Linear(double *x, double *y, double &a0, double &a1, int Anfangsindex, int Endindex)
{
	//fit der Funktion y=a0+a1x; Bestimmung von a und b im Intervall zwischen Anfangs und endindex
	a0 = 0;
	a1 = 0;
	int i;
	// benötigt werden die Mittelwerte von x,y,x*y,und x^2 ===============================
	double xsum = 0;
	double ysum = 0;
	double xysum = 0;
	double xxsum = 0;
	for (i = Anfangsindex; i <= Endindex; i++) {
		xsum += x[i];
		ysum += y[i];
		xysum += x[i] * y[i];
		xxsum += x[i] * x[i];
	}
	//Mittelwerte
	double x_m = xsum / (Endindex - Anfangsindex + 1);
	double y_m = ysum / (Endindex - Anfangsindex + 1);
	double xy_m = xysum / (Endindex - Anfangsindex + 1);
	double xx_m = xxsum / (Endindex - Anfangsindex + 1);
	//====================================================================================
	// Parameter b
	a1 = (xy_m - y_m * x_m) / (xx_m - x_m * x_m);
	//Parameter a
	a0 = y_m - a1 * x_m;
}
//////////////////////////////////////////////////
//Fit_Linear ENDE
/////////////////////////////////////////////////
//////////////////////////////////////////////////
//Fit_Peak_hyperbolic
///////////////////////////////////////////////////
void Messung_Nadir::Fit_Peak_hyperbolic(double *x, double *y, double x0, double FWHM, double &A, int Anfangsindex, int Endindex)
{
	//Folgende Funktion ist für das Integral über alles ordentlich auf 1 normiert
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
	// FWHM muss gegeben werden und wir werten die Funktion um den Mittelwert x0 aus also statt X-> X-x0

	// In dieser Funktion wird die Fläche A der Spaltfunktion bestimmt, da die Funktionwerte y=I/(piFGamma) sind
	// so ist A dann die Säulendichte

	const double pi = 3.14159265;
	double N = Endindex - Anfangsindex + 1; // Zahl der Messwertpaare //double, damit später keine Probleme beim weiterrechnen
	double sum_gy = 0;
	double sum_gg = 0;
	double g;
	double cnorm =  4.0 * pi * sqrt(2.0) / (FWHM * FWHM * FWHM); //lambda m

	for (int i = 0; i < N; i++) {
		//g berechnen
		g = 1 / (cnorm * (pow(0.5 * FWHM, 4) + pow(x0 - x[Anfangsindex + i], 4))); //eine Rechnung...nicht  Zeitkritisch
		// sum_gy erhöhen
		sum_gy += g * y[Anfangsindex + i];
		// sum_gg erhöhen
		sum_gg += g * g;
	}
	A = sum_gy / sum_gg;
}
//////////////////////////////////////////////////
//Fit_Peak_hyperbolic ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Evaluate_Error_primitive
//////////////////////////////////////////////////
double Messung_Nadir::Evaluate_Error_primitive(double *x, double *y, double a0, double a1, double A, double FWHM, double x0, int Anfangsindex, int Endindex)
{
	/*****************************************************************************************
	Wie der Name schon sagt, ist dies eine eher einfache Berechnung des Fehlers.
	Summe der Quadratischen Abweichungen-> Chi^2 hmm nicht gut... aber als Wichtungsfaktor noch akzeptabel
	******************************************************************************************/
	double Error = 0;
	const double pi = 3.14159265;
	double cnorm =  4.0 * pi * sqrt(2.0) / (FWHM * FWHM * FWHM);
	for (int i = Anfangsindex; i < Endindex + 1; i++) {
		//Funktionswert Bestimmen
		double Basis = a0 + a1 * x[i];
		double Peak = A / (cnorm * (pow(0.5 * FWHM, 4) + pow(x0 - x[i], 4)));
		double Funktionswert = Peak + Basis;
		//cout<<"Basis\t"<<Basis<<"\tPeak\t"<<Peak<<"\tFunktionswert\t"<<Funktionswert<<"\ty[i]\t"<<y[i]<<"\n";
		// Quadratische Abweichung des Funktionswerts zum Messwert Bestimmen und aufaddieren
		Error += (Funktionswert - y[i]) * (Funktionswert - y[i]);
	}
	Error /= (Endindex - Anfangsindex + 1);
	Error = sqrt(Error);
	return Error;
}
//////////////////////////////////////////////////
//Evaluate_Error_primitive ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Hilfsfunktionen  ENDE
//////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen
//
/////////////////////////////////////////////////////////////////////////////////
int Messung_Nadir::Ausgabe_in_Datei(string Dateiname)
{
	//TODO hier kann sich später auch noch was verändern
	/************************************************************
	Diese Funktion gibt den Inhalt des Objekts in eine Datei aus.
	Damit kann überprüft werden:
	a) ob das einlesen aller parameter ordentlich geklappt hat
	b) ob die Unterfunktionen das richtige errechnet haben
	c) Die Felder können geplottet werden, zusammen mit den Fitfunktionen-> Überprüfung, ob Fit sinnvoll (z.b. mit Matlab oder Gnuplot)
	************************************************************/
	//Formatierte Ausgabe
	FILE *outfile;
	//Datei öffnen
	outfile = fopen(Dateiname.c_str(), "w");
	//checken, ob Datei auch offen fehlt...aber ok Funktion wird eh beim debuggen eingesetzt...da kriegt man das schon raus..hoffentlich
	//Datei schreiben
	///////////////////////////////////////////////////////////
	//zunächst die randdaten in den header
	// Ergebnisse
	fprintf(outfile, "%12s %1.5E\n", "m_Zeilendichte: ", m_Zeilendichte);
	fprintf(outfile, "%12s %1.5E\n", "m_Fehler_Zeilendichten: ", m_Fehler_Zeilendichten);
	// Zwischenergebnisse
	fprintf(outfile, "%12s %1.5E\n", "m_Deklinationswinkel: ", m_Deklinationswinkel);
	fprintf(outfile, "%20s %1.5E\n", "m_Sonnen_Longitude: ", m_Sonnen_Longitude);
	// Herkunftsmerkmale
	fprintf(outfile, "%12s %20s\n", "m_Dateiname_L1C: ", m_Dateiname_L1C.c_str());
	fprintf(outfile, "%12s %i\n", "m_Messung_ID: ", m_Messung_ID);
	//Datum
	fprintf(outfile, "%12s %i\n", "m_Jahr: ", m_Jahr);
	fprintf(outfile, "%12s %i\n", "m_Monat: ", m_Monat);
	fprintf(outfile, "%12s %i\n", "m_Tag: ", m_Tag);
	fprintf(outfile, "%12s %i\n", "m_Stunde: ", m_Stunde);
	fprintf(outfile, "%12s %i\n", "m_Minute: ", m_Minute);
	// Geolokationen
	fprintf(outfile, "%12s %1.5E\n", "m_Lattitude_Sat: ", m_Lattitude_Sat);
	fprintf(outfile, "%12s %1.5E\n", "m_Longitude_Sat: ", m_Longitude_Sat);
	fprintf(outfile, "%12s %1.5E\n", "m_Hoehe_Sat: ", m_Hoehe_Sat);
	fprintf(outfile, "%12s %1.5E\n", "m_Lattitude_Ground: ", m_Lattitude_Ground);
	fprintf(outfile, "%12s %1.5E\n", "m_Longitude_Ground: ", m_Longitude_Ground);
	fprintf(outfile, "%12s %1.5E\n", "m_Erdradius: ", m_Erdradius);
	//Füllbare Felder
	fprintf(outfile, "%12s %i\n", "m_Number_of_Wavelength: ", m_Number_of_Wavelength);
	//Überschrift
	fprintf(outfile, "%12s %12s %12s %12s\n", "m_Wellenlaengen", "m_Intensitaeten", "m_Intensitaeten_durch_piF", "m_Intensitaeten_durch_piF_Gamma");
	//großes Feld
	for (int i = 0; i < m_Number_of_Wavelength; i++) {
		//die letzte Zeile der Datei ist so leer, das \n in der Vorletzten steht
		fprintf(outfile, "%1.5E %1.5E %1.5E %1.5E\n", m_Wellenlaengen[i], m_Intensitaeten[i], m_Intensitaeten_durch_piF[i], m_Intensitaeten_durch_piF_Gamma[i]);
	}
	///////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile);
	return 0;
}//Ausgabe_in_Datei ENDE
/////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen ENDE
//
/////////////////////////////////////////////////////////////////////////////////

