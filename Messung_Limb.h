/*
 * Messung_Limb.h
 *
 *  Created on: 12.04.2010
 *      Author: martin langowski
 */
#ifndef MESSUNG_LIMB_HH_
#define MESSUNG_LIMB_HH_

#include <vector>
#include <string>
#include "Speziesfenster.h"
#include "Ausgewertete_Messung_Limb.h"
#include "Konfiguration.h"


using namespace std;

class Messung_Limb
{
public:
	Messung_Limb();
	// copyconstructor
	Messung_Limb(const Messung_Limb &rhs);
	// Assignmentoperator Overload
	Messung_Limb &operator =(const Messung_Limb &rhs);
	//Methoden
	int Zeilendichte_Bestimmen(Speziesfenster &Spezfenst, int Index,
			string Arbeitsverzeichnis, string mache_Fit_Plots);
	int Saeulendichte_Bestimmen_MgI285nm(Speziesfenster &Spezfenst, int Index,
			string Arbeitsverzeichnis, string mache_Fit_Plots,
			double *mean_10_20);
	int Plots_der_Spektren_erzeugen(Speziesfenster &Spezfenst, int Index,
			string Arbeitsverzeichnis, string mache_Fit_Plots,
			double *mean_10_20);
	int moving_average(int window_size);
	int savitzky_golay(int window_size);
	int Intensitaeten_normieren(vector<double> &Sonnen_Intensitaet);
	//int        Intensitaeten_normieren(Sonnenspektrum Solspec, Fenster);
	// hier müsste überlegt werden, wie man mehrfachkorrekturen umgeht
	int Intensitaeten_durch_piF_Gamma_berechnen(Speziesfenster Spezfenst,
			int Index);
	// In der Formel ist piF in W/(m^2*Wellenlänge) verlangt..
	// also muss noch mit der Kanalbreite multipliziert werden
	int Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen(Speziesfenster Spezfenst, int Index);
	int Deklinationswinkel_bestimmen();
	int Sonnen_Longitude_bestimmen();

	Ausgewertete_Messung_Limb Ergebnis_Zusammenfassen();

	//Hilfsfunktionen  //nur intern aufrufen!!!
	int Get_Index(double WL);
	int sb_Get_Index(double WL);
	int sb_Get_closest_index(double WL);
	void Fit_Linear(double *x, double *y, double &a0, double &a1,
			int Anfangsindex, int Endindex);
	void Fit_Polynom_4ten_Grades(double *x, double *y, double x0,
			double *Par_a0, double *Par_a1, double *Par_a2, double *Par_a3,
			double *Par_a4, int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(double *x, double *y, double x0, double FWHM,
			double &A, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(double *x, double *y, double a0, double a1,
			double A, double FWHM, double x0, int Anfangsindex, int Endindex);
	// vectorised prototypes
	void Fit_Linear(vector<double> &x, vector<double> &y, double &a0, double &a1,
			int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(vector<double> &x, vector<double> &y, double x0,
			double FWHM, double &A, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(vector<double> &x, vector<double> &y,
			double a0, double a1, double A, double FWHM, double x0,
			int Anfangsindex, int Endindex);
	//Wartungsfunktionen
	//können und sollten sogar später auskommentiert werden,
	//und dienen im wesentlichen zum debuggen
	//zum Testen und debuggen und überprüfen, ob der fit halbwegs passt
	int Ausgabe_in_Datei(string Dateiname);

	//Membervariablen

	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	double m_Deklinationswinkel;  // aus Datum berechenbar
	double m_Sonnen_Longitude;
	//double  m_Streuwinkel;  // ist mehrere Streuwinkel auf einer Messachse...
							  //muss später berechnet werden
	//double  m_SZA;          //abhängig von Deklination Uhrzeit und Lat und Lon
							  //.....dasselbe wie beim streuwinkel
	//Dateiname
	string m_Dateiname_L1C;
	//Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	// Geolocation
	double m_Lattidude_Sat;   // hierfür muss das Programm noch geändert werden
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Lattidude_TP;
	double m_Longitude_TP;
	double m_Hoehe_TP;
	double m_Erdradius;

	double m_orbit_phase;
	double m_TP_SZA;      // alt
	//double  m_SAA_TP;     // alt
	// Datenfelder
	int m_Number_of_Wavelength;
	vector<double> m_Wellenlaengen;
	vector<double> m_Sonne;
	vector<double> m_Intensitaeten;  // genauer genommen photonen/(s cm^2nm)
	vector<double> m_Intensitaeten_relativer_Fehler; // 1=100% vom Messwert
	vector<double> m_Intensitaeten_durch_piF;
	vector<double> m_Intensitaeten_durch_piF_Gamma; // ein bisschen mehr Speicher...optimierbar
	vector<double> m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand;
};
#endif /* MESSUNG_LIMB_HH_ */
