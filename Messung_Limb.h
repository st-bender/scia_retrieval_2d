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
#include "Ausgewertete_Messung_Limb.h"

class Messung_Limb
{
public:
	Messung_Limb();
	// copyconstructor
	Messung_Limb(const Messung_Limb &rhs);
	// Assignmentoperator Overload
	Messung_Limb &operator =(const Messung_Limb &rhs);
	//Methoden
	int Zeilendichte_Bestimmen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots);
	int Saeulendichte_Bestimmen_MgI285nm(class Speziesfenster &Spezfenst,
			int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots,
			double *mean_10_20);
	int Plots_der_Spektren_erzeugen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots,
			double *mean_10_20);
	int slant_column_NO(class NO_emiss &NO, std::string mache_Fit_Plots,
			class Sonnenspektrum &sol_spec, int idx,
			class Speziesfenster &Spezfenst, std::string Arbeitsverzeichnis);
	double fit_NO_spec(class NO_emiss &NO, std::vector<double> &x,
			std::vector<double> &y, double &rms_err);
	int moving_average(int window_size);
	int savitzky_golay(int window_size);
	double msise_temperature();
	int Intensitaeten_normieren(std::vector<double> &Sonnen_Intensitaet);
	//int        Intensitaeten_normieren(Sonnenspektrum Solspec, Fenster);
	// hier müsste überlegt werden, wie man mehrfachkorrekturen umgeht
	int Intensitaeten_durch_piF_Gamma_berechnen(class Speziesfenster Spezfenst,
			double wl_gamma);
	// In der Formel ist piF in W/(m^2*Wellenlänge) verlangt..
	// also muss noch mit der Kanalbreite multipliziert werden
	int Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen(class Speziesfenster Spezfenst);
	int Deklinationswinkel_bestimmen();
	int Sonnen_Longitude_bestimmen();

	Ausgewertete_Messung_Limb Ergebnis_Zusammenfassen();

	//Hilfsfunktionen  //nur intern aufrufen!!!
	int Get_Index(double WL);
	int sb_Get_Index(double WL);
	int sb_Get_closest_index(double WL);
	void Fit_Linear(double *x, double *y, double &a0, double &a1,
			double &rms_err, int Anfangsindex, int Endindex);
	void Fit_Polynom_4ten_Grades(double *x, double *y, double x0,
			double *Par_a0, double *Par_a1, double *Par_a2, double *Par_a3,
			double *Par_a4, int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(double *x, double *y, double x0, double FWHM,
			double &A, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(double *x, double *y, double a0, double a1,
			double A, double FWHM, double x0, int Anfangsindex, int Endindex);
	// vectorised prototypes
	void Fit_Linear(std::vector<double> &x, std::vector<double> &y,
			double &a0, double &a1, double &rms_err,
			int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(std::vector<double> &x, std::vector<double> &y, double x0,
			double FWHM, double &A, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(std::vector<double> &x, std::vector<double> &y,
			double a0, double a1, double A, double FWHM, double x0,
			int Anfangsindex, int Endindex);
	//Wartungsfunktionen
	//können und sollten sogar später auskommentiert werden,
	//und dienen im wesentlichen zum debuggen
	//zum Testen und debuggen und überprüfen, ob der fit halbwegs passt
	int Ausgabe_in_Datei(std::string Dateiname);

	//Membervariablen

	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	// total number density at measurement point
	double total_number_density;
	// Zwischenergebnisse
	double m_Deklinationswinkel;  // aus Datum berechenbar
	double m_Sonnen_Longitude;
	//double  m_Streuwinkel;  // ist mehrere Streuwinkel auf einer Messachse...
							  //muss später berechnet werden
	//double  m_SZA;          //abhängig von Deklination Uhrzeit und Lat und Lon
							  //.....dasselbe wie beim streuwinkel
	//Dateiname
	std::string m_Dateiname_L1C;
	//Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	int m_Sekunde;
	// Geolocation
	double m_Latitude_Sat;   // hierfür muss das Programm noch geändert werden
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Latitude_TP;
	double m_Longitude_TP;
	double m_Hoehe_TP;
	double m_Erdradius;

	double m_orbit_phase;
	double m_TP_SZA;      // alt
	//double  m_SAA_TP;     // alt
	// Datenfelder
	int m_Number_of_Wavelength;
	std::vector<double> m_Wellenlaengen;
	std::vector<double> m_Sonne;
	std::vector<double> m_Intensitaeten;  // genauer genommen photonen/(s cm^2nm)
	std::vector<double> m_Intensitaeten_relativer_Fehler; // 1=100% vom Messwert
	std::vector<double> m_Intensitaeten_durch_piF;
	std::vector<double> m_Intensitaeten_durch_piF_Gamma; // ein bisschen mehr Speicher...optimierbar
	std::vector<double> m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand;
};

double slit_func(double fwhm, double x0, double x);
double slit_func_gauss(double fwhm, double x0, double x);
#endif /* MESSUNG_LIMB_HH_ */
