/*
 * Messung.h
 *
 * Copyright (c) 2015-2017 Stefan Bender
 *
 * Initial version created on: 27.01.2015
 *      Author: Stefan Bender
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
#ifndef MESSUNG_HH_
#define MESSUNG_HH_

#include <vector>
#include <string>

class Messung
{
public:
	explicit Messung(std::string filename = "dummy");
	// copyconstructor
	Messung(const Messung &rhs);
	// Assignmentoperator Overload
	Messung &operator =(const Messung &rhs);
	//Methoden
	void Intensitaeten_normieren(std::vector<double> &Sonnen_Intensitaet);
	//int        Intensitaeten_normieren(Sonnenspektrum Solspec, Fenster);
	// hier müsste überlegt werden, wie man mehrfachkorrekturen umgeht
	void Intensitaeten_durch_piF_Gamma_berechnen(class Speziesfenster Spezfenst,
			double wl_gamma);
	// In der Formel ist piF in W/(m^2*Wellenlänge) verlangt..
	// also muss noch mit der Kanalbreite multipliziert werden
	void Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen(class Speziesfenster Spezfenst);
	void Deklinationswinkel_bestimmen();
	void Sonnen_Longitude_bestimmen();
	void calc_SunEarthDistance();
	double fit_NO_spec(class NO_emiss &NO, std::vector<double> &x,
			std::vector<double> &y, double &rms_err);
	double fit_NO_spec_weighted(class NO_emiss &NO, std::vector<double> &x,
			std::vector<double> &y, std::vector<double> &ye, double &rms_err);
	double fit_rayleigh_and_interp_peaks(class Sonnenspektrum &sol_spec,
			double wl_min, double wl_max, bool debug = true);
	void slant_column_NO(class NO_emiss &NO, std::string mache_Fit_Plots,
			class Sonnenspektrum &sol_spec, int index,
			class Speziesfenster &Spezfenst, std::string Arbeitsverzeichnis,
			class Konfiguration &Konf,
			bool debug = true);

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
	void Ausgabe_in_Datei(std::string Dateiname);

	//Membervariablen

	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	double m_Deklinationswinkel;  // aus Datum berechenbar
	double m_Sonnen_Longitude;
	double m_SunEarthDistance;
	double m_LocalSolarTime;
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
	float m_Sekunde;
	// Geolocation
	double m_Latitude_Sat;   // hierfür muss das Programm noch geändert werden
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Erdradius;

	double m_orbit_phase;

	// limb Membervariablen
	// total number density at measurement point
	double total_number_density;
	// Geolocation
	double m_Latitude_TP;
	double m_Longitude_TP;
	double m_Hoehe_TP;
	double m_TP_SZA;      // alt
	double m_TP_rel_SAA;  // alt - relative Solar azimuth angle
	double center_lat, center_lon;

	// nadir Membervariablen
	// Herkunftsmerkmale
	int m_Messung_ID;
	// Geolokationen für Raytrace
	double m_Latitude_Ground;
	double m_Longitude_Ground;

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
#endif /* MESSUNG_HH_ */
