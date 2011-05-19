/*
 * Messung_Nadir.h
 *
 *  Created on: 28.04.2010
 *      Author: martin
 */

#ifndef MESSUNG_NADIR_HH_
#define MESSUNG_NADIR_HH_

#include <vector>
#include <string>
#include "Ausgewertete_Messung_Nadir.h"

class Messung_Nadir
{
	/****************************************************
	 Der einzige große Unterschied in den Funktionen
	 zwischen Nadir und Limb besteht hier darin,
	 dass ich bei Limb noch Vektors verwendet hab
	 ****************************************************/
public:
	Messung_Nadir();
	// copyconstructor
	Messung_Nadir(const Messung_Nadir &rhs);
	// Assignmentoperator Overload
	Messung_Nadir &operator =(const Messung_Nadir &rhs);
	//Methoden
	int Zeilendichte_Bestimmen(class Speziesfenster &Spezfenst, int Index,
			std::string Arbeitsverzeichnis, std::string mache_Fit_Plots, int MessungsNr);
	//int Sauelendichte_Bestimmen_MgI(Speziesfenster& Spezfenst, int Index,
	//  string Arbeitsverzeichnis, string mache_Fit_Plots, int MessungsNr);
	int Intensitaeten_normieren(std::vector<double> &Sonnen_Intensitaet);
	int Intensitaeten_durch_piF_Gamma_berechnen(class Speziesfenster Spezfenst, int Index);
	int Deklinationswinkel_bestimmen();
	int Sonnen_Longitude_bestimmen();
	Ausgewertete_Messung_Nadir Ergebnis_Zusammenfassen();

	//Hilfsfunktionen
	//nur intern aufrufen!!!
	// Die Hilfsfunktionen sind FAST dieselben wie bei den LIMB Messungen
	int Get_Index(double WL);
	void Fit_Linear(double *x, double *y, double &a0, double &a1,
			int Anfangsindex, int Endindex);
	void Fit_Linear(std::vector<double> &x, std::vector<double> &y, double &a0, double &a1,
			int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(double *x, double *y, double x0, double FWHM,
			double &A, int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(std::vector<double> &x, std::vector<double> &y, double x0, double FWHM,
			double &A, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(double *x, double *y, double a0, double a1,
			double A, double FWHM, double x0, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(std::vector<double> &x, std::vector<double> &y,
			double a0, double a1,
			double A, double FWHM, double x0, int Anfangsindex, int Endindex);

	//Wartungsfunktionen
	//zum Testen und debuggen und überprüfen, ob der fit halbwegs passt
	int Ausgabe_in_Datei(std::string Dateiname);

	//Membervariablen

	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	//Zwischenergebnisse
	double m_Deklinationswinkel;
	double m_Sonnen_Longitude;
	// Herkunftsmerkmale
	std::string m_Dateiname_L1C;
	int m_Messung_ID;
	//Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	float m_Sekunde;
	// Geolokationen für Raytrace
	double m_Latitude_Sat;
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Latitude_Ground;
	double m_Longitude_Ground;
	double m_Erdradius;
	double m_orbit_phase;
	//Füllbare Felder
	int m_Number_of_Wavelength;
	std::vector<double> m_Wellenlaengen;
	std::vector<double> m_Intensitaeten;
	std::vector<double> m_Intensitaeten_relativer_Fehler;
	std::vector<double> m_Intensitaeten_durch_piF;
	std::vector<double> m_Intensitaeten_durch_piF_Gamma;
};

double slit_func(double fwhm, double x0, double x);
double slit_func_gauss(double fwhm, double x0, double x);
#endif /* MESSUNG_NADIR_HH_ */
