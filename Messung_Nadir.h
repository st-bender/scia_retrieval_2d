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
#include"Speziesfenster.h"
#include"Ausgewertete_Messung_Nadir.h"



class Messung_Nadir
{
	/****************************************************
	 Der einzige große Unterschied in den Funktionen
	 zwischen Nadir und Limb besteht hier darin, dass ich bei Limb noch Vektors verwendet hab

	 ****************************************************/
public:
	Messung_Nadir();
	// copyconstructor
	Messung_Nadir(const Messung_Nadir &rhs);
	//Destructor
	~Messung_Nadir();
	// Assignmentoperator Overload
	Messung_Nadir &operator =(const Messung_Nadir &rhs);
	//Methoden
	void save_delete_all_memory();// Speicher löschen
	int Zeilendichte_Bestimmen(Speziesfenster &Spezfenst, int Index, string Arbeitsverzeichnis, string mache_Fit_Plots, int MessungsNr);
	//int Sauelendichte_Bestimmen_MgI(Speziesfenster& Spezfenst, int Index, string Arbeitsverzeichnis, string mache_Fit_Plots, int MessungsNr);
	int Intensitaeten_normieren(double *Teiler);
	int Intensitaeten_durch_piF_Gamma_berechnen(Speziesfenster Spezfenst, int Index);
	int Deklinationswinkel_bestimmen();
	int Sonnen_Longitude_bestimmen();
	Ausgewertete_Messung_Nadir Ergebnis_Zusammenfassen();

	//Hilfsfunktionen  //nur intern aufrufen!!!   // Die Hilfsfunktionen sind FAST dieselben wie bei den LIMB Messungen
	int Get_Index(double WL);
	void Fit_Linear(double *x, double *y, double &a0, double &a1, int Anfangsindex, int Endindex);
	void Fit_Peak_hyperbolic(double *x, double *y, double x0, double FWHM, double &A, int Anfangsindex, int Endindex);
	double Evaluate_Error_primitive(double *x, double *y, double a0, double a1, double A, double FWHM, double x0, int Anfangsindex, int Endindex);

	//Wartungsfunktionen
	int Ausgabe_in_Datei(string Dateiname);  //zum Testen und debuggen und überprüfen, ob der fit halbwegs passt

	//Membervariablen

	// Ergebnisse
	double m_Zeilendichte;
	double m_Fehler_Zeilendichten;
	//Zwischenergebnisse
	double m_Deklinationswinkel;
	double m_Sonnen_Longitude;
	// Herkunftsmerkmale
	string m_Dateiname_L1C;
	int m_Messung_ID;
	//Datum
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	// Geolokationen für Raytrace
	double m_Lattitude_Sat;
	double m_Longitude_Sat;
	double m_Hoehe_Sat;
	double m_Lattitude_Ground;
	double m_Longitude_Ground;
	double m_Erdradius;
	double m_orbit_phase;
	//Füllbare Felder
	int m_Number_of_Wavelength;
	double *m_Wellenlaengen;
	double *m_Intensitaeten;
	double *m_Intensitaeten_relativer_Fehler;
	double *m_Intensitaeten_durch_piF;
	double *m_Intensitaeten_durch_piF_Gamma;
};

#endif /* MESSUNG_NADIR_HH_ */
