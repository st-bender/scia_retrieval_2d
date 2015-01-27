/*
 * Matrizen_Aufbauen.h
 *
 *  Created on: 10.06.2010
 *      Author: martin
 */

//mit diesen Funktionen werden die Matrizen für das Retrieval aufgebaut
// Die Funktion Matrizen_Aufbauen ruft gleich alle diese Funktionen auf einmal
// auf


#ifndef MATRIZZEN_AUFBAUEN_HH_
#define MATRIZZEN_AUFBAUEN_HH_

#include <vector>
#include "MPL_Matrix.h"
#include "MPL_Vektor.h"

void Matrizen_Aufbauen(MPL_Matrix &S_Breite, MPL_Matrix &S_Hoehe,
						MPL_Matrix &S_letzte_Hoehe, double Lambda_letzte_Hoehe,
						MPL_Matrix &S_apriori, MPL_Matrix &S_y, MPL_Matrix &AMF,
						double Lambda_apriori,
						class Speziesfenster &Spezies_Fenster,
						class Retrievalgitter &Grid,
						std::vector<class Ausgewertete_Messung_Limb> &AM_L,
						std::vector<class Ausgewertete_Messung_Nadir> &AM_N,
						class Konfiguration &Konf, int &IERR);
void generate_Sy(MPL_Matrix &S_y, MPL_Matrix &Saeulendichten_Fehler);
MPL_Matrix Einheitsmatrix_aufbauen(int Dimension);
MPL_Matrix Werte_bei_maximaler_Hoehe_Flagmatrix_Aufbauen(class Retrievalgitter &Grid);
MPL_Matrix Differenz_von_benachbarten_Zeilenelementen_Matrix_aufbauen(int Zeilen, int Spalten);
MPL_Matrix Differenz_von_benachbarten_Spaltenelementen_Matrix_aufbauen(int Zeilen, int Spalten);
MPL_Matrix Luftmassenfaktoren_Matrix_aufbauen(/*MPL_Matrix& Zeilendichten,*/
	class Retrievalgitter &Grid,
	std::vector<class Ausgewertete_Messung_Limb> &AM_L,
	std::vector<class Ausgewertete_Messung_Nadir> &AM_N,
	class Konfiguration &Konf, class Speziesfenster &Spezies_Fenster, int &IERR);
//Hilfsfunktionen

//int interpolieren(int n,const double* x_gegebenes_Gitter,
//const double* y_gegebenes_Gitter, x_gesuchter_Wert, double& gesuchter_Wert);
//(funktioniert, wird aber nicht genutzt)
int interpolieren(MPL_Matrix &M, int x_Spalte, int y_Spalte,
		double x_Wert_des_gesuchten_Wertes, double &gesuchter_Wert);

int Pixel_finden_und_AMF_erhoehen_LOS(MPL_Matrix &AMF, class Retrievalgitter &Grid,
									  const int &MessungNr,
									  int &Pixelnummer,
									  const double &Schrittlaenge,
									  const double &Tau_LOS,
									  const double &Punkt_Hoehe,
									  const double &Erdradius,
									  const double &Punkt_Laenge,
									  const double &Punkt_Breite,
									  const double &Phasenfunktion,
									  MPL_Matrix &Tau_LOS_Limb_Matrix);

bool Punkt_Pruefen_und_ggf_AMF_erhoehen(MPL_Matrix &AMF, class Retrievalgitter &Grid,
										const int &MessungNr,
										const int &PN_Test, int &Pixelnummer,
										const double &Schrittlaenge,
										const double &Tau_LOS,
										const double &Punkt_Hoehe,
										const double &Punkt_Breite,
										const double &Phasenfunktion,
										MPL_Matrix &Tau_LOS_Limb_Matrix);
MPL_Vektor Punkt_auf_Strecke_bei_Radius(MPL_Vektor &Streckenstartpunkt,
										MPL_Vektor &Streckenvektor,
										double Radius, double Genauigkeit);

void prepare_total_density(class Retrievalgitter &grid, MPL_Matrix &dens,
		std::vector<class Ausgewertete_Messung_Limb> &aml_vec);
void SNOE_apriori_NO(class Retrievalgitter &grid,
		class Ausgewertete_Messung_Limb &aml, MPL_Matrix &apriori,
		class Konfiguration &Konf);
void regression_apriori_NO(class Retrievalgitter &grid,
		class Ausgewertete_Messung_Limb &aml, MPL_Matrix &apriori,
		class Konfiguration &Konf);
void scale_apriori(class Retrievalgitter &grid, MPL_Matrix &apriori,
		class Konfiguration &Konf);
#endif /* MATRIZZEN_AUFBAUEN_HH_ */
