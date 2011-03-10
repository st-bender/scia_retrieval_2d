/*
 * Matrizzen_Aufbauen.h
 *
 *  Created on: 10.06.2010
 *      Author: martin
 */

//mit diesen Funktionen werden die Matrizzen für das Retrieval aufgebaut
// Die Funktion Matrizzen_Aufbauen ruft gleich alle diese Funktionen auf einmal auf


#ifndef MATRIZZEN_AUFBAUEN_HH_
#define MATRIZZEN_AUFBAUEN_HH_

#include"MPL_Matrix.h"
#include"MPL_Vektor.h"

#include"Retrievalgitter.h"                          // Luftmassenfaktoren_Matrix_aufbauen
#include"Ausgewertete_Messung_Limb.h"  // Luftmassenfaktoren_Matrix_aufbauen
#include"Ausgewertete_Messung_Nadir.h" // Luftmassenfaktoren_Matrix_aufbauen
#include"Konfiguration.h"                          // Luftmassenfaktoren_Matrix_aufbauen
#include"Speziesfenster.h"                       // Luftmassenfaktoren_Matrix_aufbauen

void Matrizzen_Aufbauen(MPL_Matrix& S_Breite, MPL_Matrix& S_Hoehe,MPL_Matrix& S_letzte_Hoehe, double Lambda_letzte_Hoehe,
                                         MPL_Matrix& S_apriori,MPL_Matrix& S_y, MPL_Matrix& AMF, double Lambda_apriori,
                                         MPL_Matrix Saeulendichten_Fehler, Speziesfenster& Spezies_Fenster,
                                         Retrievalgitter& Grid,
                                         vector<Ausgewertete_Messung_Limb> AM_L,
                                         vector<Ausgewertete_Messung_Nadir> AM_N,
                                         Konfiguration& Konf, int& IERR);
MPL_Matrix Einheitsmatrix_aufbauen(int Dimension);
MPL_Matrix Werte_bei_maximaler_Hoehe_Flagmatrix_Aufbauen(Retrievalgitter Grid);
MPL_Matrix Differenz_von_benachbarten_Zeilenelementen_Matrix_aufbauen(int Zeilen, int Spalten);
MPL_Matrix Differenz_von_benachbarten_Spaltenelementen_Matrix_aufbauen(int Zeilen, int Spalten);
MPL_Matrix Luftmassenfaktoren_Matrix_aufbauen(/*MPL_Matrix& Zeilendichten,*/
                                                                             Retrievalgitter& Grid,
                                                                             vector<Ausgewertete_Messung_Limb> AM_L,
                                                                             vector<Ausgewertete_Messung_Nadir> AM_N,
                                                                             Konfiguration& Konf, Speziesfenster& Spezies_Fenster, int& IERR);
//Hilfsfunktionen

//int interpolieren(int n,const double* x_gegebenes_Gitter,const double* y_gegebenes_Gitter, x_gesuchter_Wert, double& gesuchter_Wert); (funktioniert, wird aber nicht genutzt)
int interpolieren(MPL_Matrix M,int x_Spalte,int y_Spalte, double x_Wert_des_gesuchten_Wertes, double& gesuchter_Wert);

int Pixel_finden_und_AMF_erhoehen_LOS(MPL_Matrix& AMF,Retrievalgitter& Grid,  const int& MessungNr,
                                                                    int& Pixelnummer,
                                                                    const double& Schrittlaenge,const double& Tau_LOS,
                                                                    const double& Punkt_Hoehe, const double& Erdradius,
                                                                    const double& Punkt_Laenge,const double& Punkt_Breite,
                                                                    const double& Phasenfunktion, MPL_Matrix& Tau_LOS_Limb_Matrix);

bool Punkt_Pruefen_und_ggf_AMF_erhoehen(MPL_Matrix& AMF,   Retrievalgitter& Grid , const int& MessungNr,
                                                                        const int& PN_Test                                             , int& Pixelnummer,
                                                                        const double& Schrittlaenge                   , const double& Tau_LOS,
                                                                        const double& Punkt_Hoehe                   , const double& Punkt_Breite,
                                                                        const double& Phasenfunktion, MPL_Matrix& Tau_LOS_Limb_Matrix);
MPL_Vektor Punkt_auf_Strecke_bei_Radius(MPL_Vektor Streckenstartpunkt,
                                                                    MPL_Vektor Streckenvektor, double Radius, double Genauigkeit);
#endif /* MATRIZZEN_AUFBAUEN_HH_ */
