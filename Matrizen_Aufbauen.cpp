/*
 * Matrizen_Aufbauen.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 10.06.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */


#include"Matrizen_Aufbauen.h"
#include"MPL_Matrix.h"
#include"MPL_Vektor.h"

#include"Retrievalgitter.h"            // Luftmassenfaktoren_Matrix_aufbauen
#include"Ausgewertete_Messung_Limb.h"  // Luftmassenfaktoren_Matrix_aufbauen
#include"Ausgewertete_Messung_Nadir.h" // Luftmassenfaktoren_Matrix_aufbauen
#include"Winkelstatistik.h"
#include"Speziesfenster.h"             // Luftmassenfaktoren_Matrix_aufbauen

// Luftmassenfaktoren_Matrix_aufbauen(einlesen von Atmosphärendaten)
#include"Datei_IO.h"
#include"Konfiguration.h"             // Dateinamen für die Atmosphärendaten
#include <cmath>                      // trigonometrische Funktionen

#include "Koordinatentransformation.h"
#include "Glaetten.h"
#include "NO_regression_model.h"

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iterator>
#include <numeric>
#include <functional>


using std::cerr;
using std::cout;
using std::endl;
using std::vector;

void Matrizen_Aufbauen(MPL_Matrix &S_Breite, MPL_Matrix &S_Hoehe,
						MPL_Matrix &S_letzte_Hoehe, double Lambda_letzte_Hoehe,
						MPL_Matrix &S_apriori, MPL_Matrix &S_y, MPL_Matrix &AMF,
						double Lambda_apriori,
						Speziesfenster &Spezies_Fenster,
						Retrievalgitter &Grid,
						vector<Ausgewertete_Messung_Limb> &AM_L,
						vector<Ausgewertete_Messung_Nadir> &AM_N,
						Konfiguration &Konf, int &IERR)
{
	//Fehlermatrizen
	//cerr<<"S_Breite\n";
	S_Breite = Differenz_von_benachbarten_Zeilenelementen_Matrix_aufbauen(
				Grid.m_Anzahl_Hoehen, Grid.m_Anzahl_Breiten);
	//cerr<<"S_Hoehe\n";
	S_Hoehe = Differenz_von_benachbarten_Spaltenelementen_Matrix_aufbauen(
				Grid.m_Anzahl_Hoehen, Grid.m_Anzahl_Breiten);
	//cerr<<"S_letzte_Hoehe\n";
	S_letzte_Hoehe = Werte_bei_maximaler_Hoehe_Flagmatrix_Aufbauen(Grid);
	//cerr<<"S_letzte_Hoehe_erhoehen\n";
	S_letzte_Hoehe *= Lambda_letzte_Hoehe;

	//cerr<<"S_apriori\n";
	S_apriori = Einheitsmatrix_aufbauen(Grid.m_Anzahl_Punkte);

	// Fehlermatrizen mit konstanten Wichtungsfaktoren Multiplizieren
	// nur bei apriori...die anderen erst in Normalgleichung
	S_apriori *= Lambda_apriori;
	// Flagmatrix oberste Hoehe

	//cerr<<"AMF\n";
	//Raytracing zum Aufbau der Luftmassenfaktoren Matrix durchführen
	//(sehr große Funktion)
	//Die Funktion wurde überprüft beim debuggen 28.9.2010;
	//Ende 4.10.2010...sieht erstmal fehlerfrei aus
	AMF = Luftmassenfaktoren_Matrix_aufbauen(Grid, AM_L, AM_N, Konf,
			Spezies_Fenster, IERR);
	//cerr<<"nach AMF\n";
}
void generate_Sy(MPL_Matrix &S_y, MPL_Matrix &Saeulendichten_Fehler)
{
	// Wichtungsfaktorenmatrix der Messwerte 1/Fehler^2 von y auf Diagonalen
	//cerr<<"S_y\n";
	S_y = Einheitsmatrix_aufbauen(Saeulendichten_Fehler.m_Elementanzahl);
	for (int i = 0; i < Saeulendichten_Fehler.m_Elementanzahl; i++) {
		double d = Saeulendichten_Fehler(i);
		S_y(i, i) = 1 / (d * d);
	}
}

//==============================================================================
//
//==============================================================================
////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Einheitsmatrix_aufbauen()
////////////////////////////////////////////////////////////////////////////////
MPL_Matrix Einheitsmatrix_aufbauen(int Dimension)
{
	//Einheitsmatrix muss quadratisch sein
	MPL_Matrix Eins(Dimension, Dimension);
	// Diagonalelemente 1 setzen

	for (int i = 0; i < Dimension; i++) {
		Eins(i, i) = 1;
	}
	return Eins;
} //ende Einheitsmatrix_aufbauen(int Dimension)
////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Diagonalmatrix_aufbauen()
////////////////////////////////////////////////////////////////////////////////
MPL_Matrix Diagonalmatrix_aufbauen(MPL_Matrix &V)
{
	// Copy all elements of V (which is thereby flattened).
	// Ideally, V is a row or column vector.
	MPL_Matrix Diag{V.m_Elementanzahl, V.m_Elementanzahl};

	for (int i = 0; i < V.m_Elementanzahl; ++i) {
		Diag(i, i) = V(i);
	}
	return Diag;
} //ende Diagonalmatrix_aufbauen
//==============================================================================
//
//==============================================================================
////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Werte_bei_maximaler_Hoehe_Flagmatrix_Aufbauen
////////////////////////////////////////////////////////////////////////////////
MPL_Matrix Werte_bei_maximaler_Hoehe_Flagmatrix_Aufbauen(Retrievalgitter &Grid)
{
	MPL_Matrix Flagmatrix(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// Auf den Diagonalen dort eine 1 schreiben, wo die höchste Höhe liegt....
	// das sollten die letzten Punkte der Matrix sein
	//cerr<<"Grid.m_Anzahl_Punkte: "<<Grid.m_Anzahl_Punkte<<"\n";
	for (int i = 0; i < Grid.m_Anzahl_Breiten; i++) {
		//cerr<<"i: "<<i<<"\n";
		//cerr<<"Grid.m_Anzahl_Punkte-1-i: "<<Grid.m_Anzahl_Punkte-1-i<<"\n";
		Flagmatrix(Grid.m_Anzahl_Punkte - 1 - i, Grid.m_Anzahl_Punkte - 1 - i)
			= 1;
	}
	//cerr<<"fertig\n";
	return Flagmatrix;
}
////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Werte_bei_maximaler_Hoehe_Flagmatrix_Aufbauen
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Differenz_von_benachbarten_Zeilenelementen_Matrix_aufbauen()
////////////////////////////////////////////////////////////////////////////////
MPL_Matrix Differenz_von_benachbarten_Zeilenelementen_Matrix_aufbauen(
		int Zeilenzahl, int Spaltenzahl)
{
	// mit Zeilenzahl und Spaltenzahl sind die ursprünglichen Dimensionen der
	// Matrix gemeint, die Zeilenweise in einen Vektor gepresst wurde...Vektor
	// und Matrix haben Zeilenzahl*Spaltenzahl viele Elemente...um nun
	// benachbarte Zeilenelemente der Matrix zu erkennen, müssen die
	// urdimensionen angegeben werden

	// Input war die ursprüngliche Zeilen und Spaltenzahl, aus der der Vektor
	// der Länge Elementzahl aufgebaut ist Matrix ist quadratisch
	int Elementzahl = Zeilenzahl * Spaltenzahl;
	// Quadratische Matrix mit der 2Dimensionen der Länge des Vektors
	MPL_Matrix Mat(Elementzahl, Elementzahl);

	// Matrix wird Zeilenweise in einen Vektor gepackt, d.h. aus x1 x2 x3 x4
	//                                                           x5 x6 x7 x8
	//                                      wird x1 x2 x3 x4 x5 x6 x7 x8 gemacht
	// Zeilennachbarn sind also in den meisten Fällen direkt benachbarte
	// Elemente Für das ite Element des neuen Vektors wird immer das i-te
	// Element des alten Vektors vom i+1 Elementer der Zeile abgezogen. Das
	// letzte Element der Zeile, was übrig bleibt muss 0 sein(Später wird eh
	// das Quadrat des Vektors als Skalarprodukt gebildet, da soll dann in dem
	// Element auch +0 rauskommen ).

	for (int i = 0; i < Elementzahl; i++) {
		Mat(i, i) = -1;
		//Prüfen, ob x nicht das Maximum der ursprünglichen Zeile war
		if (((i + 1) % Spaltenzahl) != 0) {
			Mat(i, i + 1) = 1;
		}
	}
	return Mat;
}// Ende Differenz_von_benachbarten_Zeilenelementen_Matrix_aufbauen
//==============================================================================
//
//==============================================================================
////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Differenz_von_benachbarten_Spaltenelementen_Matrix_aufbauen()
////////////////////////////////////////////////////////////////////////////////
MPL_Matrix Differenz_von_benachbarten_Spaltenelementen_Matrix_aufbauen(
		int Zeilenzahl, int Spaltenzahl)
{
	// Input war die ursprüngliche Zeilen und Spaltenzahl, aus der der Vektor
	// der Länge Elementzahl aufgebaut ist
	int Elementzahl = Zeilenzahl * Spaltenzahl;
	// Quadratische Matrix mit der 2Dimensionen der Länge des Vektors
	MPL_Matrix Mat(Elementzahl, Elementzahl);
	// Matrix wird Zeilenweise in einen Vektor gepackt, d.h. aus x1 x2 x3 x4
	//                                                           x5 x6 x7 x8
	//                                      wird x1 x2 x3 x4 x5 x6 x7 x8 gemacht
	// Spaltennachbarn sind immer um eine Spaltenbreite versetzt
	// Wieder wird für das i-te Element der Spalte das i-te Elemente vom i+1ten
	// Element der Spalte abgezogen die Elemente, die 0 sein sollen sind in
	// dieser Reihenfolge am Ende der Matrix zu finden, sodass man sich die if
	// Abfrage sparen kann
	for (int i = 0; i < Elementzahl; i++) {
		Mat(i, i) = -1;
		if (i < Elementzahl - Spaltenzahl)
			Mat(i, i + Spaltenzahl) = 1;
	}
	return Mat;
}// Ende  Differenz_von_benachbarten_Spaltenelementen_Matrix_aufbauen
//==============================================================================
//
//==============================================================================

////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Luftmassenfaktoren_Matrix_aufbauen()
////////////////////////////////////////////////////////////////////////////////
MPL_Matrix Luftmassenfaktoren_Matrix_aufbauen(/*MPL_Matrix& Zeilendichten,*/
	Retrievalgitter &Grid,
	vector<Ausgewertete_Messung_Limb> &AM_L,
	vector<Ausgewertete_Messung_Nadir> &AM_N,
	Konfiguration &Konf, Speziesfenster &Spezies_Fenster, int &IERR)
{
	IERR = 0; //0 ist ok
	//const double pi=3.1415926535897;
	// Die Luftmassenfaktoren(Die Matrixelemente) haben die EINHEIT cm
	//Hier werden nun die Luftmassenfaktoren unseres Vorwärtsmodells bestimmt,
	//dabei steht jede Zeile der Matrix für eine unabhängige Gleichung zur
	//Errechnung einer Zeilendichte Jede Zeile hat also soviele Elemente wie
	//das Retrievalgitter Punkte hat
	//Die Anzahl der Zeilen entspricht der Anzahl der Messungen
	//Luftmassenfaktormatrix
	MPL_Matrix AMF(AM_L.size() + AM_N.size(), Grid.m_Anzahl_Punkte);

	////////////////////////////////////////////////////////////////////////////
	//Zum Test wird noch die Tau_LOS_Limb Matrix eingeführt, die später mal
	//geplottet werden kann Diese beinhaltet den letzten Eintrag von Tau in dem
	//Gitterelement(nicht das durchschnittliche)
	//Zum plot mal die Gittergrenzen in den Breitengraden fein machen
	//An Tau sollte dann der Strahlengang besser sichtbar werden
	// Kann später entfernt werden...(ist aber auch nicht langsam)
	//cerr<<"TAU.Null_Initialisierung\n";
	MPL_Matrix Tau_LOS_Matrix(AM_L.size() + AM_N.size(), Grid.m_Anzahl_Punkte);
	////////////////////////////////////////////////////////////////////////////

	//Zunächst müssen die Athmosphärendaten eingelesen werden
	// Es werden also die Dichten der Atmosphärengase Luft und Ozon bei
	// verschiedenen Höhen eingelesen und deren Wirkungsquerschnitte bei
	// verschiedenen Wellenlängen
	// Dann wird für die jeweilige Wellenlänge des Metallemissionsübergangs der
	// jeweilige Wirkungsquerschnitt herausgesucht /Damit die Suche schneller
	// geht, werden die wenigen interessanten Wellenlänge und
	// Wirkungsquerschnittspaare in einen kleineren Vektor geschrieben.
	////// Atmosphärendaten einlesen////////////////////////////////////////////
	//cerr<<"Atmosphärendaten einlesen\n";
	MPL_Matrix M_Atmo_Dichten;
	M_Atmo_Dichten = Read_Atmodatei(Konf.m_Pfad_Dichten_der_Atmosphaerengase);
	MPL_Matrix M_Atmo_Wirkungsquerschnitte;
	M_Atmo_Wirkungsquerschnitte =
		Read_Atmodatei(Konf.m_Pfad_Wirkungsquerschnitte_der_Atmosphaerengase);
	// hier kann man noch Geschwindigkeit rausholen, wenn man nur die nur die 2
	// oder etwas mehr benötigten Querschnitte verwendet ( Funktion bremst aber
	// nicht)
	//////Ende Atmosphärendaten einlesen////////////////////////////////////////
	//Für das Raytracing müssen nun einige Skizzen bertrachtet werden, die ich
	//in eine PDF-Datei-Plotten werde Irgendwann mal machen ;) : PDF-Datei mit
	//Skizzen erzeugen
	//cerr<<"Limb Raytracing\n";
	////////////////////////////////////////////////////////////////////////////
	// Zunächst die Limbmessungen in Eintragen  ////////////////////////////////
	Winkelstatistik Wstat;  //kleine Winkelstatistik initialisieren
	//TODO früher im Programm anbringen

	vector<Ausgewertete_Messung_Limb>::iterator aml_it;

	for (aml_it = AM_L.begin(); aml_it != AM_L.end(); ++aml_it) {
		unsigned int MessungNR = distance(AM_L.begin(), aml_it);

		// Für jede Teilchensorte können mehrere Linien ausgewertet werden.
		// Welche gerade verwendet wird, wird ermittelt, Aus der Wellenlänge
		// des Übergangs sind die (beiden) Wirkungsquerschnitte für die beiden
		// Atmosphärengase zu bestimmen
		MPL_Vektor V_Atmo_Wirkungsquerschnitte(2);
		MPL_Vektor V_Atmo_Wqs_abs(2);
		//Dichten sind Höhenabhängig, Bestimmung dort
		MPL_Vektor V_Atmo_Dichten(2);

		// Wellenlängen sind in Atmo-Datei monoton ansteigend-> Quicksearch
		// Algorithmus für Interpolation(schnell)
		//int interpolieren(MPL_Matrix M,int x_Spalte,int y_Spalte,
		//   double //x_Wert_des_gesuchten_Wertes, double& gesuchter_Wert);
		interpolieren(M_Atmo_Wirkungsquerschnitte, 0, 1,
				aml_it->m_Wellenlaenge, V_Atmo_Wirkungsquerschnitte(0));
		interpolieren(M_Atmo_Wirkungsquerschnitte, 0, 2,
				aml_it->m_Wellenlaenge, V_Atmo_Wirkungsquerschnitte(1));
		// absorption cross section, might be different from the emission
		interpolieren(M_Atmo_Wirkungsquerschnitte, 0, 1,
				aml_it->m_Wellenlaenge_abs, V_Atmo_Wqs_abs(0));
		interpolieren(M_Atmo_Wirkungsquerschnitte, 0, 2,
				aml_it->m_Wellenlaenge_abs, V_Atmo_Wqs_abs(1));
		//std::cerr << "# MessungNR = " << MessungNR << ", wl = "
		//	<< aml_it->m_Wellenlaenge << ", wl_abs = "
		//	<< aml_it->m_Wellenlaenge_abs << std::endl;
		//std::cerr << "# WQs: emiss " << V_Atmo_Wirkungsquerschnitte(0)
		//	<< ", " << V_Atmo_Wirkungsquerschnitte(1) << "; abs: "
		//	<< V_Atmo_Wqs_abs(0) << ", " << V_Atmo_Wqs_abs(1) << std::endl;
		//zwei Wege müssen betrachtet werden:
		//Satellit-Punkt
		//Punkt-Sonne
		//und zwar in der Richtung und der Reihenfolge,
		//falls man aufaddieren will
		////////////////////////
		//Satellit-Punkt //
		////////////////////////
		////////////////////////////////////////////////////////////////////////
		// LOS Raytracing für Limbmessungen
		////////////////////////////////////////////////////////////////////////
		// Funktion
		//   Umwandlung_Kugel_in_Karthesisch(double r,double phi,double theta,
		//     double& x,double& y,double& z);
		// Funktion
		//   Umwandlung_Karthesisch_in_Kugel(double x,double y,double z,
		//     double& r,double& phi,double& theta);
		// Position des Satelliten in kartesischen Koordinaten
		MPL_Vektor Sat_Pos(3);
		Umwandlung_Kugel_in_Karthesisch(aml_it->m_Erdradius
											+ aml_it->m_Hoehe_Sat,
										aml_it->m_Longitude_Sat,
										aml_it->m_Latitude_Sat,
										Sat_Pos(0), Sat_Pos(1), Sat_Pos(2));
		// Position des Tangentenpunkts in kartesischen Koordinaten
		MPL_Vektor TP_Pos(3);
		Umwandlung_Kugel_in_Karthesisch(aml_it->m_Erdradius
											+ aml_it->m_Hoehe_TP,
										aml_it->m_Longitude_TP,
										aml_it->m_Latitude_TP,
										TP_Pos(0), TP_Pos(1), TP_Pos(2));
		// Verbindungsvektor Sat-TP startend vom Satelliten
		MPL_Vektor Verbindungsvektor(TP_Pos - Sat_Pos);
		MPL_Vektor Verbindungsvektor_normal(Verbindungsvektor);
		Verbindungsvektor_normal.Normieren();
		//START kleine Statistik für Winkelabweichungen bei Limb   /////////////
		int Winkel_OK = 3;
		Wstat.Winkel_berechnen_und_einordnen(Verbindungsvektor,
				TP_Pos, Winkel_OK);
		if (Winkel_OK != 0) {
			cout << "Winkel in Messung " << MessungNR
				 << "zu groß... überspringe Messung\n";
			continue;
			//Führt zu leerer Zeile in Matrix....das ist nicht gut....
			//deshalb vor dem Matrixaufbau rausschmeißen
		}
		//TODO Messungen mit zu großen Abweichungen rausschmeißen...
		// nicht hier...früher
		// ENDE kleine Statistik für Winkelabweichungen bei Limb   /////////////

		//Die Höchste Höhe steckt in der letzten Zeile, also im Zweifel im
		//letzten Element Achtung hier muss man immer aufpassen:es gibt 2 Höhen
		//die vom Erdkern und die von der Erdoberfläche
		double Hoehe_TOA = aml_it->m_Erdradius
			+ Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1].m_Max_Hoehe;
		//cerr<<"Hoehe_TOA: "<<Hoehe_TOA<<"\n";
		//cout<<"Maxhoehe: "<<Grid.m_Gitter[Grid.m_Anzahl_Punkte-1].m_Max_Hoehe<<"\n";
//        if (Hoehe_TOA>(200.0+aml_it->m_Erdradius))
//        {    Hoehe_TOA=200.0+aml_it->m_Erdradius;}
		// GENAUIGKEIT prüfen, ob 200 nicht sinnvoller ist...bzw das aus der
		// Datei geschrieben wird...eventuell exponentialfunktion testen
		double Max_Hoehe_Absorption = Konf.m_TOA + aml_it->m_Erdradius;
		// das ist was anderes als die maxhoehe...
		// sonst leere(singuläre) Spalten in Matrix
		//  genauer ist bis 200
		//cerr<<"Max_Hoehe_Absorption: "<<Max_Hoehe_Absorption<<"\n";
		// cout<<"Max_Hoehe_Absorption:"<<Max_Hoehe_Absorption<<"\n";
		if (TP_Pos.Betrag_ausgeben() > Hoehe_TOA) {
			// prevents inifite loops below
			cout << "TP außerhalb des Grids, Anzahl_Hoehen ist zu erhöhen."
				 << endl;
			continue;
		}
		double Sehnenlaenge_in_Atmosphaere = 2. * sqrt(Hoehe_TOA*Hoehe_TOA
				- std::pow(aml_it->m_Erdradius + aml_it->m_Hoehe_TP, 2));
		//cerr<<"Sehnenlaenge_in_Atmosphaere: "<<Sehnenlaenge_in_Atmosphaere<<"\n";
		//cout<<"Weglaenge in Atmosphäre: "<<Sehnenlaenge_in_Atmosphaere<<" km.\n";
		////////////////////////////////////////////////////////////////////////
		// GESCHWINDIGKEITS/GENAUIGKEITSBESTIMMENDER PARAMETER
		int Schrittzahl = 10000;
		// 40000 ist etwa alle 110m     ; 10000 ist wie in Marcos Programm
		////////////////////////////////////////////////////////////////////////
		//cout<<Schrittzahl<<" Schritte der Länge"
		//  <<Sehnenlaenge_in_Atmosphaere/Schrittzahl<<".\n";
		int Winkelberechnungsfrequenz = 100; //d.h. alle x Schritte
		//cout<<"Winkelberechnung alle "<< Winkelberechnungsfrequenz
		//  <<" Schritte, also "<<Schrittzahl/Winkelberechnungsfrequenz<<"mal.\n";
		MPL_Vektor TOA_Start(3);
		TOA_Start = TP_Pos
			- 0.5 * (Verbindungsvektor_normal * Sehnenlaenge_in_Atmosphaere);
		MPL_Vektor TOA_Schritt(3);
		//   Einheitsvektor in Richtung      Gesamtlänge     Anzahl Teilstücke
		TOA_Schritt = (Verbindungsvektor_normal * Sehnenlaenge_in_Atmosphaere)
			/ (double)Schrittzahl;
		//Sehnenlänge ist auf 100m genau...deshalb später überschüssige Punkte
		//abfangen

		//cout<<"Verbindungsvektor: "<<Verbindungsvektor(0)<<"\t"
		//  <<Verbindungsvektor(1)<<"\t"<<Verbindungsvektor(2)<<"\n";
		//cout<<"Betrag_Verbindungsvektor: "<<Verbindungsvektor.Betrag_ausgeben()<<"\n";
		//cout<<"Sehnenlänge in atmo: "<<Sehnenlaenge_in_Atmosphaere<<"\n";
		//cout<<"schrittzahl: "<<Schrittzahl<<"\n";
		//cout<<"TOA_Schritt.Betrag: "<<TOA_Schritt.Betrag_ausgeben()<<"\n";
		//sleep(2);
		//cout<<"TOA_Schritt: "<<TOA_Schritt(0)<<"\t"
		//  <<TOA_Schritt(1)<<"\t"<<TOA_Schritt(2)<<"\n";
		//initialisierung...keine Absorption bei TOA auf der Satellitenseite
		double Tau_LOS = 0.;
		// Startpixel festlegen... im alten Algorithmus werden alle ca. 200
		// Pixel untersucht, ob der neu berechnete Punkt drin liegt.  Das ist
		// nicht nötig. Da es reicht, wenn man den alten Punkt kennt, nur
		// dessen 8 Nachbarn zu untersuchen,
		// Beschleunigung  200/8rund=20 mal so schnell  bei 10001 Punkten für
		// jede Messung ist das vermutlcih sogar Zeitkritisch
		int Pixelnummer = -1; // nicht alle Pixel müssen definiert sein
		//IST DAS WICHTIG? liegen alle Punkte,die in keinem Pixel liegen so
		//hoch, das Absorption vernachlässigbar?  AN den POLEN macht das
		//Probleme ..deshalb besser gitter mit Orbitphase
		// Vor dem Raytracing, sollten auch alle Durchstoßpunkte Null gesetzt
		// werden
		Grid.Alle_Durchstosspunkte_Null_setzen();

		// RAYTRACINGSCHLEIFE LOS///////////////////////////////////////////////
		double Cos_Streuwinkel = -1.;
		for (int aktuelle_Schritt_Nr = 0;
				aktuelle_Schritt_Nr < Schrittzahl; aktuelle_Schritt_Nr++) {
			//cerr<<"aktuelle_Schritt_Nr: "<<aktuelle_Schritt_Nr<<"\n";
			//cout<<"aktuelle Schrittnummer: "<<aktuelle_Schritt_Nr<<"\n";
			// Schritt machen
			MPL_Vektor aktueller_Punkt(3);
			aktueller_Punkt = TOA_Start
				+ (double)aktuelle_Schritt_Nr * TOA_Schritt;

			double Punkt_Hoehe;
			double Punkt_Radius;
			double Punkt_Laenge;
			double Punkt_Breite;
			// Längen und Breitengrad, Höhe und Radius berechnen
			Umwandlung_Karthesisch_in_Kugel(aktueller_Punkt(0),
					aktueller_Punkt(1), aktueller_Punkt(2),
					Punkt_Radius, Punkt_Laenge, Punkt_Breite);
			//cerr<<"Punkt_Radius: "<<Punkt_Radius<<"\n";
			if (Punkt_Radius > Hoehe_TOA) {
				//cout<<"Hoehe_TOA:"<<Hoehe_TOA<<"\n";
				//cout<<"Punkt_Radius:"<<Punkt_Radius<<"\n";
				continue;
			}
			// Am Ende gibts aufgrund der Toleranz in der Sehnenlängenbestimmung
			// mal ein par Punkte mehr, die nicht beitragen
			Punkt_Hoehe = Punkt_Radius - aml_it->m_Erdradius;
			//cerr<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";

			//Testen ob Streuwinkel berechnet werden sollen
			if (aktuelle_Schritt_Nr % Winkelberechnungsfrequenz == 0) {

				//cout<<"aktuelle Schrittnummer: "<<aktuelle_Schritt_Nr<<"\n";
				//sleep(1);
				// Streuwinkel berechnen...hier nur Streuwinkel interessant
				// nicht SZA Benötigt wird die Sonnenposition....r wird nicht
				// gemessen, ist aber für die relevanten Messungen unendlich
				// groß also nur eine ungefähre Distanzangabe nötig
				// Position des Tangentenpunkts in kartesischen Koordinaten
				MPL_Vektor Sonne_Pos(3);
				Umwandlung_Kugel_in_Karthesisch(149.6E6,
												aml_it->m_Sonnen_Longitude,
												aml_it->m_Deklination,
												Sonne_Pos(0), Sonne_Pos(1), Sonne_Pos(2));
				// Das normerte Skalarprodukt der Sonnenposition und des Line
				// of Sight-Vektors liefert gerade den cosinus des Winkels, der
				// mit dem Streuwinkel zusammen 180° ergibt das entpricht
				// gerade einer Spiegelung am Einheitskreis an der y achse
				// cos(A)=cos(180-B)=cacb+sasb=-1cb=-cos(B)
				// da nur das Quadrat benötigt wird, ist cos(A)^2=cos(B)^2
				// da nicht der Winkel sondern eh nur der Cosinus gesucht wird,
				// gibts auch keine Verwirrungen, ob der Stumpfe oder der
				// Spitze Winkel gemeint ist
				//
				MPL_Vektor Verbindung(Sonne_Pos - aktueller_Punkt);
				MPL_Vektor LOS_normal(Verbindungsvektor_normal);
				Verbindung.Normieren();
				Cos_Streuwinkel = -1 * (Verbindung * LOS_normal); //skalarprodukt
				//cerr<<"Streuwinkel berechnet\n";
			}

			//TAU+=n sigma erhöhen
			// Wirkungsquerschnitte schon interpoliert....
			// nun n aus n-höhen information interpolieren
			//cout<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
			V_Atmo_Dichten.Null_Initialisierung();
			//time_t Zeit1,Zeit2,deltaZeit;
			//time(&Zeit1);
			if (Punkt_Radius <= Max_Hoehe_Absorption) {
				//cout<<"Max_Hoehe_Absorption:"<<Max_Hoehe_Absorption<<"\n";
				//cout<<"Punkt_Radius:"<<Punkt_Radius<<"\n";
				//verhindert dass interpolation fehler ausgibt...
				//über 200km eh keine Absorption mehr
				// SPEEDUP!!!!
				interpolieren(M_Atmo_Dichten, 0, 1, Punkt_Hoehe, V_Atmo_Dichten(0));
				interpolieren(M_Atmo_Dichten, 0, 2, Punkt_Hoehe, V_Atmo_Dichten(1));
				// SPEED Die Funktion hier braucht fast 10 Sekunden...
				// und das 2 mal
				//Eventuell interpolierte Datei erzeugen und den nächsten Punkt
				//laden oder lookup table erzeugen damit (oder Schrittzalh
				//erniedrigen

				// Datei sieht ok aus
				//M_Atmo_Dichten.in_Datei_speichern("/home/martin/TEMP/Atmodichten.txt");
			} //falls diese Bedingung nicht erfüllt sind,
			  // sind die dichten null und tau wird um 0 erhöht
			//Woher kommt die 100000-> Schrittlänge wird in km angegeben die
			//Dichten in 1/cm^3 und die WirkungsQS in cm^2 Tau ist Einheitenlos
			////Skalarprodukt beider Vektoren
			double Schrittlaenge
				= Sehnenlaenge_in_Atmosphaere / (double)Schrittzahl;
			//cerr<<"Schrittlaenge: "<<Schrittlaenge<<"\n";
			Tau_LOS += Schrittlaenge * 100000.0
					  * (V_Atmo_Dichten * V_Atmo_Wirkungsquerschnitte);
			//cout<<"V_Atmodichten(0)"<<V_Atmo_Dichten(0)<<"\n";
			//cout<<"V_Atmodichten(1)"<<V_Atmo_Dichten(1)<<"\n";
			//cout<<"V_Wirkungsquerschnitte(0)"<<V_Atmo_Wirkungsquerschnitte(0)<<"\n";
			//cout<<"V_Wirkungsquerschnitte(1)"<<V_Atmo_Wirkungsquerschnitte(1)<<"\n";
			//cout<<"Tau_LOS: "<<Tau_LOS<<"\n";
			// Überprüfen, ob Punkt in Dunkelheit liegt
			// (Nord ist +pi/2 süd -pi/2)
			// (90.0 ist standardmäßig double literal)
			//cerr<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
			//cerr<<"aml_it->m_Deklination: "<<aml_it->m_Deklination<<"\n";
			if ((Punkt_Breite > 90.0 + aml_it->m_Deklination)
					|| (Punkt_Breite < -90.0 + aml_it->m_Deklination)) {
				if (aktuelle_Schritt_Nr == Schrittzahl - 1) {
					//cout<<"letzter Punkt und im dunkeln\n";
					if (Pixelnummer != -1 &&
						Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt.Betrag_ausgeben() == 0) {
							cout << "MessungNR: " << MessungNR << "\n";
							cout << "fehlenden hinteren Durchstoßpunkt setzen\n";
							Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + aml_it->m_Erdradius,
									Punkt_Laenge,
									Punkt_Breite,
									Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(0),
									Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(1),
									Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(2));
					}
				}
				continue;
			}
			//Falls die alte Pixelnummer bekannt ist,
			//brauchen wir nur die Nachbarn absuchen
			// Funktion int Find_Pixel_and_increase_AMF
			//Parameter E1 und E2 für Phasenfunktion
			//cerr<<"bla\n";
			double E1 = 0., E2 = 1.;
			unsigned int PeakNr;
			for (PeakNr = 0; PeakNr < Spezies_Fenster.m_Wellenlaengen.size(); PeakNr++) {
				//cerr<<"Spezies_Fenster.m_Wellenlaengen.size(): "
				//  <<Spezies_Fenster.m_Wellenlaengen.size()<<"\n";
				//cerr<<"PeakNr: "<<PeakNr<<"\n";
				//cerr<<"aml_it->m_Wellenlaenge: "
				//  <<aml_it->m_Wellenlaenge<<"\n";
				//cerr<<"Spezies_Fenster.m_Wellenlaengen[PeakNr]: "
				//  <<Spezies_Fenster.m_Wellenlaengen[PeakNr]<<"\n";
				if (aml_it->m_Wellenlaenge == Spezies_Fenster.m_Wellenlaengen[PeakNr]) {
					E1 = Spezies_Fenster.m_Liniendaten[PeakNr].m_E1;
					E2 = Spezies_Fenster.m_Liniendaten[PeakNr].m_E2;
					//cerr<<"E1: "<<E1<<"\n";
					//cerr<<"E2: "<<E2<<"\n";
					break;
				}
			}
			if (PeakNr == Spezies_Fenster.m_Wellenlaengen.size()) {
				//eigentlich kann das gar nicht passieren
				//cout << "E1 und E2 konnten nicht gefunden werden..."
				//	 << "mysteriöser bug\n";
				E1 = 0.0;
				E2 = 1.0;
			}
			//Phasenfunktion aus Streuwinkel bestimmen
			//und mit AMF multiplizieren...in Funktion
			//cerr<<"vor Phasenfunktion\n";
			double Phasenfunktion = 0.75 * E1
				* (Cos_Streuwinkel * Cos_Streuwinkel + 1) + E2;
			//cerr<<"Phasenfunktion: "<<Phasenfunktion<<"\n";
			//cout<<"Phasenfunktion: "<<Phasenfunktion<<"\n";
			//cout<<"E1: "<<E1<<"\n";
			//cout<<"E2: "<<E2<<"\n";
			//cout<<"Cos_Streuwinkel: "<<Cos_Streuwinkel<<"\n";
			//cout<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
			//cout<<"Punkt_Laenge: "<<Punkt_Laenge<<"\n";
			//cout<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
			//cout<<"Tau_LOS: "<<Tau_LOS<<"\n";

			int myerr = Pixel_finden_und_AMF_erhoehen_LOS(AMF, Grid, MessungNR,
					Pixelnummer,
					Schrittlaenge, Tau_LOS,
					Punkt_Hoehe, aml_it->m_Erdradius,
					Punkt_Laenge, Punkt_Breite,
					Phasenfunktion, Tau_LOS_Matrix);
			if (myerr == 4) {
				IERR = 1;
				cout << "mysteriöser fall beendet Programm\n";
				MPL_Matrix dummy(1, 1);
				return dummy;
			}

			//letzten Punkt als hinteren Durchstoßpunkt nutzen, des letzen
			//Feldes nutzen, falls zufällig erster und letzter Punkt gleich,
			//wird das nicht gemacht
			if (aktuelle_Schritt_Nr == Schrittzahl - 1) {
				//cout<<"letzter Punkt\n";
				//cout<<"Pixelnummer:  "<<Pixelnummer<<"\n\n";
				if (Pixelnummer != -1 &&
					Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt.Betrag_ausgeben() == 0) {
						//cout<<"MessungNR: "<<MessungNR<<"\n";
						//cout<<"fehlenden hinteren Durchstoßpunkt setzen\n";
						Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + aml_it->m_Erdradius,
								Punkt_Laenge,
								Punkt_Breite,
								Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(0),
								Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(1),
								Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(2));
				}
			}
			//sleep(1);
		}
		//if(MessungNR==AM_L.size()-1) //nach allen Limbmessungen
		//{
		//    AMF.in_Datei_speichern("/tmp/mlangowski/0/AMF_nach_LOS_Limb.txt");
		//    Tau_LOS_Limb_Matrix.in_Datei_speichern("/tmp/mlangowski/0/TAU_nach
		//}
		//cout<<"Raytracing Limb LOS zu Ende.\n";
		// RAYTRACINGSCHLEIFE LOS ENDE /////////////////////////////////////////

		//cout<<t_Limb_LOS_delta<<"\n";

		////////////////////////////////////////////////////////////////////////
		// LOS Raytracing für Limbmessungen abgeschlossen   9.9.2010
		// Debugging 30.9.2010 Plots für Tau und AMF sehen ok aus
		////////////////////////////////////////////////////////////////////////
		// zusätzliche Abdämpfung durch LFS (Fehler, die Hier gefunden werden,
		// sind auch in Nadir zu berichtigen(beide fast gleich)
		////////////////////////////////////////////////////////////////////////
		//cerr<<"Limb LFS\n";

		//Nun sind die Durchstoßpunkte ja bekannt...
		//also kann grid nochmal ausgeben werden
		//if(MessungNR==210)
		//{  Grid.In_Datei_Ausgeben("/tmp/mlangowski/0/Gitter_nach_LOS.txt");
		//}

		for (int GitterpunktNR = 0; GitterpunktNR < Grid.m_Anzahl_Punkte; GitterpunktNR++) {
			// Schleife über alle Gitterpunkte
			const double Epsilon_double_precision = 2.e-14; //etwa 100 mal größer
			//AMF wird nun durch *=exp(-Tau_LFS) erhöht,
			//falls AMF==0, muss AMF nicht betrachtet werden
			if (abs(AMF(MessungNR, GitterpunktNR)) < Epsilon_double_precision) {
				// == 0 bei double, prüfen, ob korrekt
				continue;
			}
			//cout<<"GitterpunktNR"<<GitterpunktNR<<"\n";
			//cout<<"AMF(MessungNR,GitterpunktNR)"<<AMF(MessungNR,GitterpunktNR)<<"\n";
			//cout<<"MessungNR:"<<MessungNR<<"\n";
			// Alle hier noch übrigen Gitterelemente werden
			// von der LOS geschnitten
			//km  soviel, wie gerade Werte in Tabelle da sind.....
			double TOA_LFS = Konf.m_TOA;
			//Werte, die drüber liegen werden ignoriert
			// hier geht wohl auch 100.....
			// werte zwischen 50 und 200 sind vorhanden
			//Der Startpunkt liegt zwischen den beiden Durchstoßpunkten
			MPL_Vektor Start_Punkt_Polar(3);
			//Wahl des Mittelpunkts des Gitterelements zwischen beiden
			//Durchstoßpunkten der LOS wir mitteln i Kugelkoordinaten und
			//rechnen dann in Karthesische zurück
			MPL_Vektor VordererPunkt_Polar(3);
			MPL_Vektor HintererPunkt_Polar(3);
			Umwandlung_Karthesisch_in_Kugel(Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt(0),
					Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt(1),
					Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt(2),
					VordererPunkt_Polar(0),
					VordererPunkt_Polar(1),
					VordererPunkt_Polar(2));
			Umwandlung_Karthesisch_in_Kugel(Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt(0),
					Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt(1),
					Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt(2),
					HintererPunkt_Polar(0),
					HintererPunkt_Polar(1),
					HintererPunkt_Polar(2));
			Start_Punkt_Polar
				= (HintererPunkt_Polar + VordererPunkt_Polar) * 0.5;
			MPL_Vektor Start_Punkt(3);
			Umwandlung_Kugel_in_Karthesisch(Start_Punkt_Polar(0),
					Start_Punkt_Polar(1), Start_Punkt_Polar(2),
					Start_Punkt(0), Start_Punkt(1), Start_Punkt(2));
			// cout<<"Startpunkt(0): "<<Start_Punkt(0)<<"\t"
			//       <<"Startpunkt(1): "<<Start_Punkt(1)<<"\t"
			//       <<"Startpunkt(2): "<<Start_Punkt(2)<<"\n";
			//MPL_Vektor Start_Punkt(0.5 * (Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt
			//				+ Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt));
			//vom Startpunkt aus gehts in Richtung Sonne voran
			MPL_Vektor Sonne_normal(3); //wie bei LOS
			Umwandlung_Kugel_in_Karthesisch(1,       // 1 km lang damit normiert
					aml_it->m_Sonnen_Longitude,
					aml_it->m_Deklination,
					Sonne_normal(0), Sonne_normal(1), Sonne_normal(2));
			//Einmal den Sonnenzenitwinkel ausrechnen. ist dieser größer als
			//90°, so liegt der Gitterpunkt im Dunkeln, und muss 0 Gesetzt
			//werden
			MPL_Vektor Startpunkt_normiert(Start_Punkt);
			Startpunkt_normiert.Normieren();
			double Cos_SZA_LFS = Startpunkt_normiert * Sonne_normal;
			double b = Start_Punkt.Betrag_ausgeben(); // = r_E + H
			if (Cos_SZA_LFS < 0. &&
					/* check if LFS goes below BOA (default 50 km) */
					b * std::sqrt(1. - Cos_SZA_LFS*Cos_SZA_LFS)
						< aml_it->m_Erdradius + Konf.m_BOA) {
				//cout<<"Limb LFS Gitterpunkt im Dunkeln (SZA>90°)\n";
				AMF(MessungNR, GitterpunktNR) = 0;
				continue;
			}

			// double Punkt_Hoehe=Start_Punkt.Betrag_ausgeben()-aml_it->m_Erdradius;
			// ??? wofür stand das mal

			//Abschätzung der Strecke zwischen Startpunkt und TOA_LFS in km
			// hier kommt der Sonnenzenitwinkel ins Spiel
			// mit Cosinussatz a^2=b^2+c^2-2bc cos(alpha)
			// a ist gesucht alpha=SZA, b=r_E+H, c=r_E+TOA
			double c = TOA_LFS + aml_it->m_Erdradius;
			// Cos(SZA) gerade aus Skalarprodukt des Einheitsvektors in
			// Sonnenrichtung und des Einheitsvektors in Startpunktrichtung
			// bestimmbar
			double Sehnenlaenge_LFS = sqrt(b * b + c * c + 2.0 * b * c * Cos_SZA_LFS);
			////////////////////////////////////////////////////////////////////
			// GESCHWINDIGKEITS/GENAUIGKEITSBESTIMMENDER PARAMETER
			int Schrittzahl2 = 1000;
			// bei TOA_LFS=100.0 und 1000 Schritten ist der
			// Höhenunterschied zweier Punkte etwa 300m
			////////////////////////////////////////////////////////////////////
			double Schrittlaenge = Sehnenlaenge_LFS / ((double)Schrittzahl2);
			//cout<<"Schrittlaenge_LFS_LIMB: "<<Schrittlaenge<<"\n";
			double Tau_LFS = 0.0;
			// RAYTRACINGSCHLEIFE LIMB LFS        /////////////////////////////

			// SPEED Diese Schleife ist ein erheblicher Zeitschritt des Programms
			for (int Schritt = 0; Schritt < Schrittzahl2; Schritt++) {
				MPL_Vektor aktueller_Punkt(3);  //Bei Schrittzahl 10000 3s

				//Bei Schrittzahl 10000 11s
				aktueller_Punkt = Start_Punkt
					+ Sonne_normal * Schrittlaenge * (double) Schritt;
				double Punkt_Hoehe;
				double Punkt_Radius;
				double Punkt_Laenge;
				double Punkt_Breite;
				// Längen und Breitengrad, Höhe und Radius berechnen
				//Bei Schrittzahl 10000 4s
				Umwandlung_Karthesisch_in_Kugel(aktueller_Punkt(0),
						aktueller_Punkt(1), aktueller_Punkt(2),
						Punkt_Radius, Punkt_Laenge, Punkt_Breite);

				// Berechne Hoehe
				Punkt_Hoehe = Punkt_Radius - aml_it->m_Erdradius;
				if (Punkt_Hoehe > TOA_LFS) {
					//keine zusätzliche Absorption von LFS für diese Punkte....
					//(Nur Emission)
					break;
					//Da vom Punkt zur Sonne gelaufen wird, sond alle anderen
					//Punkte ab hier auch höher, deshalb break
				}

				// Testen auf Sonnenzenizenitwinkel unter 90 grad
				if ((Punkt_Breite > 90.0 + aml_it->m_Deklination)
						|| (Punkt_Breite < -90.0 + aml_it->m_Deklination)) {
					continue;
				}
				double BOA_LFS = Konf.m_BOA; //bottom of atmosphere
				if (Punkt_Hoehe < BOA_LFS) {
					cout << "50 km unterschritten bei LFS, "
						 << "das Licht wird schon bei LFS absorbiert\n";
					AMF(MessungNR, GitterpunktNR) = 0;
					break; //weniger al 0 geht nicht
				}
				//cout<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
				//usleep(300000);

				// Not To do Berechne Breite ...wozu auch...
				// brauch ja keinen Streuwinkel
				// Wirkungsquerschnitte schon interpoliert....
				// nun n aus n-höhen information interpolieren
				//  wieder interpolieren....
				// (langsam, hier aber nicht so sehr geschwindigkeitsbestimmend)
				//time(&t_interpolieren_start);
				//Bei Schrittzahl 10000 1,5s
				interpolieren(M_Atmo_Dichten, 0, 1, Punkt_Hoehe, V_Atmo_Dichten(0));
				//Bei Schrittzahl 10000 1,5s
				interpolieren(M_Atmo_Dichten, 0, 2, Punkt_Hoehe, V_Atmo_Dichten(1));
				//time(&t_interpolieren_ende);
				//t_interpolieren_delta_gesamt +=
				//   t_interpolieren_ende - t_interpolieren_start;
				// Berechne Tau=ds*n*sigma * Faktor
				//Woher kommt die 100000-> Schrittlänge wird in km angegeben
				//die Dichten in 1/cm^3 und die WirkungsQS in cm^2
				//Tau ist Einheitenlos
				//Skalarprodukt beider Vektoren

				Tau_LFS += Schrittlaenge * 100000.0
					* (V_Atmo_Dichten * V_Atmo_Wqs_abs);
				//std::cerr << "WQ = " << V_Atmo_Dichten * V_Atmo_Wqs_abs << std::endl;
				//std::cerr << "Tau = " << 100000.0 * Schrittlaenge * V_Atmo_Dichten * V_Atmo_Wqs_abs << std::endl;
			} // Ende for Schritt
			// ENDE RAYTRACINGSCHLEIFE LIMB LFS /////////////////////////////
			// Multipliziere AMF mit exp(-Tau)
			// Bei LFS zählt nur Tau_gesamt, deshalb ist der weg,
			// wie tau erhöht wird egal
			//std::cerr << "Tau_LFS = " << Tau_LFS << std::endl;
			AMF(MessungNR, GitterpunktNR) *= exp(-Tau_LFS);

		} // ende Schleife über alle Gitterpunkte für zusätzliche Abdämpfung LFS
		  // (for Gitterpunktnummer)

	}//ende for MessungNR
	// Ende Limb  RAYTRACING////////////////////////////////////////////////////
	//START kleine Statistik für Winkelabweichungen bei Limb   /////////////////
	Wstat.Statistik_auf_Bildschirm_ausgeben();
	// ENDE kleine Statistik für Winkelabweichungen bei Limb   /////////////////

	////////////////////////////////////////////////////////////////////////////
	// START NADIR RAYTRACING //////////////////////////////////////////////////
	//
	//cerr<<"NADIR Raytracing\n";
	vector<Ausgewertete_Messung_Nadir>::iterator amn_it;

	for (amn_it = AM_N.begin(); amn_it != AM_N.end(); ++amn_it) {
		//cout<<"MessungNR: "<<MessungNR<<"\n";
		// ACHTUNG 2 zählweisen
		// in Matrizen Messungnummer in Quellvektoren NadirmessungNr
		unsigned int MessungNR = AM_L.size() + distance(AM_N.begin(), amn_it);
		// DAS STAMMT HIER AUS DER FUNKTION FÜR LIMB(einfach nur copy+paste) ///
		// Für jede Teilchensorte können mehrere Linien ausgewertet werden.
		// Welche gerade verwendet wird, wird ermittelt, Aus der Wellenlänge
		// des Übergangs sind die (beiden) Wirkungsquerschnitte für die beiden
		// Atmosphärengase zu bestimmen
		MPL_Vektor V_Atmo_Wirkungsquerschnitte(2);
		//Dichten sind Höhenabhängig, Bestimmung dort
		MPL_Vektor V_Atmo_Dichten(2);
		// Wellenlängen sind in Atmo-Datei monoton ansteigend
		// -> Quicksearch Algorithmus für Interpolation(schnell)
		//int interpolieren(MPL_Matrix M,int x_Spalte,int y_Spalte,
		//  double x_Wert_des_gesuchten_Wertes, double& gesuchter_Wert);
		interpolieren(M_Atmo_Wirkungsquerschnitte, 0, 1,
				amn_it->m_Wellenlaenge,
				V_Atmo_Wirkungsquerschnitte(0));
		interpolieren(M_Atmo_Wirkungsquerschnitte, 0, 2,
				amn_it->m_Wellenlaenge,
				V_Atmo_Wirkungsquerschnitte(1));
		//ENDE DAS STAMMT HIER AUS DER FUNKTION FÜR LIMB ///////////////////////


		//Die Nadirmessungen in die Matrix eintragen
		//zwei Wege müssen betrachtet werden:
		MPL_Vektor Sat_POS(3);
		MPL_Vektor GP_POS(3);     //GP für Grundpunkt
		MPL_Vektor Sonne_POS(3);
		//  ermitteln der 3 characteristischen Punkte
		Umwandlung_Kugel_in_Karthesisch(
				amn_it->m_Hoehe_Sat + amn_it->m_Erdradius,
				amn_it->m_Longitude_Sat,
				amn_it->m_Latitude_Sat,
				Sat_POS(0), Sat_POS(1), Sat_POS(2));
		Umwandlung_Kugel_in_Karthesisch(amn_it->m_Erdradius,
				amn_it->m_Longitude_Ground,
				amn_it->m_Latitude_Ground,
				GP_POS(0), GP_POS(1), GP_POS(2));
		Umwandlung_Kugel_in_Karthesisch(149.6E6,
				amn_it->m_Sonnen_Longitude,
				amn_it->m_Deklination,
				Sonne_POS(0), Sonne_POS(1), Sonne_POS(2));

		//////////////////////////
		//Satellit-Punkt  //
		/////////////////////////
		//Raytracing LOS-Vom Satelliten bis zum tiefsten Boxenpunkt
		// Zunächst die ganze LOS
		MPL_Vektor Nadir_LOS(GP_POS - Sat_POS);
		// Startpunkt des Strahls finden

		// Der Startpunkt liegt bei der höchsten Boxenhoehe, somit wird AMF für
		// die hohen Boxen nicht gedämpft, da Tau==0
		MPL_Vektor Startpunkt(3);
		Startpunkt = Punkt_auf_Strecke_bei_Radius(Sat_POS, Nadir_LOS,
					 Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1].m_Max_Hoehe
					 + amn_it->m_Erdradius, 0.1);

		//oder auch 200.0 bei Nadir wird nicht noch +Erdradius genommen
		double Max_Hoehe_Absorption = Konf.m_TOA;
		// Endpunkt des Strahls finden
		MPL_Vektor Endpunkt(3);
		// der 0te Gitterpunkt hat die niedrigste Höhe
		double R_Min
			= Grid.m_Gitter[0].m_Min_Hoehe + amn_it->m_Erdradius;

		Endpunkt = Punkt_auf_Strecke_bei_Radius(Sat_POS, Nadir_LOS, R_Min, 0.1);
		//LOS anpassen auf relevantes Höhenintervall
		Nadir_LOS = Endpunkt - Startpunkt;
		MPL_Vektor Nadir_LOS_Einheitsvektor(Nadir_LOS);
		Nadir_LOS_Einheitsvektor.Normieren();
		// Distanz zwischen Start und Endpunkt bestimmen
		double Distanz_in_Atmosphaere = Nadir_LOS.Betrag_ausgeben();
		// Schrittzahl festlegen
		int Schrittzahl = 1000; //<----- bestimmt GESCHWINDIGKEIT des Programms
		// Alle x Schritte soll Streuwinkel neu bestimmt werden....x festlegen
		int Winkelberechnungsfrequenz = 100;
		// Schrittweite bestimmen
		double Schrittweite = Distanz_in_Atmosphaere / (double)Schrittzahl;

		////////////////////////////////////////////////////////////////////////
		//Nadir LOS Raytracingschleife /////////////////////////////////////////
		double Tau_Nadir_LOS = 0.0;
		double Cos_Streuwinkel = -1.0;
		Grid.Alle_Durchstosspunkte_Null_setzen();
		//Startparameter fürs Raytracing, -1 sagt,
		//dass der aktuelle Gitterpunkt nicht bekannt ist
		int Pixelnummer = -1;
		for (int aktueller_Schritt = 0; aktueller_Schritt < Schrittzahl; aktueller_Schritt++) {
			//   TODOLIST //////////////////////////////////////
			// aktuellen Punkt berechnen
			// Prüfen ob Streuwinkel berechnet sollen...und falls ja Berechnen
			// Hoehe, Länge,Breite bestimmen
			// Dichten der Atmosphärengase bestimmen
			// Tau erhoehen
			// Phasenfunktion bestimmen
			// Gitterpunkt finden und AMF erhöhen
			////////////////////////////////////////////////////////////

			// aktuellen Punkt berechnen
			MPL_Vektor aktueller_Punkt(3);
			aktueller_Punkt = Startpunkt
				+ aktueller_Schritt * Schrittweite * Nadir_LOS_Einheitsvektor;

			// Prüfen ob Streuwinkel berechnet sollen...und falls ja Berechnen
			if (aktueller_Schritt % Winkelberechnungsfrequenz == 0) {
				//Erklärung siehe Limb
				MPL_Vektor Sonne_normal(Sonne_POS);
				Sonne_normal.Normieren();
				Cos_Streuwinkel = -1 * Sonne_normal * Nadir_LOS_Einheitsvektor;
			}

			// Hoehe, Länge,Breite bestimmen  //aktueller Punkt AP
			double AP_R, AP_Hoehe, AP_Laenge, AP_Breite;
			Umwandlung_Karthesisch_in_Kugel(aktueller_Punkt(0),
											aktueller_Punkt(1),
											aktueller_Punkt(2),
											AP_R, AP_Laenge, AP_Breite);
			AP_Hoehe = AP_R - amn_it->m_Erdradius;
			if (AP_Hoehe >= Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1].m_Max_Hoehe) {
				//Punkt ausserhalb der Maximalen Gitterhoehe
				continue;
			}

			if (AP_Hoehe <= Max_Hoehe_Absorption) {
				// Dichten der Atmosphärengase bestimmen
				interpolieren(M_Atmo_Dichten, 0, 1, AP_Hoehe, V_Atmo_Dichten(0));
				interpolieren(M_Atmo_Dichten, 0, 2, AP_Hoehe, V_Atmo_Dichten(1));
			}
			// Tau erhoehen (Tau ist Einheitenlos)  (100000 aus Umwandlung km zu cm)
			Tau_Nadir_LOS += Schrittweite * 100000
				* (V_Atmo_Dichten * V_Atmo_Wirkungsquerschnitte);

			// Check, ob Punkt in Dunkelheit(extrem unwahrscheinlich bei Nadir,
			// aber sicher ist sicher)
			if ((AP_Breite > 90.0 + amn_it->m_Deklination)
					|| (AP_Breite < -90.0 + amn_it->m_Deklination)) {
				if (aktueller_Schritt == Schrittzahl - 1) {
					//cout<<"MessungNR: "<<MessungNR<<"\n";
					//cout<<"letzter Punkt und im dunkeln(Nadir)\n";
					if (Pixelnummer != -1) {
						if (Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt.Betrag_ausgeben() == 0) {
							cout << "MessungNR: " << MessungNR << "\n";
							cout << "fehlenden hinteren Durchstoßpunkt setzen\n";
							Umwandlung_Kugel_in_Karthesisch(
									AP_Hoehe + amn_it->m_Erdradius,
									AP_Laenge,
									AP_Breite,
									Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(0),
									Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(1),
									Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(2));
						}
					}
				}
				continue;
			}

			// Phasenfunktion bestimmen
			double E1 = 0., E2 = 1.;
			unsigned int PeakNr;
			for (PeakNr = 0; PeakNr < Spezies_Fenster.m_Wellenlaengen.size(); PeakNr++) {
				if (amn_it->m_Wellenlaenge == Spezies_Fenster.m_Wellenlaengen[PeakNr]) {
					E1 = Spezies_Fenster.m_Liniendaten[PeakNr].m_E1;
					E2 = Spezies_Fenster.m_Liniendaten[PeakNr].m_E2;
					break;
				}
			}
			if (PeakNr == Spezies_Fenster.m_Wellenlaengen.size()) {
				//eigentlich kann das gar nicht passieren
				cout << "E1 und E2 konnten nicht gefunden werden..."
					 << "mysteriöser bug\n";
				E1 = 0.0;
				E2 = 1.0;
			}
			//Phasenfunktion aus Streuwinkel bestimmen und
			//mit AMF multiplizieren...in Funktion
			double Phasenfunktion = 0.75 * E1
				* (Cos_Streuwinkel * Cos_Streuwinkel + 1) + E2;
			// Gitterpunkt finden und AMF erhöhen;
			// Durchstoßpunkte der Gitterelemente finden (siehe LFS)

			int myerr = Pixel_finden_und_AMF_erhoehen_LOS(AMF, Grid, MessungNR,
					Pixelnummer,
					Schrittweite, Tau_Nadir_LOS,
					AP_Hoehe, amn_it->m_Erdradius,
					AP_Laenge, AP_Breite,
					Phasenfunktion, Tau_LOS_Matrix);
			if (myerr == 4) {
				IERR = 1;
				cout << "mysteriöser fall beendet Programm\n";
				MPL_Matrix dummy(1, 1);
				return dummy;
			}

			if (aktueller_Schritt == Schrittzahl - 1) {
				//cout<<"letzter Punkt\n";
				//cout<<"Pixelnummer:  "<<Pixelnummer<<"\n\n";
				if (Pixelnummer != -1) {
					if (Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt.Betrag_ausgeben() == 0) {
						cout << "MessungNR: " << MessungNR << "\n";
						cout << "fehlenden hinteren Durchstoßpunkt setzen\n";
						Umwandlung_Kugel_in_Karthesisch(
								AP_Hoehe + amn_it->m_Erdradius,
								AP_Laenge,
								AP_Breite,
								Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(0),
								Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(1),
								Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(2));
					}
				}
			}
		}
		//Ende Nadir LOS Raytracingschleife ////////////////////////////////////
		////////////////////////////////////////////////////////////////////////
		// nun die gesamte Tau matrix und AMF MATRIx plotten
		// (achtung AMF für Nadir noch nicht gedämpft)
		//if(MessungNR==AMF.m_Zeilenzahl-1)
		//{
		//    AMF.in_Datei_speichern("/tmp/mlangowski/0/AMF_nach_LOS_Nadir.txt");
		//    Tau_LOS_Matrix.in_Datei_speichern("/tmp/mlangowski/0/TAU_nach_LOS_Nadir.txt");
		//    //schaut gut aus
		//}
		/////////////////////////
		//Punkt-Sonne  //
		/////////////////////////
		//if(MessungNR==260)
		//{
		//    Grid.In_Datei_Ausgeben("/tmp/mlangowski/0/Gitter_260.txt");
		//}

		////////////////////////////////////////////////////////////////////////
		// zusätzliche Abdämpfung durch LFS
		// DIESER ABSCHNITT IST BIS AUF DIE VARIABLENNAMEN IDENTISCH MIT
		// DEM LFS ABSCHNITT BEI LIMB ( eventuell Funktion schreiben)
		////////////////////////////////////////////////////////////////////////
		for (int GitterpunktNR = 0; GitterpunktNR < Grid.m_Anzahl_Punkte; GitterpunktNR++) {
			// Schleife über alle Gitterpunkte
			const double Epsilon_double_precision = 2e-14; //etwa 100 mal größer
			//AMF wird nun durch *=exp(-Tau_LFS) erhöht,
			//falls AMF==0, muss AMF nicht betrachtet werden
			if (abs(AMF(MessungNR, GitterpunktNR)) < Epsilon_double_precision) {
				// == 0 bei double, prüfen, ob korrekt
				continue;
			}

			// Alle hier noch übrigen Gitterelemente werden von der LOS
			// geschnitten
			double TOA_LFS = Konf.m_TOA;
			//km  soviel, wie gerade Werte in Tabelle da sind
			//Der Startpunkt liegt zwischen den beiden Durchstoßpunkten
			MPL_Vektor Start_Punkt_Polar(3);
			//Wahl des Mittelpunkts des Gitterelements zwischen
			//beiden Durchstoßpunkten der LOS
			// wir mitteln i Kugelkoordinaten und rechnen dann in Karthesische
			// zurück
			MPL_Vektor VordererPunkt_Polar(3);
			MPL_Vektor HintererPunkt_Polar(3);
			Umwandlung_Karthesisch_in_Kugel(
					Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt(0),
					Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt(1),
					Grid.m_Gitter[GitterpunktNR].m_vorderer_Durchstosspunkt(2),
					VordererPunkt_Polar(0),
					VordererPunkt_Polar(1),
					VordererPunkt_Polar(2));
			Umwandlung_Karthesisch_in_Kugel(
					Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt(0),
					Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt(1),
					Grid.m_Gitter[GitterpunktNR].m_hinterer_Durchstosspunkt(2),
					HintererPunkt_Polar(0),
					HintererPunkt_Polar(1),
					HintererPunkt_Polar(2));
			Start_Punkt_Polar
				= (HintererPunkt_Polar + VordererPunkt_Polar) * 0.5;
			MPL_Vektor Start_Punkt(3);
			Umwandlung_Kugel_in_Karthesisch(
					Start_Punkt_Polar(0), Start_Punkt_Polar(1), Start_Punkt_Polar(2),
					Start_Punkt(0), Start_Punkt(1), Start_Punkt(2));
			MPL_Vektor Sonne_normal(Sonne_POS); //wie bei LOS
			Sonne_normal.Normieren();
			//Einmal den Sonnenzenitwinkel ausrechnen. ist dieser größer als
			//90° so liegt der Gitterpunkt im Dunkeln, und muss 0 Gesetzt
			//werden
			MPL_Vektor Startpunkt_normiert(Start_Punkt);
			Startpunkt_normiert.Normieren();
			double Cos_SZA_LFS = Startpunkt_normiert * Sonne_normal;
			//if(MessungNR<500)
			//{cout<<Cos_SZA_LFS<<"\n";}
			if (Cos_SZA_LFS < 0) {
				//cout<<"Limb LFS Gitterpunkt im Dunkeln (SZA>90°)\n";
				AMF(MessungNR, GitterpunktNR) = 0;
				continue;
			}
			// hier kommt der Sonnenzenitwinkel ins Spiel
			// mit Cosinussatz a^2=b^2+c^2-2bc cos(alpha)
			// a ist gesucht
			// alpha=SZA, b=r_E+H, c=r_E+TOA
			double b = Start_Punkt.Betrag_ausgeben();
			double c = TOA_LFS + amn_it->m_Erdradius;
			double Sehnenlaenge_LFS = sqrt(b * b + c * c + 2.0 * b * c * Cos_SZA_LFS);
			////////////////////////////////////////////////////////////////////
			// GESCHWINDIGKEITS/GENAUIGKEITSBESTIMMENDER PARAMETER
			int Schrittzahl2 = 1000;
			////////////////////////////////////////////////////////////////////
			double Schrittlaenge = Sehnenlaenge_LFS / ((double)Schrittzahl2);
			double Tau_LFS = 0.0;
			// RAYTRACINGSCHLEIFE NADIR LFS          ///////////////////////////
			for (int Schritt = 0; Schritt < Schrittzahl2; Schritt++) {
				MPL_Vektor aktueller_Punkt(3);
				aktueller_Punkt = Start_Punkt
					+ Sonne_normal * Schrittlaenge * (double) Schritt;
				double Punkt_Hoehe;
				double Punkt_Radius;
				double Punkt_Laenge;
				double Punkt_Breite;
				// Längen und Breitengrad, Höhe und Radius berechnen
				//Bei Schrittzahl 10000 4s
				Umwandlung_Karthesisch_in_Kugel(aktueller_Punkt(0),
						aktueller_Punkt(1), aktueller_Punkt(2),
						Punkt_Radius, Punkt_Laenge, Punkt_Breite);

				// Berechne Hoehe
				Punkt_Hoehe = Punkt_Radius - amn_it->m_Erdradius;
				if (Punkt_Hoehe > TOA_LFS) {
					//keine zusätzliche Absorption von LFS für diese Punkte...
					//(Nur Emission)
					break;
					//Da vom Punkt zur Sonne gelaufen wird, sond alle anderen
					//Punkte ab hier auch höher, deshalb break
				}
				// Testen auf Sonnenzenizenitwinkel unter 90 grad
				if ((Punkt_Breite > 90.0 + amn_it->m_Deklination) ||
						(Punkt_Breite < -90.0 + amn_it->m_Deklination)) {
					continue;
				}
				double BOA_LFS = Konf.m_BOA; //bottom of atmosphere
				if (Punkt_Hoehe < BOA_LFS) {
					cout << "50 km unterschritten bei LFS, "
						 << "das Licht wird schon bei LFS absorbiert\n";
					AMF(MessungNR, GitterpunktNR) = 0;
					break; // weniger al 0 geht nicht
				}
				//time(&t_interpolieren_start);
				//Bei Schrittzahl 10000 1,5s
				interpolieren(M_Atmo_Dichten, 0, 1, Punkt_Hoehe, V_Atmo_Dichten(0));
				//Bei Schrittzahl 10000 1,5s
				interpolieren(M_Atmo_Dichten, 0, 2, Punkt_Hoehe, V_Atmo_Dichten(1));
				//time(&t_interpolieren_ende);
				//t_interpolieren_delta_gesamt +=
				//  t_interpolieren_ende - t_interpolieren_start;
				// Berechne Tau=ds*n*sigma *Faktor
				//Woher kommt die 100000-> Schrittlänge wird in km angegeben
				//die Dichten in 1/cm^3 und die WirkungsQS in cm^2 Tau ist
				//Einheitenlos //Skalarprodukt beider Vektoren

				Tau_LFS += Schrittlaenge * 100000.0
					* (V_Atmo_Dichten * V_Atmo_Wirkungsquerschnitte);
			} // Ende for Schritt
			// Multipliziere AMF mit exp(-Tau)
			AMF(MessungNR, GitterpunktNR) *= exp(-Tau_LFS);
			// ENDE RAYTRACINGSCHLEIFE NADIR LFS /////////////////////////////
		}   // ende Schleife über alle Gitterpunkte
			//für zusätzliche Abdämpfung LFS Nadir

	}// Ende for MessungNR
	// ENDE NADIR RAYTRACING ///////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// LIMB UND NADIR RAYTRACING ABGESCHLOSSEN AMF ZURÜCKGEBEN

//    AMF.in_Datei_speichern("/tmp/mlangowski/0/AMF_fertig.txt");
//    Tau_LOS_Matrix.in_Datei_speichern("/tmp/mlangowski/0/TAU_fertig.txt");

	// char buf[256];
	// string dum;
	// sprintf(buf,"Laufzeit für LIMB_LOS in Sekunden:\t %ld",t_Limb_LOS_delta);
	// dum=buf;
	// cout<<dum<<"\n";
	// sprintf(buf,"Laufzeit für interpolationen:\t %ld",t_interpolieren_delta_gesamt);
	// dum=buf;
	// cout<<dum<<"\n";
	return AMF;
}// Ende Luftmassenfaktoren_Matrix_aufbauen()



////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  interpolieren()
////////////////////////////////////////////////////////////////////////////////
int interpolieren(MPL_Matrix &M, int x_Spalte, int y_Spalte,
		double x_Wert_des_gesuchten_Wertes, double &gesuchter_Wert)
{
	// TODO Die Funktion ist, falls sie sehr oft aufgerufen wird eine ziemliche
	// Bremse; x werte müssen monoton steigen


	//Die x und y Werte stehen in den ersten beiden Spalten der Matrix, die
	//oben stehende Funktion wird 1:1 übernommen, sodass die kommentierung hier
	//sparsamer ausfällt (0te Spalte x-Werte, 1 bzw 2ten Spalte y-Werte)
	int n = M.m_Zeilenzahl;

	int Index_Anfang = 0;
	int Index_Ende = n - 1;

	/* constant at boundaries */
	if (x_Wert_des_gesuchten_Wertes < M(Index_Anfang, x_Spalte)) {
		gesuchter_Wert = M(Index_Anfang, y_Spalte);
		return 1;
	}
	if (x_Wert_des_gesuchten_Wertes > M(Index_Ende, x_Spalte)) {
		gesuchter_Wert = M(Index_Ende, y_Spalte);
		return 1;
	}
	//cout<<"M(0,x_Spalte): "<<M(0,x_Spalte)<<"\n";
	//cout<<"M(n-1,x_Spalte): "<<M(n-1,x_Spalte)<<"\n";
	//cout<<"x_Wert_des_gesuchten_Wertes: "<<x_Wert_des_gesuchten_Wertes<<"\n";
	while (Index_Anfang + 1 < Index_Ende) {

		int Index = (Index_Anfang + Index_Ende) / 2;
		//cout<<"Index:"<<Index<<"\n";
		//cout<<"Index_anfang:"<<Index_Anfang<<"\n";
		//cout<<"Index_ende:"<<Index_Ende<<"\n";

		if (M(Index, x_Spalte) > x_Wert_des_gesuchten_Wertes) {
			Index_Ende = Index;
		} else if (M(Index, x_Spalte) < x_Wert_des_gesuchten_Wertes) {
			Index_Anfang = Index;
		} else {
			gesuchter_Wert = M(Index, y_Spalte);
			return 0;
		}
	}// while(true)


	//Anfangs und Endindex sind nun bekannt und um 1 unterschiedlich,
	//sodass die eigentliche Interpolation druchgeführt werden kann
	double I2 = ((double)(x_Wert_des_gesuchten_Wertes - M(Index_Anfang, x_Spalte)))
				/ ((double)(M(Index_Ende, x_Spalte) - M(Index_Anfang, x_Spalte)));
	double I1 = 1.0 - I2;
	//cout<<"Index_Anfang:"<<Index_Anfang<<"   Index_Ende:"<<Index_Ende<<"\n";
	//cout<<"x_Anfang:"<<M(Index_Anfang, x_Spalte)<<"   x_Ende:"<<M(Index_Ende, x_Spalte)<<"\n";
	//cout<<"I2: "<<I2<<"   I1: "<<I1<<"\n";
	//Bsp x_ges=5,1 ,x_A=5 x_E=6  ->I2=0,1, I1=0,9

	//Ergebnisrückgabe
	gesuchter_Wert = I1 * M(Index_Anfang, y_Spalte) + I2 * M(Index_Ende, y_Spalte);
	return 0;
}//Ende interpolieren

////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Pixel_finden_und_AMF_erhoehen_LOS();
////////////////////////////////////////////////////////////////////////////////
int Pixel_finden_und_AMF_erhoehen_LOS(MPL_Matrix &AMF, Retrievalgitter &Grid,
		const int &MessungNr, int &Pixelnummer, const double &Schrittlaenge,
		const double &Tau_LOS, const double &Punkt_Hoehe,
		const double &Erdradius,
		const double &Punkt_Laenge, const double &Punkt_Breite,
		const double &Phasenfunktion, MPL_Matrix &Tau_LOS_Limb_Matrix)
{
	/***************************************************************************
	1.Die Funktion ruft zuerst für den aktuellen Gitterpunkt und dessen 8
	Nachbarn die Suchfunktion (Punkt_Pruefen_und_ggf_AMF_erhoehen) auf, welche
	prüft, ob sich der aktuelle Punkt(Hoehe,Breite) innerhalb eines dieser
	Gitterelemente befindet.
	2.Falls der letzte gesuchte Punkt in keinem Gitterelement lag, so ist es
	ziemlich wahrscheinlich, dass der aktuelle Punkt wieder nicht in einem
	Gitterelement liegt, sodass zunächst die Grenzen des Gitters überprüft
	werden, da der nächste Schritt viel Zeit benötigt
	3.Als letzten Schritt werden alle Gitterelemente abgesucht

	Die Fälle 1,2 und 3 sollten alle Möglichkeiten abdecken
	Der Rückgabewert der Funktion kann zum debuggen und erheben von Statistiken
	genutzt werden(z.b. einfach zählen, wie oft jeder Fall so vorkommt) Es gibt
	also die Rückgabewerte 1,2,3 für die 3 Fälle
	Falls keiner der Fälle auftritt, so wird 4 zurückgegeben....Dies weist
	darauf hin, das etwas nicht stimmt (z.B, das Gitter am Randf nicht dicht
	ist) und DARF eigentlich NICHT vorkommen

	Der Einzige IO Parameter, der verändert wird ist ist die Pixelnummer(ein
	des alten GP ;aus:des aktuellen Gitterpunktes)
	 **************************************************************************/
	// Die Funktion ist leider unschön redundant und lang, ob der langen
	// funktionsparameterlisten..(funktioniert aber) scheint aber hinreichend
	// schnell zu sein

	// Die Sache mit den Durchstoßpunkten scheint noch nicht richtig
	// zu funktionieren...bzw muss besser getestet werden
	// Evtl. Vektor bauen, in für einen Strahlengang beim wechsel der Punkte
	// nur die Gitterpunkte notiert werden und verfolgen Der sollte dann auch
	// die Richtung des Strahlengangs richtig angeben

	// FALL 1
	if (Pixelnummer != -1) {
		//cout<<"möglicher Fall 1\n";
		//sleep(1);
		//cout<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
		//cout<<"Punkt_Laenge: "<<Punkt_Laenge<<"\n";
		//cout<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
		// Suche zunächst Punkt, dann 8 Nachbarn des Punktes ab
		//Punkt selbst
		//FALL 1a
		if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
					Pixelnummer, Pixelnummer, Schrittlaenge, Tau_LOS,
					Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
					Tau_LOS_Limb_Matrix) == true) {
			// hier keine weiteren Veränderungen...Fall 1a kommt sehr häufig vor
			// jetzt muss das hier doch stehen, weil manchmal hintere
			// Durchstoßpunkte fehlen wenn man da den Spezialfall findet,
			// könnte man das hier wegnehmen....(scheint aber auch nicht zu
			// sehr zu bremsen)
			Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
					Punkt_Laenge,
					Punkt_Breite,
					Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(0),
					Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(1),
					Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(2));

			return 1;
		}
		//FALL 1b
		//8 Nachbarn absuchen
		if ((Grid.m_Gitter[Pixelnummer].m_Index_oberer_Nord_Nachbar != -1)) {
			//Nord oben
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_oberer_Nord_Nachbar;
			// Pixelnummer wird ja im nächsten Schritt angepasst
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {
				//if(MessungNr==28)
				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_Nord_Nachbar != -1)) {
			//Nord
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_Nord_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {

				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_unterer_Nord_Nachbar != -1)) {
			//Nord unten
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_unterer_Nord_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {

				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_oberer_Nachbar != -1)) {
			// oben
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_oberer_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {

				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_unterer_Nachbar != -1)) {
			// unten
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_unterer_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {

				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_oberer_Sued_Nachbar != -1)) {
			//Süd oben
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_oberer_Sued_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {
				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_Sued_Nachbar != -1)) {
			//Süd
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_Sued_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {
				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		if ((Grid.m_Gitter[Pixelnummer].m_Index_unterer_Sued_Nachbar != -1)) {
			//Süd unten
			int PN_Test = Grid.m_Gitter[Pixelnummer].m_Index_unterer_Sued_Nachbar;
			int Pixno_old = Pixelnummer;
			if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr,
						PN_Test, Pixelnummer, Schrittlaenge, Tau_LOS,
						Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
						Tau_LOS_Limb_Matrix) == true) {
				// Für den Testpunkt ist das vermutlich der vordere
				// Durchstoßpunkt...also der erste Punkt auf der LOS
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(0),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(1),
						Grid.m_Gitter[PN_Test].m_vorderer_Durchstosspunkt(2));
				//Für den alten Gitterpunkt ist das quasi einer hinter dem
				//letzten Punkt also der hintere Durchstoßpunkt;
				Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
						Punkt_Laenge,
						Punkt_Breite,
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(0),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(1),
						Grid.m_Gitter[Pixno_old].m_hinterer_Durchstosspunkt(2));
				return 1;
			}
		}
		// FALL 1c Vorgänger war bekannt, aber Punkt trotzdem nicht gefunden
		// Na dann ist der Vorgänger erstmal zu Ende
		Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
				Punkt_Laenge,
				Punkt_Breite,
				Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(0),
				Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(1),
				Grid.m_Gitter[Pixelnummer].m_hinterer_Durchstosspunkt(2));
		//cout<<"Fall 1c\n";
		//cout<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
		//cout<<"Punkt_Laenge: "<<Punkt_Laenge<<"\n";
		//cout<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
	}//ENDE if (Pixelnummer!=-1)
	// FALL 2 Punkt garnicht im Pixelgitter?
	if (Punkt_Hoehe > Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1].m_Max_Hoehe) {
		Pixelnummer = -1;
		cout << "Fall2a\n";
		cout << "Punkt_Hoehe: " << Punkt_Hoehe << "\n";
		cout << "Grid.m_Gitter[Grid.m_Anzahl_Punkte].m_Max_Hoehe: "
			 << Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1].m_Max_Hoehe << "\n";
		return 2;
	}  // kein Punkt des Gitters
	// the grid is sorted from north to south,
	// Grid.m_Gitter[0] is the northernmost point
	if (Punkt_Breite > Grid.m_Gitter[0].m_Max_Breite) {
		Pixelnummer = -1;
		//cout<<"Fall2b\n";
		//cout<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
		//cout<<"Grid.m_Gitter[0].m_Max_Breite: "<<Grid.m_Gitter[0].m_Max_Breite<<"\n";
		return 2;
	}  // kein Punkt des Gitters
	// Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1] is the southernmost point
	if (Punkt_Breite < Grid.m_Gitter[Grid.m_Anzahl_Punkte - 1].m_Min_Breite) {
		Pixelnummer = -1;
		//cout<<"Fall2c\n";
		//cout<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
		//cout<<"Grid.m_Gitter[Grid.m_Anzahl_Punkte-1].m_Min_Breite: "
		//  <<Grid.m_Gitter[Grid.m_Anzahl_Punkte-1].m_Min_Breite<<"\n";
		return 2;
	}  // kein Punkt des Gitters
	if (Punkt_Hoehe <= Grid.m_Gitter[0].m_Min_Hoehe) {
		Pixelnummer = -1;
		/*
		cout << "Fall2d\n";
		cout << "Achtung...der Fall verfälscht das Ergebnis stark\n";
		cout << "MessungNr: " << MessungNr << "\n";
		cout << "Punkt_Hoehe: " << Punkt_Hoehe << "\n";
		cout << "Grid.m_Gitter[0].m_Min_Hoehe: "
			 << Grid.m_Gitter[0].m_Min_Hoehe << "\n";
		// */
		return 2;
	}  // kein Punkt des Gitters
	//FALL3 Alle Punkte des Pixelgitters absuchen
	for (int j = 0; j < Grid.m_Anzahl_Punkte; j++) {
		//cout<<"j: "<<j<<"\n";
		if (Punkt_Pruefen_und_ggf_AMF_erhoehen(AMF, Grid, MessungNr, j,
					Pixelnummer, Schrittlaenge, Tau_LOS,
					Punkt_Hoehe, Punkt_Breite, Phasenfunktion,
					Tau_LOS_Limb_Matrix) == true) {
			// ab hier sind j und Pixelnummer gleich!
			//cout<<"Fall3\n";
			//cout<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
			//cout<<"Punkt_Laenge: "<<Punkt_Laenge<<"\n";
			//cout<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
			//cout<<"Pixelnummer:"<<Pixelnummer<<"\n";
			// Für den Testpunkt ist das vermutlich der vordere Durchstoßpunkt
			// ...also der erste Punkt auf der LOS
			Umwandlung_Kugel_in_Karthesisch(Punkt_Hoehe + Erdradius,
					Punkt_Laenge,
					Punkt_Breite,
					Grid.m_Gitter[j].m_vorderer_Durchstosspunkt(0),
					Grid.m_Gitter[j].m_vorderer_Durchstosspunkt(1),
					Grid.m_Gitter[j].m_vorderer_Durchstosspunkt(2));
			//cout<<"Gitterpunkt: "<<j<<"\n";
			//einen alten Messpunkt gibts ja quasi nicht...
			//also auch kein hinterer Durchstoßpunkt
			return 3;
		}
	}// Ende for j
	//FALL 4 Mysteriöse Effekte
	cout << "Mysteriöser Fall in Pixel_finden_und_AMF_erhoehen_LOS\n";
	return 4;
}//Ende Pixel_finden_und_AMF_erhoehen();
////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Punkt_Pruefen_und_ggf_AMF_erhoehen()
////////////////////////////////////////////////////////////////////////////////
bool Punkt_Pruefen_und_ggf_AMF_erhoehen(MPL_Matrix &AMF, Retrievalgitter &Grid,
		const int &MessungNR, const int &PN_Test, int &Pixelnummer,
		const double &Schrittlaenge, const double &Tau_LOS,
		const double &Punkt_Hoehe, const double &Punkt_Breite,
		const double &Phasenfunktion, MPL_Matrix &Tau_LOS_Limb_Matrix)
{
	//ofstream outfile;
	//outfile.open("/home/martin/TEMP/Punkt_pruefen_output.txt",ios_base::out | ios_base::app);
	//outfile<<"Start Punkt_Pruefen_und_ggf_AMF_erhoehen\n";
	//outfile<<"PN_Test: "<<PN_Test<<"\n";
	//outfile<<"Punkt_Hoehe: "<<Punkt_Hoehe<<"\n";
	//outfile<<"Punkt_Breite: "<<Punkt_Breite<<"\n";
	//outfile<<"Grid.m_Gitter[PN_Test].m_Max_Breite: "<<Grid.m_Gitter[PN_Test].m_Max_Breite<<"\n";
	//outfile<<"Grid.m_Gitter[PN_Test].m_Min_Breite: "<<Grid.m_Gitter[PN_Test].m_Min_Breite<<"\n";
	//outfile<<"Grid.m_Gitter[PN_Test].m_Max_Hoehe: "<<Grid.m_Gitter[PN_Test].m_Max_Hoehe<<"\n";
	//outfile<<"Grid.m_Gitter[PN_Test].m_Min_Hoehe: "<<Grid.m_Gitter[PN_Test].m_Min_Hoehe<<"\n";
	//outfile.close();

	// Die Funktion hat einen etwas großen Funktionskopf da könnte man
	// optimieren, falls das Geschwindigkeit bringt Die Funktion wird in
	// Pixel_finden_und_AMF_erhoehen() 10 mal aufgerufen

	if (Grid.m_Gitter[PN_Test].Punkt_in_Gitterpunkt(Punkt_Breite, Punkt_Hoehe) == true) {
		//cout<<"Punkt gefunden\n";
		//Luftmassenfaktorenmatrix erhöhen in cm...Schrittlänge in km gegeben
		//also 10^5 UMRECHNUNGSFAKTOR die ite Limbmessung bildet die ite Zeile
		//der Matrix
		Pixelnummer = PN_Test;
		//cout<<"AMF(MessungNR,Pixelnummer)"<<AMF(MessungNR,Pixelnummer)<<"\n";
		// 100000.0 * Schrittlaenge ~ s_{ij} = part of the LOS within the pixel
		AMF(MessungNR, Pixelnummer) += 100000.0 * Schrittlaenge
			* exp(-Tau_LOS) * Phasenfunktion;
		Tau_LOS_Limb_Matrix(MessungNR, Pixelnummer) = (Tau_LOS);
		//cout<<"MessungNR: "<<MessungNR<<"\n";
		//cout<<"Tau_LOS:"<<Tau_LOS<<"\n";
		//cout<<"exp(-Tau_LOS):"<<exp(-Tau_LOS)<<"\n";
		//cout<<"Schrittlaenge:"<<Schrittlaenge<<"\n";
		//cout<<"Phasenfunktion:"<<Phasenfunktion<<"\n";
		//cout<<"AMF(MessungNR,Pixelnummer)"<<AMF(MessungNR,Pixelnummer)<<"\n";
		//sleep(2);
		return true;
	}
	return false;
}// Ende Punkt_Pruefen_und_ggf_AMF_erhoehen

////////////////////////////////////////////////////////////////////////////////
//Funktionsstart  Punkt_auf_Strecke_bei_Radius
////////////////////////////////////////////////////////////////////////////////
MPL_Vektor Punkt_auf_Strecke_bei_Radius(MPL_Vektor &Streckenstartpunkt,
		MPL_Vektor &Streckenvektor, double Radius, double Genauigkeit)
// int max iterationen...lass ich weg
{
	// ACHTUNG...die Funktion muss richtig aufgerufen werden
	// Die gesuchte Höhe muss auch auf dem Vektor zu finden sein,
	// sonst endlosschleife
	// Iterationsstart
	double Startpunkt_Faktor = 0.5;
	double Veraenderung = 0.5;
	MPL_Vektor aktueller_Vektor(3);
	aktueller_Vektor = Streckenstartpunkt + Startpunkt_Faktor * Streckenvektor;
	double Eps = Genauigkeit; // in km (100m=0.1)
	//Iterationsschleife
	while (!(((aktueller_Vektor.Betrag_ausgeben() - Eps) < Radius)
				&& ((aktueller_Vektor.Betrag_ausgeben() + Eps) > Radius))) {
		//cout<<"Startpunkt_Faktor: "<<Startpunkt_Faktor<<"\n";
		//cout<<"Veraenderung: "<<Veraenderung<<"\n";
		//cout<<"aktueller_Vektor.Betrag_ausgeben()"<<aktueller_Vektor.Betrag_ausgeben()<<"\n";
		//cout<<"Radius "<<Radius<<"\n";
		Veraenderung /= 2;
		if ((aktueller_Vektor.Betrag_ausgeben() + Eps) <= Radius) {
			Startpunkt_Faktor -= Veraenderung;
		} else {
			Startpunkt_Faktor += Veraenderung;
		}
		aktueller_Vektor = Streckenstartpunkt
			+ Startpunkt_Faktor * Streckenvektor;
	}
	return aktueller_Vektor;
}
////////////////////////////////////////////////////////////////////////////////
//ENDE  Punkt_auf_Strecke_bei_Radius
////////////////////////////////////////////////////////////////////////////////

// Prepares a (column) vector to hold the total number density at the grid
// points using the number densities from the measurement points that fall
// into the corresponding grid point.
void prepare_total_density(Retrievalgitter &grid, MPL_Matrix &dens,
		std::vector<Ausgewertete_Messung_Limb> &aml_vec)
{
	std::vector<Gitterpunkt> zero_dens_pts;
	std::vector<Gitterpunkt>::iterator gp;

	for (gp = grid.m_Gitter.begin(); gp != grid.m_Gitter.end(); ++gp) {
		int i = std::distance(grid.m_Gitter.begin(), gp);
		std::vector<double> densities;
		std::vector<Ausgewertete_Messung_Limb>::iterator aml_it;
		/*
		 * finds measurement points that fall into the grid point
		 * and averages their already calculated number densities.
		 */
		for (aml_it = aml_vec.begin(); aml_it != aml_vec.end(); ++aml_it) {
			if (gp->Punkt_in_Gitterpunkt(aml_it->m_Latitude_TP,
						aml_it->m_Hoehe_TP)) {
				densities.push_back(aml_it->total_number_density);
			}
		}
		dens(i) = std::accumulate(densities.begin(), densities.end(), 0.);
		if (densities.size() > 0)
			dens(i) /= densities.size();
		else
			zero_dens_pts.push_back(*gp);
	}

	for (gp = zero_dens_pts.begin(); gp != zero_dens_pts.end(); ++gp) {
			int N_d = 0;
			double dd = 0., d[8];
			// neighbourhood indices
			// closest neighbours
			int idx1[4] = { gp->m_Index_unterer_Nachbar,
				gp->m_Index_oberer_Nachbar, gp->m_Index_Nord_Nachbar,
				gp->m_Index_Sued_Nachbar };
			// diagonal neighbours
			int idx2[4] = { gp->m_Index_unterer_Nord_Nachbar,
				gp->m_Index_unterer_Sued_Nachbar,
				gp->m_Index_oberer_Nord_Nachbar, gp->m_Index_oberer_Sued_Nachbar};

			// save neighbourhood densities...
			for (int j = 0; j < 4; j++) {
				// closest neighbours
				if (idx1[j] != -1)
					d[j] = dens(idx1[j]);
				else
					d[j] = 0.;
				// weight diagonal neighbours with 1/sqrt(2)
				if (idx2[j] != -1)
					d[j + 4] = M_SQRT1_2 * dens(idx2[j]);
				else
					d[j + 4] = 0.;
			}

			// ...and average them if they are larger than zero
			for (int j = 0; j < 8; j++)
				if (d[j] > 0.) {
					dd += d[j];
					N_d++;
				}
			if (N_d > 0)
				dens(gp->m_eigener_Index) = dd / N_d;
			else
				// default if no other data is available
				// sensible for low altitudes
				dens(gp->m_eigener_Index) = 2.e16;
	}
}

#ifdef HAVE_NOEM
/*
 * interface to the SNOEM model for the SNOE NO data as apriori input
 * for the retrieval.
 * nx, ny, nz: dimensions of the longitude, latitude, altitude arrays
 * glon, glat, zkm: arrays for the longitude, latitude, and altitudes
 * kp, f107: solar data input at the respective day
 * snoe_3d(): the actual model procedure
 * *snoe_no: pointer to the output array (of dim. nx, ny, nz)
 */
extern "C" {
	int __params_MOD_nx, __params_MOD_ny, __params_MOD_nz;
	float __params_MOD_kp, __params_MOD_f107;
	float *__dynam_MOD_glon, *__dynam_MOD_glat, *__dynam_MOD_zkm;
	void __snoe_MOD_snoe_3d(int *doy, float *snoe_no);
}

/*
 * runs the model for the retrieval grid, for each latitude we have one
 * longitude and the model returns the number densities for the requested
 * altitudes, found in zkm[].
 * the analysed limb data is used to determine the date for the model run.
 *
 * this code must be linked with
 *   -lNOEM -lnetcdf
 * libNOEM schould be compiled with gfortran for the correct symbol names
 * and the model expects to find a file "noem_eof.nc" in the "input/" subdir.
 */
void SNOE_apriori_NO(Retrievalgitter &grid, Ausgewertete_Messung_Limb &aml,
		MPL_Matrix &apriori, Konfiguration &Konf)
{
	// to get the day of the year
	int days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	if (aml.m_Jahr % 4 == 0 &&
			!(aml.m_Jahr % 100 == 0 && aml.m_Jahr % 400 != 0))
		days[1] = 29;

	int doy = 0;
	for (int i = 0; i < (aml.m_Monat - 1); i++) doy += days[i];
	doy += aml.m_Tag;

	// get solar data from the spidr input files
	__params_MOD_f107 = spidr_value_from_file(aml.m_Jahr, aml.m_Monat,
			aml.m_Tag, Konf.m_Pfad_f107_adj_index, 150, 1);
	__params_MOD_kp = spidr_value_from_file(aml.m_Jahr, aml.m_Monat,
			aml.m_Tag, Konf.m_Pfad_Kp_index, 3, 1);
	std::cout << "# snoe parameters: f10.7 = " << __params_MOD_f107
		<< ", kp = " << __params_MOD_kp << std::endl;

	__params_MOD_nx = 1;
	__params_MOD_ny = 1;
	__params_MOD_nz = grid.m_Anzahl_Hoehen;

	__dynam_MOD_glon = new float[__params_MOD_nx];
	__dynam_MOD_glat = new float[__params_MOD_ny];
	__dynam_MOD_zkm = new float[__params_MOD_nz];
	float *snoe_no = new float[__params_MOD_nz];

	// initialise altitude array
	for (int i = 0; i < __params_MOD_nz; i++)
		__dynam_MOD_zkm[i] =
			my_clamp(grid.m_Gitter[i * grid.m_Anzahl_Breiten].m_Hoehe,
				Konf.NO_apriori_bottom, Konf.NO_apriori_top);

	// run the model for each latitude
	for (int i = 0; i < grid.m_Anzahl_Breiten; i++) {
		Gitterpunkt gp = grid.m_Gitter[i];
		__dynam_MOD_glon[0] = std::fmod(gp.longitude, 360.);
		if (grid.m_Anzahl_Breiten > 1)
			__dynam_MOD_glat[0] = gp.m_Breite;
		else
			__dynam_MOD_glat[0] = aml.m_Latitude_TP;
		__snoe_MOD_snoe_3d(&doy, snoe_no);
		for (int j = 0; j < __params_MOD_nz; j++) {
			apriori(j * grid.m_Anzahl_Breiten + i) = snoe_no[j];
		}
	}

	/* "Automatic" a priori scaling:
	 * In the case of a negative scale factor, we simply overwrite it with
	 * F10.7 / 150. 150 is the mean F10.7 during the SNOE measurement time
	 * period used for NOEM. This should avoid the "fixed scale factor"
	 * problem. We set the scale factor here because we already read the F10.7
	 * value for the NOEM model run. */
	if (Konf.NO_apriori_scale == -1) {
		// simply overwrite the configured scale factor
		Konf.NO_apriori_scale = __params_MOD_f107 / 150.;
		std::cout << "# snoe a priori auto scale: "
			<< Konf.NO_apriori_scale << std::endl;
	}

	// cleanup
	delete[] snoe_no;
	delete[] __dynam_MOD_zkm;
	delete[] __dynam_MOD_glat;
	delete[] __dynam_MOD_glon;
}
#else /* HAVE_NOEM */
void SNOE_apriori_NO(Retrievalgitter &grid, Ausgewertete_Messung_Limb &aml,
		MPL_Matrix &apriori, Konfiguration &Konf)
{
	// not available
	std::cerr << "NOEM model not available." << std::endl;
}
#endif /* HAVE_NOEM */

/*
 * NO apriori derived from the regression model.
 */
void regression_apriori_NO(Retrievalgitter &grid, Ausgewertete_Messung_Limb &aml,
		MPL_Matrix &apriori, Konfiguration &Konf)
{
	std::vector<double> alt_vec;
	for (int i = 0; i < grid.m_Anzahl_Hoehen; i++)
		alt_vec.push_back(
			my_clamp(grid.m_Gitter[i * grid.m_Anzahl_Breiten].m_Hoehe,
				Konf.NO_apriori_bottom, Konf.NO_apriori_top));

	for (int j = 0; j < grid.m_Anzahl_Breiten; j++) {
		double lat = grid.m_Gitter[j].m_Breite;
		std::vector<double> NO_model
				= NO_regress_model_python(aml, Konf, alt_vec, lat);
		for (int i = 0; i < grid.m_Anzahl_Hoehen; i++)
			apriori(i * grid.m_Anzahl_Breiten + j) = NO_model.at(i);
	}

	/* "Automatic" a priori scaling:
	 * In the case of a negative scale factor, we simply overwrite it with
	 * F10.7 / 150. 150 is the mean F10.7 during the SNOE measurement time
	 * period used for NOEM. This should avoid the "fixed scale factor"
	 * problem. We read the F10.7 in the same way as we do for the NOEM
	 * model run. */
	if (Konf.NO_apriori_scale == -1) {
		// get solar data from the spidr input files
		double f107 = spidr_value_from_file(aml.m_Jahr, aml.m_Monat,
				aml.m_Tag, Konf.m_Pfad_f107_adj_index, 150, 1);
		// simply overwrite the configured scale factor
		Konf.NO_apriori_scale = f107 / 150.;
		std::cout << "# regression a priori auto scale: "
			<< Konf.NO_apriori_scale << std::endl;
	}
}

/* arbitrary transition using given function */
double transition_func(double (*trans01)(double), double z, Konfiguration &Konf)
{
	return Phi_func(trans01, Konf.NO_apriori_bottom, Konf.NO_apriori_top,
			Konf.NO_apriori_smoothness, z);
}

/*
 * Scales the apriori values by a constant factor given by
 * Konf.NO_apriori_scale and an (optional) transition function.
 * The transition function restricts the altitudes at which the
 * apriori values are taken into account for the retrieval.
 */
void scale_apriori(Retrievalgitter &grid, MPL_Matrix &apriori,
		Konfiguration &Konf)
{
	/*
	 * Chose the transition (step) function from zero to one here.
	 * The current choices are:
	 *
	 * hardstep: = 0, x < 1; = 1, x >= 1
	 * my_phi: = 0, x < 0, = 1, x > 1, (in C^infinity)
	 * smoothstep: = 0, x < 0, = 1, x > 1, (with C^1 endpoints)
	 * smootherstep: = 0, x < 0, = 1, x > 1, (with C^2 endpoints)
	 *
	 * The default is my_phi.
	 */
	auto trans01 = my_phi;
	auto transition = std::bind(transition_func, trans01,
			std::placeholders::_1, std::placeholders::_2);
	double NO_apriori_scale = Konf.NO_apriori_scale;
	/* Calculate only the transition (i.e. don't scale)
	 * if the scale is negative. */
	if (NO_apriori_scale < 0)
		NO_apriori_scale = 1.;

	for (int i = 0; i < grid.m_Anzahl_Hoehen; i++) {
		double z = grid.m_Gitter[i * grid.m_Anzahl_Breiten].m_Hoehe;
		double f = transition(z, Konf);
		for (int j = 0; j < grid.m_Anzahl_Breiten; j++)
			apriori(i * grid.m_Anzahl_Breiten + j) *= f * NO_apriori_scale;
	}
}
