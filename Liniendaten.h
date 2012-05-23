/*
 * Liniendaten.h
 *
 *  Created on: 13.04.2010
 *      Author: martin
 */

#ifndef  LINIENDATEN_HH_
#define LINIENDATEN_HH_

#include <string>

class Liniendaten
{
	//Daten für eine Linie //Da nur wenig Daten-> nicht zeitkritisch
public:
	//Konstruktor
	Liniendaten();
	//Methoden
	void Einlesen(std::string Dateiname, double Wellenlaenge);  // aus Datei einlesen
	void Emissivitaet_ermitteln(); // aus Rohdaten Gamma berechnen
								   // !! ACHTUNG HIER GIBTS DISKUSSIONSBEDARF!!
								   // siehe Code
	void Auf_Bildschirm_Ausgeben();


	//Membervariablen
	// eingelesene
	double m_Wellenlaenge;  //in nm
	double m_rel_Einstein;
	double m_f_Wert;            // Oszillatorstärke
	double m_E1;
	double m_E2;
	//
	//double m_theta;               //Streuwinkel Theta ist Lichtwegabhängig, also nicht konstant
	//errechnete
	double m_Gamma;          //Emissivität
};


#endif /* LINIENDATEN_HH_ */
