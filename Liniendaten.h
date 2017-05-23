/*
 * Liniendaten.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 13.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
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
