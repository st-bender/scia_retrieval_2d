/*
 * Winkelstatistik.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 11.10.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 *
 *      Die 3 Teilabschnitte, in Luftmassenfaktoren_aufbauen in
 *      Matrizen_Aufbauen.cpp, in denen eine Statistik über die
 *      Winkelabweichungen zwischen Tangentialvektor und Zenitvektor von 90°
 *      erfasst werden wird, damit es dort übersichtlicher
 *      aussieht hier implementiert..... Die Winkelstatistik könnte schon
 *      früher durchgeführt werden(vermutlich noch vor den
 *      Säulendichten...Funktionen hier erleichern die spätere Migration
 */

#ifndef WINKELSTATISTIK_HH_
#define WINKELSTATISTIK_HH_

#include"MPL_Vektor.h"
#include<string>

class Winkelstatistik
{
public:
	//Konstruktor
	Winkelstatistik();
	//Methoden
	int Winkel_berechnen_und_einordnen(MPL_Vektor Verbindungsvektor,
									   MPL_Vektor Tangentenpunkt,
									   int &Winkel_OK);
	void Statistik_auf_Bildschirm_ausgeben();
	int Statistik_in_Datei_ausgeben(std::string Dateiname);

	int m_counter_0_001;  //winkel zwischen 0 und 0,001Grad
	int m_counter_0_002;
	int m_counter_0_005;
	int m_counter_0_01;
	int m_counter_0_02;    //Alles bis hier ist OK
	int m_counter_0_05;    //Das geht noch
	int m_counter_0_1;     // Ab hier wird so ziemlich alles falsch
	int m_counter_0_2;
	int m_counter_0_5;
	int m_counter_1;
	int m_counter_2;
	int m_counter_5;
	int m_counter_mehr;
};



#endif /* WINKELSTATISTIK_HH_ */
