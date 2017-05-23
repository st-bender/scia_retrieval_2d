/*
 * Orbitliste.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 15.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
// In dieser Liste stehen alle Dateinamen, die zu einem Orbit gehören
// Diese wird aus einer Datei geladen -> Die Orbitliste wird dort erzeugt
#ifndef ORBITLISTE_HH_
#define ORBITLISTE_HH_

#include<fstream>
#include<string>
#include<iostream>
#include <cstdlib>
#include<vector>

class Orbitliste
{
public:
	int Liste_Laden(std::string Dateiname);
	bool Ist_Messung_Limbmessung(uint Index);
	bool Ist_Messung_Nadirmessung(uint Index);
	std::vector<std::string> m_Dateinamen;

	//Wartungsfunktion
	void In_Datei_Speichern(std::string Dateiname);
};


#endif /* ORBITLISTE_HH_ */
