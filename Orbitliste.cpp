/*
 * Orbitliste.cpp
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

#include<fstream>
#include<string>
#include<iostream>
#include <cstdlib>
#include<vector>
#include <stdexcept>

#include"Orbitliste.h"

using std::cerr;
using std::endl;
using std::string;

//*****************************************************************
int Orbitliste::Liste_Laden(string Dateiname)
{
	std::ifstream infile(Dateiname.c_str());
	if (!(infile.is_open())) {
		std::cout << Dateiname << " kann nicht gefunden werden\n";
		return -1;
	}
	while (!(infile.eof())) {
		string Zeile;
		getline(infile, Zeile);
		m_Dateinamen.push_back(Zeile);
	}
	return 0;
} //ende Orbitliste::Liste_Laden
//*****************************************************************
//*****************************************************************
bool Orbitliste::Ist_Messung_Limbmessung(uint Index)
{
	//zuerstmal prüfen ob m_Dateinamen[Index] überhaupt existiert
	try {
		string tmp = m_Dateinamen.at(Index);
	} catch(std::out_of_range &err) {
		cerr << "Achtung Zugriff auf nicht vorhandenes Listenelement in "
			 << "Funktion Orbitliste::Ist_Messung_Limbmessung\n";
		cerr << "Exception: " << err.what() << endl;
		return false;
	}
	// charakteristische Zeichenkette suchen
	size_t Stringindex = m_Dateinamen[Index].find("SCIA_limb");
	// Überprüfen, ob Zeichenkette überhaupt vorhanden
	if (Stringindex == std::string::npos) {
		return false; //nicht gefunden
	}
	return true; //gefunden
}//ende Ist_Messung_Limbmessung
//*****************************************************************
//*****************************************************************
bool Orbitliste::Ist_Messung_Nadirmessung(uint Index)
{
	//zuerstmal prüfen ob m_Dateinamen[Index] überhaupt existiert
	try {
		string tmp = m_Dateinamen.at(Index);
	} catch(std::out_of_range &err) {
		cerr << "Achtung Zugriff auf nicht vorhandenes Listenelement in "
			 << "Funktion Orbitliste::Ist_Messung_Nadirmessung\n";
		cerr << "Exception: " << err.what() << endl;
		return false;
	}
	// charakteristische Zeichenkette suchen
	size_t Stringindex = m_Dateinamen[Index].find("SCIA_nadir");
	// Überprüfen, ob Zeichenkette überhaupt vorhanden
	if (Stringindex == std::string::npos) {
		return false;//nicht gefunden
	}
	return true;  //gefunden
}//ende Ist_Messung_Nadirmessung
//*****************************************************************

//Wartungsfunktion ///////////////////////////////////////////////////////
void Orbitliste::In_Datei_Speichern(string Dateiname)
{
	std::ofstream outfile(Dateiname.c_str());
	//ich verzichte hier mal auf weitere sicherheitsfragen
	for (unsigned int i = 0; i < this->m_Dateinamen.size(); i++) {
		outfile << m_Dateinamen[i] << "\n";
	}
}//Ende in_Datei_Speichern
////////////////////////////////////////////////////////////////////////////////
