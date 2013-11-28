/*
 * Orbitliste.cpp
 *
 *  Created on: 15.04.2010
 *      Author: martin
 */

#include<fstream>
#include<string>
#include<iostream>
#include <cstdlib>
#include<vector>
#include <stdexcept>

#include"Orbitliste.h"

using namespace std;

//*****************************************************************
int Orbitliste::Liste_Laden(string Dateiname)
{
	ifstream infile(Dateiname.c_str());
	if (!(infile.is_open())) {
		cout << Dateiname << " kann nicht gefunden werden\n";
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
	ofstream outfile(Dateiname.c_str());
	//ich verzichte hier mal auf weitere sicherheitsfragen
	for (unsigned int i = 0; i < this->m_Dateinamen.size(); i++) {
		outfile << m_Dateinamen[i] << "\n";
	}
}//Ende in_Datei_Speichern
////////////////////////////////////////////////////////////////////////////////
