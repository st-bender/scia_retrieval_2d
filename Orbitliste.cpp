/*
 * Orbitliste.cpp
 *
 *  Created on: 15.04.2010
 *      Author: martin
 */

#include<fstream>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<vector>

#include"Orbitliste.h"

using namespace std;

//*****************************************************************
int Orbitliste::Liste_Laden(string Dateiname)
{
	ifstream infile;
	infile.open(Dateiname.c_str());
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
int Orbitliste::Ist_Messung_Limbmessung(uint Index)    //0 ja alles andere nein
{
	//zuerstmal prüfen ob m_Dateinamen[Index] überhaupt existiert
	if (Index >= m_Dateinamen.size()) {
		cout << "Achtung Zugriff auf nicht vorhandenes Listenelement in Funktion Orbitliste::Ist_Messung_Limbmessung\n";
		return 1;
	}
	// charakteristische Zeichenkette suchen
	int Stringindex = m_Dateinamen[Index].find("SCIA_limb");
	// Überprüfen, ob Zeichenkette überhaupt vorhanden
	if (Stringindex == -1) {
		return 2; //nichtr gefunden
	}
	return 0; //gefunden
}//ende Ist_Messung_Limbmessung
//*****************************************************************
//*****************************************************************
int Orbitliste::Ist_Messung_Nadirmessung(uint Index)   //0 ja alles andere nein
{
	//zuerstmal prüfen ob m_Dateinamen[Index] überhaupt existiert
	if (Index >= m_Dateinamen.size()) {
		cout << "Achtung Zugriff auf nicht vorhandenes Listenelement in Funktion Orbitliste::Ist_Messung_Nadirmessung\n";
		return 1;
	}
	// charakteristische Zeichenkette suchen
	int Stringindex = m_Dateinamen[Index].find("SCIA_nadir");
	// Überprüfen, ob Zeichenkette überhaupt vorhanden
	if (Stringindex == -1) {
		return 2;//nicht gefunden
	}
	return 0;  //gefunden
}//ende Ist_Messung_Nadirmessung
//*****************************************************************

//Wartungsfunktion ///////////////////////////////////////////////////////
void Orbitliste::In_Datei_Speichern(string Dateiname)
{
	ofstream outfile;
	outfile.open(Dateiname.c_str());
	//ich verzichte hier mal auf weitere sicherheitsfragen
	for (unsigned int i = 0; i < this->m_Dateinamen.size(); i++) {
		outfile << m_Dateinamen[i] << "\n";
	}
	outfile.close();
}//Ende in_Datei_Speichern
///////////////////////////////////////////////////////////////////////////////////
