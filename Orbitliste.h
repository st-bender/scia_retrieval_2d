/*
 * Orbitliste.h
 *
 *  Created on: 15.04.2010
 *      Author: martin
 */
// In dieser Liste stehen alle Dateinamen, die zu einem Orbit gehören
// Diese wird aus einer Datei geladen -> Die Orbitliste wird dort erzeugt


#include<fstream>
#include<string>
#include<iostream>
#include<stdlib.h>
#include<vector>


using namespace std;



#ifndef ORBITLISTE_HH_
#define ORBITLISTE_HH_


class Orbitliste
{
public:
	int Liste_Laden(string Dateiname);
	int Ist_Messung_Limbmessung(uint Index);   //0 ja alles andere nein
	int Ist_Messung_Nadirmessung(uint Index);  //0 ja alles andere nein
	vector<string> m_Dateinamen;

	//Wartungsfunktion
	void In_Datei_Speichern(string Dateiname);
};


#endif /* ORBITLISTE_HH_ */
