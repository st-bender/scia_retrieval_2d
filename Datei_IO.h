/*
 * Datei_IO.h
 *
 *  Created on: 19.04.2010
 *      Author: martin
 */

#include<string>
#include<vector>
#include"Messung_Limb.h"
#include"Messung_Nadir.h"
#include"Ausgewertete_Messung_Limb.h"
#include"Ausgewertete_Messung_Nadir.h"
#include"MPL_Matrix.h"
#include"Retrievalgitter.h"

#ifndef DATEI_IO_HH_
#define DATEI_IO_HH_
// Einlesen
// Alte Ascii einlesefunktionen......
// später löschen... sind im Feuerballprojekt viel besser
//vector<Messung_Limb> ReadL1C_Limb(string Dateiname);
//Messung_Nadir* ReadL1C_Nadir(string Dateiname, int& Anzahl_Messungen);
// //Besser ein dynamisches Array einlesen, schneller als Vektor
// EINLESEN DER BINÄREN DATEN/////////////////////////////////////////////////
// Beide Funktionen sehen aus entwicklungsgeschichtlichen Gründen
// unterschiedlich aus (vektor vs array) zeittechnisch ist die nadirvariante
// günstiger...kann später mal gefixt werden...erstmal programm fertig kriegen
vector<Messung_Limb> ReadL1C_Limb_mpl_binary(string Dateiname,
		Messung_Limb &Troposphaerische_Saeule, Messung_Limb &mean_10_20);
vector<Messung_Limb> ReadL1C_Limb_meso_thermo_mpl_binary(string Dateiname,
		Messung_Limb &niedrigste_Hoehe);
vector<Messung_Limb>
ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(string Dateiname,
		Messung_Limb &niedrigste_Hoehe, int Anzahl_Hoehen);
vector<Messung_Nadir> ReadL1C_Nadir_mpl_binary(string Dateiname, int &Anzahl_Messungen);
//Besser ein dynamisches Array einlesen, schneller als Vektor
////////////////////////////////////////////////////////////////////////////////
int Ausgabe_Saeulendichten(string Dateiname,
		vector<Ausgewertete_Messung_Limb> &A_Messung_L);
int Ausgabe_Saeulendichten(string Dateiname,
		vector<Ausgewertete_Messung_Nadir> &A_Messung_N);
//funktion ist ja überladbar

MPL_Matrix Read_Atmodatei(string Dateiname);
//Ausgeben
int Ausgabe_Dichten(string Dateiname_out, Retrievalgitter &Grid,
		MPL_Matrix &Dichten, MPL_Matrix &S_x, MPL_Matrix &AKM);

#endif /* DATEI_IO_HH_ */
