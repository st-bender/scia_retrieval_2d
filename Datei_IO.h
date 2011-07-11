/*
 * Datei_IO.h
 *
 *  Created on: 19.04.2010
 *      Author: martin
 */
#ifndef DATEI_IO_HH_
#define DATEI_IO_HH_

#include <string>
#include <vector>
#include "MPL_Matrix.h"

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
std::vector<class Messung_Limb> ReadL1C_Limb_mpl_binary(std::string Dateiname,
		class Messung_Limb &Troposphaerische_Saeule, class Messung_Limb &mean_10_20);
std::vector<Messung_Limb> ReadL1C_Limb_meso_thermo_mpl_binary(std::string Dateiname,
		class Messung_Limb &niedrigste_Hoehe, class Messung_Limb &space);
std::vector<class Messung_Limb>
ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(std::string Dateiname,
		class Messung_Limb &niedrigste_Hoehe, class Messung_Limb &space, int Anzahl_Hoehen);
std::vector<class Messung_Nadir> ReadL1C_Nadir_mpl_binary(std::string Dateiname, int &Anzahl_Messungen);
//Besser ein dynamisches Array einlesen, schneller als Vektor
////////////////////////////////////////////////////////////////////////////////
int Ausgabe_Saeulendichten(std::string Dateiname,
		std::vector<class Ausgewertete_Messung_Limb> &A_Messung_L);
int Ausgabe_Saeulendichten_back(std::string Dateiname,
		std::vector<class Ausgewertete_Messung_Limb> &aml_vec, MPL_Matrix &y);
int Ausgabe_Saeulendichten(std::string Dateiname,
		std::vector<class Ausgewertete_Messung_Nadir> &A_Messung_N);
//funktion ist ja überladbar

MPL_Matrix Read_Atmodatei(std::string Dateiname);
//Ausgeben
int Ausgabe_Dichten(std::string Dateiname_out, class Retrievalgitter &Grid,
		MPL_Matrix &Dichten, MPL_Matrix &Dichten_tot,
		MPL_Matrix &S_x, MPL_Matrix &S_x_meas, MPL_Matrix &AKM);

#endif /* DATEI_IO_HH_ */
