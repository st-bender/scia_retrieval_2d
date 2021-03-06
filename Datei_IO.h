/*
 * Datei_IO.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 19.04.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
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
		class Messung_Limb &Troposphaerische_Saeule, class Messung_Limb &mean_10_20,
		int Anzahl_Hoehen = -1, double dark_bg = 3.9e9);
std::vector<Messung_Limb> ReadL1C_Limb_meso_thermo_mpl_binary(std::string Dateiname,
		class Messung_Limb &niedrigste_Hoehe, class Messung_Limb &space);
std::vector<class Messung_Limb>
ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(std::string Dateiname,
		class Messung_Limb &niedrigste_Hoehe, class Messung_Limb &space,
		int Anzahl_Hoehen = -1, double dark_bg = 3.9e9);
std::vector<class Messung_Nadir> ReadL1C_Nadir_mpl_binary(std::string Dateiname, int &Anzahl_Messungen);
//Besser ein dynamisches Array einlesen, schneller als Vektor
////////////////////////////////////////////////////////////////////////////////
void Ausgabe_Saeulendichten(std::string Dateiname,
		std::vector<class Ausgewertete_Messung_Limb> &A_Messung_L);
void Ausgabe_Saeulendichten_back(std::string Dateiname,
		std::vector<class Ausgewertete_Messung_Limb> &aml_vec, MPL_Matrix &y);
void Ausgabe_Saeulendichten(std::string Dateiname,
		std::vector<class Ausgewertete_Messung_Nadir> &A_Messung_N);
//funktion ist ja überladbar

MPL_Matrix Read_Atmodatei(std::string Dateiname);
//Ausgeben
void Ausgabe_Dichten(std::string Dateiname_out, class Retrievalgitter &Grid,
		MPL_Matrix &Dichten, MPL_Matrix &Dichten_tot, MPL_Matrix &apriori,
		MPL_Matrix &S_x, MPL_Matrix &S_x_meas, MPL_Matrix &AKM,
		bool save_sx = true, bool save_akm = true);

#endif /* DATEI_IO_HH_ */
