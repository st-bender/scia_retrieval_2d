/*
 * Dateinamensteile_Bestimmen.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 23.09.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include<string>

std::string sb_basename(std::string filename)
{
	std::string bname = filename;
	std::string::size_type pos = filename.find_last_of("/\\");

	if (pos != std::string::npos)
		bname = filename.substr(pos + 1);

	return bname;
}
////////////////////////////////////////////////////
// Funktionsstart  xxxxx_Bestimmen
////////////////////////////////////////////////////
std::string xxxxx_Bestimmen(std::string Orbitlistenpfad)
{
	std::string bname = sb_basename(Orbitlistenpfad);
	std::string::size_type pos = bname.find_first_of("0123456789");

	if (pos != std::string::npos)
		return bname.substr(pos, 5);
	else
		return "xxxxx";
}
////////////////////////////////////////////////////
// ENDE xxxxx_Bestimmen
////////////////////////////////////////////////////

////////////////////////////////////////////////////
// Funktionsstart yyyymmdd_hhmm
////////////////////////////////////////////////////
std::string yyyymmdd_hhmm_Bestimmen(std::string Name_erste_Limbdatei)
{
	// Die Datei ist im MPL_binary Format
	// z.b. SCIA_limb_20040111_084344_1_0_09752.dat.l_mpl_binary
	std::string bname = sb_basename(Name_erste_Limbdatei);
	std::string::size_type pos = bname.find_first_of("0123456789");

	//TODO Das eventuell ganz umbennnen...
	// suchen eigentlich nur Datum...Uhrzeit Wuascht
	if (pos != std::string::npos)
		return bname.substr(pos, 8);
	else
		return "yyyymmdd";
}
////////////////////////////////////////////////////
// ENDE yyyymmdd_hhmm
////////////////////////////////////////////////////
