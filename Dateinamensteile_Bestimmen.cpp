/*
 * Dateinamensteile_Bestimmen.cpp
 *
 *  Created on: 23.09.2010
 *      Author: martin
 */

#include<string>
#include"Dateinamensteile_Bestimmen.h"

using namespace std;

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
string xxxxx_Bestimmen(string Orbitlistenpfad)
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
string yyyymmdd_hhmm_Bestimmen(string Name_erste_Limbdatei)
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
