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
	int lang = Orbitlistenpfad.size();
	//string xxxxx=Orbitlistenpfad.substr(lang-20,5);  old
	string xxxxx = Orbitlistenpfad.substr(lang - 15, 5);
	// cout<<xxxxx<<"\n";
	return xxxxx;
}
////////////////////////////////////////////////////
// ENDE xxxxx_Bestimmen
////////////////////////////////////////////////////

////////////////////////////////////////////////////
// Funktionsstart yyyymmdd_hhmm
////////////////////////////////////////////////////
string yyyymmdd_hhmm_Bestimmen(string Name_erste_Limbdatei)
{
	string yyyymmdd_hhmm;
	// Die Datei ist im MPL_binary Format
	// z.b. SCIA_limb_20040111_084344_1_0_09752.dat.l_mpl_binary
	int lang = Name_erste_Limbdatei.size();
	yyyymmdd_hhmm = Name_erste_Limbdatei.substr(lang - 42, 13);
	//TODO Das eventuell ganz umbennnen...
	// suchen eigentlich nur Datum...Uhrzeit Wuascht
	yyyymmdd_hhmm = Name_erste_Limbdatei.substr(lang - 42, 8);
	return yyyymmdd_hhmm;
}
////////////////////////////////////////////////////
// ENDE yyyymmdd_hhmm
////////////////////////////////////////////////////
