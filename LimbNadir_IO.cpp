/*
 * LimbNadir_IO.cpp
 *
 *  Created on: 10.08.2010
 *      Author: martin
 */

#include"LimbNadir_IO.h"

#include<string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <netcdf.h>

using std::cout;
using std::ios;
using std::ios_base;
using std::ifstream;
using std::ofstream;
using std::string;

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                          LADEN                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


template<typename T>
std::istream& binary_read(std::istream* stream, T& value, size_t N = 1) {
	return stream->read(reinterpret_cast<char *>(&value), N * sizeof(T));
}

inline void Limb_Datensatz::read_from_mpl_binary(std::ifstream *stream,
		int no_of_pix)
{
	m_N_radiances = no_of_pix;
	m_radiance.resize(no_of_pix);
	m_error.resize(no_of_pix);
	binary_read(stream, m_Sub_Sat_Lat);
	binary_read(stream, m_Sub_Sat_Lon);
	binary_read(stream, m_TP_Lat);
	binary_read(stream, m_TP_Lon);
	binary_read(stream, m_Tangentenhoehe);
	binary_read(stream, m_TP_SZA);
	binary_read(stream, m_TP_SAA);
	binary_read(stream, m_TP_LOS_Zenit);
	binary_read(stream, m_TOA_SZA);
	binary_read(stream, m_TOA_SAA);
	binary_read(stream, m_TOA_LOS_Zenit);
	binary_read(stream, m_Sat_SZA);
	binary_read(stream, m_Sat_SAA);
	binary_read(stream, m_Sat_LOS_Zenit);
	binary_read(stream, m_Sat_Hoehe);
	binary_read(stream, m_Erdradius);
	//cout<<"Lese feld\n";
	binary_read(stream, m_radiance[0], no_of_pix);
	binary_read(stream, m_error[0], no_of_pix);
}

////////////////////////////////////////////////////////////////////////////////
//int Load_Limb_Ascii Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Load_Limb_Ascii(string Datei_in, string textheader[31], int &no_of_alt,
					int &no_of_pix, int Orbitstate[5], int Datum[6],
					float Center_Lat_Lon[10], float &orbit_phase,
					float*& Wellenlaengen, Limb_Datensatz*& Limbdaten)
{
	int lang_textheader = 31;
	////////////////////////////////////////////////////////////////////////////
	//
	//Schritt 1 Laden der Datei in Datenstrukturen
	//
	////////////////////////////////////////////////////////////////////////////
	ifstream infile(Datei_in.c_str());
	if (!infile.is_open()) {
		cout << Datei_in << " konnte nicht zum lesen geöffnet werden\n";
		return 1;
	}
	infile >> lang_textheader; // number of comment lines
	// Datenstrukturen einlesen
	for (int i = 0; i <= lang_textheader; i++) {
		getline(infile, textheader[i]);
	}
	infile >> no_of_alt;
	infile >> no_of_pix;
	for (int i = 0; i < 5; i++) {
		infile >> Orbitstate[i];
	}
	for (int i = 0; i < 6; i++) {
		infile >> Datum[i];
	}
	if (no_of_alt != 31) {
		cout << "no_of_alt:" << no_of_alt
			 << "\n no_of_alt ist nicht 31...Fehler\n";
		return 2;
	}
	Limbdaten = new Limb_Datensatz[no_of_alt];
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Sub_Sat_Lat;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Sub_Sat_Lon;
	}
	infile >> orbit_phase;
	for (int i = 0; i < 10; i++) {
		infile >> Center_Lat_Lon[i];
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TP_Lat;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TP_Lon;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Tangentenhoehe;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TP_SZA;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TP_SAA;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TP_LOS_Zenit;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TOA_SZA;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TOA_SAA;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_TOA_LOS_Zenit;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Sat_SZA;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Sat_SAA;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Sat_LOS_Zenit;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Sat_Hoehe;
	}
	for (int i = 0; i < no_of_alt; i++) {
		infile >> Limbdaten[i].m_Erdradius;
	}
	// Großes Feld einlesen
	//!!!!!!!!!!!!!!!!! Gehe von 31 Höhen aus!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	int Anzahl_Big_Array = (no_of_alt + 1) * no_of_pix;
	float *Big_Array;
	Big_Array = new float[Anzahl_Big_Array];
	for (int i = 0; i < Anzahl_Big_Array; i++) {
		infile >> Big_Array[i];
	}
	string errortext;
	infile >> errortext;
	if (errortext != "ERRORS") {
		// Fehlermeldung und Funktion verlassen
		cout << "errortext: " << errortext << "\n";
		cout << "Da sollte ERRORS stehen, sonst ist was falsch\n";
		delete[] Big_Array;
		return 3;
	}
	float *Big_Array2;  // Feld mit Fehlern
	Big_Array2 = new float[Anzahl_Big_Array];
	for (int i = 0; i < Anzahl_Big_Array; i++) {
		infile >> Big_Array2[i];
	}

	// einlesen beendet
	infile.close();
	////////////////////////////////////////////////////////////////////////////
	//
	//Schritt 2 Speichern des Arrays in Spaltenarrays
	//
	////////////////////////////////////////////////////////////////////////////
	// Das ist zwar länglich und redundandant, aber vermutlich schneller als
	// immer das ganze Feld zu betrachten und hier wird das nur einmal gemacht,
	// sodass es nicht bei jedem laden neu gemacht werden muss
	Wellenlaengen = new float[no_of_pix];
	for (int i = 0; i < no_of_alt; i++) {
		Limbdaten[i].m_N_radiances = no_of_pix;
		Limbdaten[i].m_radiance.resize(no_of_pix);
		Limbdaten[i].m_error.resize(no_of_pix);
	}
	for (int i = 0; i < no_of_pix; i++) {
		Wellenlaengen[i] = Big_Array[i * (no_of_alt + 1)];
		for (int j = 0; j < no_of_alt; j++) {
			Limbdaten[j].m_radiance[i] = Big_Array[i * (no_of_alt + 1) + j + 1];
			Limbdaten[j].m_error[i] = Big_Array2[i * (no_of_alt + 1) + j + 1];
		}
	}
	delete[] Big_Array;
	delete[] Big_Array2;
	//cout<<Wellenlaengen[30]<<"\n";
	//cout<<"Laden fertig\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Load_Limb_Ascii
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Load_Limb_l_mpl_binary
////////////////////////////////////////////////////////////////////////////////
int Load_Limb_l_mpl_binary(string Datei_in,
						   string textheader[31], int &no_of_alt,
						   int &no_of_pix, int Orbitstate[5], int Datum[6],
						   float Center_Lat_Lon[10], float &orbit_phase,
						   std::vector<float> &Wellenlaengen,
						   std::vector<Limb_Datensatz> &Limbdaten)
{
	int lang_textheader = 31;
	////////////////////////////////////////////////////////////////////////////
	//
	//Schritt 1 Laden der Datei in Datenstrukturen
	//
	////////////////////////////////////////////////////////////////////////////
	ifstream infile(Datei_in.c_str(), ios_base::binary);
	if (!infile.is_open()) {
		cout << Datei_in << " konnte nicht zum lesen geöffnet werden\n";
		return 1;
	}
	//cout<<"Lese header\n";
	std::stringstream ss;
	char c_textheader[100];
	infile.read((char *) c_textheader, sizeof(char) * 100);
	ss << c_textheader;
	ss >> lang_textheader;
	for (int i = 0; i < lang_textheader; i++) {
		int pos = infile.tellg();
		infile.read((char *) c_textheader, sizeof(char) * 100);
		if (c_textheader[0] != '#') {
			infile.seekg(pos);
			break;
		}
		textheader[i] = c_textheader;
		//cout<<textheader[i]<<"\n";
	}
	//cout<<"Lese header2\n";
	infile.read((char *) &no_of_alt, sizeof(int));
	infile.read((char *) &no_of_pix, sizeof(int));
	infile.read((char *) Orbitstate, sizeof(int) * 5);
	infile.read((char *) Datum, sizeof(int) * 6);
	infile.read((char *) Center_Lat_Lon, sizeof(float) * 10);
	// the orbit phase is only there if the textheader is longer than 29 lines
	if (lang_textheader > 29)
		infile.read((char *) &orbit_phase, sizeof(float));
	// no of pixrel bekannt....speicher reservieren
	Wellenlaengen.resize(no_of_pix);
	//cout<<"Lese Wellenlängen\n";
	infile.read((char *) Wellenlaengen.data(), sizeof(float)*no_of_pix);
	Limbdaten.reserve(no_of_alt);
	for (int i = 0; i < no_of_alt; i++) {
		//cout<<"Header2.2\n";
		Limb_Datensatz lds;
		lds.read_from_mpl_binary(&infile, no_of_pix);
		Limbdaten.push_back(lds);
	}
	//cout<<"Laden fertig\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Load_Limb_l_mpl_binary
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Load_Limb_l_nc
////////////////////////////////////////////////////////////////////////////////
int Load_Limb_l_nc(string Datei_in,
				   string textheader[31], int &no_of_alt,
				   int &no_of_pix, int Orbitstate[5], int Datum[6],
				   float Center_Lat_Lon[10], float &orbit_phase,
				   std::vector<float> &Wellenlaengen,
				   std::vector<Limb_Datensatz> &Limbdaten)
{
	int infile = ncopen(Datei_in.c_str(), NC_NOWRITE);
	if (!infile) {
		cout << Datei_in << " konnte nicht zum lesen geöffnet werden\n";
		return 1;
	}
	/* Read the textheader first.
	 * We save the netcdf attribute in a temporary char array first and
	 * then manually convert it to the requested array of strings. */
	int th_len;
	nc_type th_type;
	ncattinq(infile, NC_GLOBAL, "textheader", &th_type, &th_len);
	char* th_read = new char[th_len];
	ncattget(infile, NC_GLOBAL, "textheader", th_read);
	std::stringstream th_ss(th_read);
	string th_item;
	int th_cnt = 0;
	while (std::getline(th_ss, th_item, '\n')) {
		textheader[th_cnt] = th_item;
		++th_cnt;
	}
	delete[] th_read;
	/* Directly reading into orbit_phase does not seem to work properly.
	 * Thus, we also use a temporary variable for it. The rest of the
	 * netcdf reading seems to work as expected, though. */
	float ophase;
	ncattget(infile, NC_GLOBAL, "orbit_phase", &ophase);
	ncattget(infile, NC_GLOBAL, "orbit_state", Orbitstate);
	ncattget(infile, NC_GLOBAL, "date", Datum);
	ncattget(infile, NC_GLOBAL, "cent_lat_lon", Center_Lat_Lon);
	int lid = ncdimid(infile, "limb");
	int wlid = ncdimid(infile, "wavelength");
	ncdiminq(infile, lid, 0, (long *)&no_of_alt);
	ncdiminq(infile, wlid, 0, (long *)&no_of_pix);
	long start[] = {0};
	long acount[] = {no_of_alt};
	long wlcount[] = {no_of_pix};
	long rstart[] = {0, 0};
	long rcount[] = {no_of_alt, no_of_pix};
	Wellenlaengen.resize(no_of_pix);
	ncvarget(infile, ncvarid(infile, "wavelength"), start, wlcount, Wellenlaengen.data());
	std::vector<float> v_Sub_Sat_Lat(no_of_alt);
	std::vector<float> v_Sub_Sat_Lon(no_of_alt);
	std::vector<float> v_TP_Lat(no_of_alt);
	std::vector<float> v_TP_Lon(no_of_alt);
	std::vector<float> v_Tangentenhoehe(no_of_alt);
	std::vector<float> v_TP_SZA(no_of_alt);
	std::vector<float> v_TP_SAA(no_of_alt);
	std::vector<float> v_TP_LOS_Zenit(no_of_alt);
	std::vector<float> v_TOA_SZA(no_of_alt);
	std::vector<float> v_TOA_SAA(no_of_alt);
	std::vector<float> v_TOA_LOS_Zenit(no_of_alt);
	std::vector<float> v_Sat_SZA(no_of_alt);
	std::vector<float> v_Sat_SAA(no_of_alt);
	std::vector<float> v_Sat_LOS_Zenit(no_of_alt);
	std::vector<float> v_Sat_Hoehe(no_of_alt);
	std::vector<float> v_Erdradius(no_of_pix);
	ncvarget(infile, ncvarid(infile, "sub_sat_lat"), start, acount, v_Sub_Sat_Lat.data());
	ncvarget(infile, ncvarid(infile, "sub_sat_lon"), start, acount, v_Sub_Sat_Lon.data());
	ncvarget(infile, ncvarid(infile, "TP latitude"), start, acount, v_TP_Lat.data());
	ncvarget(infile, ncvarid(infile, "TP longitude"), start, acount, v_TP_Lon.data());
	ncvarget(infile, ncvarid(infile, "TP altitude"), start, acount, v_Tangentenhoehe.data());
	ncvarget(infile, ncvarid(infile, "TP SZA"), start, acount, v_TP_SZA.data());
	ncvarget(infile, ncvarid(infile, "TP SAA"), start, acount, v_TP_SAA.data());
	ncvarget(infile, ncvarid(infile, "TP LOS Zenit"), start, acount, v_TP_LOS_Zenit.data());
	ncvarget(infile, ncvarid(infile, "TOA SZA"), start, acount, v_TOA_SZA.data());
	ncvarget(infile, ncvarid(infile, "TOA SAA"), start, acount, v_TOA_SAA.data());
	ncvarget(infile, ncvarid(infile, "TOA LOS Zenit"), start, acount, v_TOA_LOS_Zenit.data());
	ncvarget(infile, ncvarid(infile, "SAT SZA"), start, acount, v_Sat_SZA.data());
	ncvarget(infile, ncvarid(infile, "SAT SAA"), start, acount, v_Sat_SAA.data());
	ncvarget(infile, ncvarid(infile, "SAT LOS Zenit"), start, acount, v_Sat_LOS_Zenit.data());
	ncvarget(infile, ncvarid(infile, "SAT altitude"), start, acount, v_Sat_Hoehe.data());
	ncvarget(infile, ncvarid(infile, "earthradius"), start, acount, v_Erdradius.data());
	std::vector<float> rads(no_of_alt * no_of_pix);
	std::vector<float> errs(no_of_alt * no_of_pix);
	ncvarget(infile, ncvarid(infile, "radiance"), rstart, rcount, rads.data());
	ncvarget(infile, ncvarid(infile, "radiance errors"), rstart, rcount, errs.data());
	for (int i = 0; i < no_of_alt; ++i) {
		Limb_Datensatz lds;
		lds.m_N_radiances = no_of_pix;
		lds.m_Sub_Sat_Lat = v_Sub_Sat_Lat.at(i);
		lds.m_Sub_Sat_Lon = v_Sub_Sat_Lon.at(i);
		lds.m_TP_Lat = v_TP_Lat.at(i);
		lds.m_TP_Lon = v_TP_Lon.at(i);
		lds.m_Tangentenhoehe = v_Tangentenhoehe.at(i);
		lds.m_TP_SZA = v_TP_SZA.at(i);
		lds.m_TP_SAA = v_TP_SAA.at(i);
		lds.m_TP_LOS_Zenit = v_TP_LOS_Zenit.at(i);
		lds.m_TOA_SZA = v_TOA_SZA.at(i);
		lds.m_TOA_SAA = v_TOA_SAA.at(i);
		lds.m_TOA_LOS_Zenit = v_TOA_LOS_Zenit.at(i);
		lds.m_Sat_SZA = v_Sat_SZA.at(i);
		lds.m_Sat_SAA = v_Sat_SAA.at(i);
		lds.m_Sat_LOS_Zenit = v_Sat_LOS_Zenit.at(i);
		lds.m_Sat_Hoehe = v_Sat_Hoehe.at(i);
		lds.m_Erdradius = v_Erdradius.at(i);
		lds.m_radiance.assign(rads.begin() + i * no_of_pix, rads.begin() + (i + 1) * no_of_pix);
		lds.m_error.assign(errs.begin() + i * no_of_pix, errs.begin() + (i + 1) * no_of_pix);
		Limbdaten.push_back(lds);
	}
	orbit_phase = ophase;

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Load_Limb_l_nc
////////////////////////////////////////////////////////////////////////////////

inline void Nadir_Datensatz::read_from_mpl_binary(std::ifstream *stream,
		int no_of_pix)
{
	m_N_radiances = no_of_pix;
	m_radiance.resize(no_of_pix);
	m_error.resize(no_of_pix);
	binary_read(stream, m_Messung_ID);
	binary_read(stream, m_state_ID);
	binary_read(stream, m_Jahr);
	binary_read(stream, m_Monat);
	binary_read(stream, m_Tag);
	binary_read(stream, m_Stunde);
	binary_read(stream, m_Minute);
	binary_read(stream, m_Sekunde);
	binary_read(stream, m_SZA_TOA[0], 3);
	binary_read(stream, m_SAA_TOA[0], 3);
	binary_read(stream, m_LOS_Zenit_Winkel[0], 3);
	binary_read(stream, m_LOS_Azimut_Winkel[0], 3);
	binary_read(stream, m_Hoehe);
	binary_read(stream, m_Sat_Lat);
	binary_read(stream, m_Sat_Lon);
	binary_read(stream, m_Sat_Erdradius);
	binary_read(stream, m_orbit_phase);
	binary_read(stream, m_geo_nadir_corner_lat[0], 4);
	binary_read(stream, m_geo_nadir_corner_lon[0], 4);
	binary_read(stream, m_geo_nadir_center_lat);
	binary_read(stream, m_geo_nadir_center_lon);
	binary_read(stream, m_Integrationszeit);
	binary_read(stream, m_radiance[0], no_of_pix);
	binary_read(stream, m_error[0], no_of_pix);
}

////////////////////////////////////////////////////////////////////////////////
//int Load_Nadir_Ascii();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Load_Nadir_Ascii(string Datei_in,
					 string textheader[7], int &No_of_Messungen, int &No_of_Pix,
					 int*& Kanal_Nr, float*& Wellenlaenge,
					 Nadir_Datensatz*& Nadirdaten)
{
	int Anzahl_Textheaderzeilen = 7;

	string s_dummy;
	int i_dummy;
	float f_dummy;

	//////////////////////////////////////////////////////////
	// Ascii Datei Laden
	/////////////////////////////////////////////////////////
	// ERSTER DURCHLAUF ARRAYDIMENSIONEN ERMITTELN
	ifstream infile(Datei_in.c_str());
	if (!infile.is_open()) {
		cout << Datei_in << " kann nicht geladen werden in Nadirumwandlung\n";
		return 1;
	}
	for (int i = 0; i < Anzahl_Textheaderzeilen; i++) {
		getline(infile, textheader[i]);
	}
	infile >> No_of_Messungen;
	getline(infile, s_dummy); //ende der Zeile
	//wir brauchen noch die Anzahl der Pixel...
	//die für alle Messungen gleich sind
	for (int i = 0; i < 8; i++) {
		getline(infile, s_dummy);
	}
	infile >> No_of_Pix;
	//cout<<"No_of_Pix: "<<No_of_Pix<<"\n";
	//speicher anmelden
	Kanal_Nr = new int[No_of_Pix];
	Wellenlaenge = new float[No_of_Pix];
	Nadirdaten = new Nadir_Datensatz[No_of_Messungen];
	// ZWEITER DURCHLAUF Kanal und WL Felder füllen
	infile.clear();
	infile.seekg(0, ios::beg);
	if (!infile.is_open()) {
		cout << Datei_in << " kann nicht geladen werden in Nadirumwandlung\n";
		return 1;
	}

	for (int i = 0; i < 17; i++) {
		getline(infile, s_dummy);
	}
	for (int i = 0; i < No_of_Pix; i++) {
		infile >> Kanal_Nr[i]
			   >> Wellenlaenge[i];
		getline(infile, s_dummy);
	}
	//DRITTER DURCHLAUF großes Feld Füllen
	infile.clear();
	infile.seekg(0, ios::beg);
	if (!infile.is_open()) {
		cout << Datei_in << " kann nicht geladen werden in Nadirumwandlung\n";
		return 1;
	}
	for (int i = 0; i < 8; i++) {
		getline(infile, s_dummy);   //header von datei
	}
	for (int i = 0; i < No_of_Messungen; i++) {
		Nadirdaten[i].m_N_radiances = No_of_Pix;
		Nadirdaten[i].m_radiance.resize(No_of_Pix);
		Nadirdaten[i].m_error.resize(No_of_Pix);
		//Zeile mit Messungsnummer...von 0 an durchnummeriert
		infile >> s_dummy >> Nadirdaten[i].m_Messung_ID;
		infile >> Nadirdaten[i].m_state_ID;
		infile >> Nadirdaten[i].m_Jahr >> Nadirdaten[i].m_Monat
			   >> Nadirdaten[i].m_Tag >> Nadirdaten[i].m_Stunde
			   >> Nadirdaten[i].m_Minute >> Nadirdaten[i].m_Sekunde;
		infile >> Nadirdaten[i].m_SZA_TOA[0] >> Nadirdaten[i].m_SAA_TOA[0]
			   >> Nadirdaten[i].m_SZA_TOA[1] >> Nadirdaten[i].m_SAA_TOA[1]
			   >> Nadirdaten[i].m_SZA_TOA[2] >> Nadirdaten[i].m_SAA_TOA[2];
		infile >> Nadirdaten[i].m_LOS_Zenit_Winkel[0] >> Nadirdaten[i].m_LOS_Azimut_Winkel[0]
			   >> Nadirdaten[i].m_LOS_Zenit_Winkel[1] >> Nadirdaten[i].m_LOS_Azimut_Winkel[1]
			   >> Nadirdaten[i].m_LOS_Zenit_Winkel[2] >> Nadirdaten[i].m_LOS_Azimut_Winkel[2];
		infile >> Nadirdaten[i].m_Hoehe >> Nadirdaten[i].m_Sat_Lat
			   >> Nadirdaten[i].m_Sat_Lon >> Nadirdaten[i].m_Sat_Erdradius
			   >> Nadirdaten[i].m_orbit_phase;
		infile >> Nadirdaten[i].m_geo_nadir_corner_lat[0] >> Nadirdaten[i].m_geo_nadir_corner_lon[0]
			   >> Nadirdaten[i].m_geo_nadir_corner_lat[1] >> Nadirdaten[i].m_geo_nadir_corner_lon[1]
			   >> Nadirdaten[i].m_geo_nadir_corner_lat[2] >> Nadirdaten[i].m_geo_nadir_corner_lon[2]
			   >> Nadirdaten[i].m_geo_nadir_corner_lat[3] >> Nadirdaten[i].m_geo_nadir_corner_lon[3]
			   >> Nadirdaten[i].m_geo_nadir_center_lat >> Nadirdaten[i].m_geo_nadir_center_lon;
		infile >> Nadirdaten[i].m_Integrationszeit;
		infile >> i_dummy;
		for (int j = 0; j < No_of_Pix; j++) {
			infile >> f_dummy >> f_dummy >> Nadirdaten[i].m_radiance[j];
		}
		string errortext;
		infile >> errortext;
		if (errortext != "ERRORS") {
			cout << "errortext: " << errortext << "\n";
			cout << " Fehler in Load_Nadir_Ascii\n"
				 << "errortext ist nicht ERRORS... "
				 << "Datei in undefinierten Zustand\n";
			return 1;
		}
		// Fehlerfeld einlesen
		for (int j = 0; j < No_of_Pix; j++) {
			infile >> f_dummy >> f_dummy >> Nadirdaten[i].m_error[j];
		}

	}// ende for i
	//cout<<"ASCII Nadir Laden fertig\n";
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//ENDE int Load_Nadir_Ascii();
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//int Load_Nadir_n_mpl_binary();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Load_Nadir_n_mpl_binary(string Datei_in,
							string textheader[7], int &No_of_Messungen,
							int &No_of_Pix, std::vector<int> &Kanal_Nr,
							std::vector<float> &Wellenlaenge,
							std::vector<Nadir_Datensatz> &Nadirdaten)
{
	int Anzahl_Textheaderzeilen = 7;
	////////////////////////////////////////////////////////////////
	// Binäre Datei Laden
	////////////////////////////////////////////////////////////////
	ifstream infile(Datei_in.c_str(), ios_base::binary);
	if (!infile.is_open()) {
		cout << Datei_in
			 << "kann nicht zur Umwandlung in ASCII geöffnet werden\n";
		return 1;
	}
	// Zuerst alles was nur einmal geschrieben ist lesen
	char c_textheader[100];
	for (int i = 0; i < Anzahl_Textheaderzeilen; i++) {
		int pos = infile.tellg();
		infile.read((char *) c_textheader, sizeof(char) * 100);
		if (c_textheader[0] != '#') {
			infile.seekg(pos);
			break;
		}
		textheader[i] = c_textheader;
	}
	infile.read((char *) &No_of_Messungen, sizeof(int));
	infile.read((char *) &No_of_Pix, sizeof(int));
	// Die wichtigsten Dimensionen der Felder sind bekannt....
	// speicher allokieren
	Kanal_Nr.resize(No_of_Pix);
	Wellenlaenge.resize(No_of_Pix);
	////////////////////////////////////////////////////////////////////////////
	// Datensätze lesen
	infile.read((char *) Kanal_Nr.data(), sizeof(int)*No_of_Pix);
	infile.read((char *) Wellenlaenge.data(), sizeof(float)*No_of_Pix);
	// Messungsspezifische Daten hinschreiben
	Nadirdaten.reserve(No_of_Messungen);
	for (int i = 0; i < No_of_Messungen; i++) {
		Nadir_Datensatz nds;
		nds.read_from_mpl_binary(&infile, No_of_Pix);
		Nadirdaten.push_back(nds);
	}
	//cout<<"Laden fertig\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Load_Nadir_l_mpl_binary();
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                          SPEICHERN                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//int Save_Limb_Ascii();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Save_Limb_Ascii(string Datei_out,
					string textheader[31], int &no_of_alt, int &no_of_pix,
					int Orbitstate[5], int Datum[6],
					float Center_Lat_Lon[10], float &orbit_phase,
					std::vector<float> &Wellenlaengen,
					std::vector<Limb_Datensatz> &Limbdaten)
{
	int lang_textheader = 31;
	////////////////////////////////////////////////////////////////////////////
	//
	//Schritt 2 Speichern der Daten in Ausgabedatei
	//
	////////////////////////////////////////////////////////////////////////////
	//cout<<"speichern anfang\n";
	ofstream outfile(Datei_out.c_str());
	if (!outfile.is_open()) {
		cout << Datei_out << " konnte nicht zum schreiben geöffnet werden\n";
		return 1;
	}
	// Datenstrukturen einlesen
	for (int i = 0; i < lang_textheader; i++) {
		//cout<<textheader[i]<<"\n";
		outfile << textheader[i]
				<< "\n";
	}
	outfile << no_of_alt;
	outfile << "\t";
	outfile << no_of_pix;
	outfile << "\n";
	for (int i = 0; i < 5; i++) {
		outfile << Orbitstate[i] << "\t";
	}
	outfile << "\n";
	for (int i = 0; i < 6; i++) {
		outfile << Datum[i] << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Sub_Sat_Lat << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Sub_Sat_Lon << "\t";
	}
	outfile << "\n";
	outfile << "                               " << orbit_phase << "\n";
	outfile << "                               ";
	for (int i = 0; i < 10; i++) {
		outfile << Center_Lat_Lon[i] << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TP_Lat << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TP_Lon << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Tangentenhoehe << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TP_SZA << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TP_SAA << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TP_LOS_Zenit << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TOA_SZA << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TOA_SAA << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_TOA_LOS_Zenit << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Sat_SZA << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Sat_SAA << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Sat_LOS_Zenit << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Sat_Hoehe << "\t";
	}
	outfile << "\n";
	outfile << "                               ";
	for (int i = 0; i < no_of_alt; i++) {
		outfile << Limbdaten[i].m_Erdradius << "\t";
	}
	outfile << "\n";
	// Großes Feld ausschreiben  ( letzte Zeile immer ohne \n
	int l = no_of_pix - 1;
	for (int i = 0; i < l; i++) {
		outfile << Wellenlaengen[i];
		for (int j = 0; j < no_of_alt; j++) {
			outfile << "\t" << Limbdaten[j].m_radiance[i];
		}
		outfile << "\n";
	}
	outfile << Wellenlaengen[l];
	for (int j = 0; j < no_of_alt; j++) {
		outfile << "\t" << Limbdaten[j].m_radiance[l];
	}
	outfile << "\nERRORS\n";
	for (int i = 0; i < l; i++) {
		outfile << Wellenlaengen[i];
		for (int j = 0; j < no_of_alt; j++) {
			outfile << "\t" << Limbdaten[j].m_error[i];
		}
		outfile << "\n";
	}
	outfile << Wellenlaengen[l];
	for (int j = 0; j < no_of_alt; j++) {
		outfile << "\t" << Limbdaten[j].m_error[l];
	}
	//ausgabe ende
	//cout<<"Speichern fertig\n";
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
// ENDE int Save_Limb_Ascii();
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//int Save_Limb_l_mpl_binary();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Save_Limb_l_mpl_binary(string Datei_out,
						   string textheader[31], int &no_of_alt,
						   int &no_of_pix, int Orbitstate[5], int Datum[6],
						   float Center_Lat_Lon[10], float &orbit_phase,
						   float*& Wellenlaengen, Limb_Datensatz*& Limbdaten)
{
	int lang_textheader = 31;
	////////////////////////////////////////////////////////////////////////////
	//
	//Schritt 3 Speichern der Daten in Ausgabedatei
	//
	////////////////////////////////////////////////////////////////////////////
	ofstream outfile(Datei_out.c_str(), ios_base::binary);
	for (int i = 0; i < lang_textheader; i++) {
		outfile.write((char *) textheader[i].c_str(), sizeof(char) * 100);
		//   cout<<textheader[i].c_str()<<"\n";
		//   cout<<i<<"\n";
	}
	outfile.write((char *) &no_of_alt, sizeof(int));
	outfile.write((char *) &no_of_pix, sizeof(int));
	outfile.write((char *) Orbitstate, sizeof(int) * 5);
	outfile.write((char *) Datum, sizeof(int) * 6);
	outfile.write((char *) Center_Lat_Lon, sizeof(float) * 10);
	outfile.write((char *) &orbit_phase, sizeof(float));
	outfile.write((char *) Wellenlaengen, sizeof(float)*no_of_pix);
	for (int i = 0; i < no_of_alt; i++) {
		outfile.write((char *) &Limbdaten[i].m_Sub_Sat_Lat, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Sub_Sat_Lon, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TP_Lat, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TP_Lon, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Tangentenhoehe, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TP_SZA, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TP_SAA, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TP_LOS_Zenit, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TOA_SZA, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TOA_SAA, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_TOA_LOS_Zenit, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Sat_SZA, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Sat_SAA, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Sat_LOS_Zenit, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Sat_Hoehe, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_Erdradius, sizeof(float));
		outfile.write((char *) &Limbdaten[i].m_radiance[0], sizeof(float)*no_of_pix);
		outfile.write((char *) &Limbdaten[i].m_error[0], sizeof(float)*no_of_pix);
	}
	//cout<<"Speichern fertig\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Save_Limb_l_mpl_binary();
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//int Save_Nadir_Ascii();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Save_Nadir_Ascii(string Datei_out,
					 string textheader[7], int No_of_Messungen, int No_of_Pix,
					 std::vector<int> &Kanal_Nr, std::vector<float> &Wellenlaenge,
					 std::vector<Nadir_Datensatz> &Nadirdaten)
{
	int Anzahl_Textheaderzeilen = 7;
	////////////////////////////////////////////////////////////
	// ASCII Datei schreiben
	////////////////////////////////////////////////////////////
	ofstream outfile(Datei_out.c_str());
	if (!outfile.is_open()) {
		cout << Datei_out
			 << " kann nicht erzeugt werden in Nadir Rückumwandlung in Ascii\n";
		return 1;
	}
	for (int i = 0; i < Anzahl_Textheaderzeilen; i++) {
		outfile << textheader[i].c_str() << "\n";   //header von datei
	}
	outfile << No_of_Messungen << "\n";
	for (int i = 0; i < No_of_Messungen; i++) {
		//Zeile mit Messungsnummer...von 0 an durchnummeriert
		outfile << "# \t" << Nadirdaten[i].m_Messung_ID << "\n";
		outfile << Nadirdaten[i].m_state_ID << "\n";
		outfile << Nadirdaten[i].m_Jahr << "\t" << Nadirdaten[i].m_Monat
				<< "\t" << Nadirdaten[i].m_Tag << "\t" << Nadirdaten[i].m_Stunde
				<< "\t" << Nadirdaten[i].m_Minute << "\t"
				<< Nadirdaten[i].m_Sekunde << "\n";
		outfile << Nadirdaten[i].m_SZA_TOA[0] << "\t" << Nadirdaten[i].m_SAA_TOA[0]
				<< "\t" << Nadirdaten[i].m_SZA_TOA[1] << "\t" << Nadirdaten[i].m_SAA_TOA[1]
				<< "\t" << Nadirdaten[i].m_SZA_TOA[2] << "\t" << Nadirdaten[i].m_SAA_TOA[2] << "\n";
		outfile << Nadirdaten[i].m_LOS_Zenit_Winkel[0] << "\t" << Nadirdaten[i].m_LOS_Azimut_Winkel[0]
				<< "\t" << Nadirdaten[i].m_LOS_Zenit_Winkel[1] << "\t" << Nadirdaten[i].m_LOS_Azimut_Winkel[1]
				<< "\t" << Nadirdaten[i].m_LOS_Zenit_Winkel[2] << "\t" << Nadirdaten[i].m_LOS_Azimut_Winkel[2] << "\n";
		outfile << Nadirdaten[i].m_Hoehe << "\t" << Nadirdaten[i].m_Sat_Lat
				<< "\t" << Nadirdaten[i].m_Sat_Lon << "\t"
				<< Nadirdaten[i].m_Sat_Erdradius << "\t"
				<< Nadirdaten[i].m_orbit_phase << "\n";
		outfile << Nadirdaten[i].m_geo_nadir_corner_lat[0] << "\t" << Nadirdaten[i].m_geo_nadir_corner_lon[0]
				<< "\t" << Nadirdaten[i].m_geo_nadir_corner_lat[1] << "\t" << Nadirdaten[i].m_geo_nadir_corner_lon[1]
				<< "\t" << Nadirdaten[i].m_geo_nadir_corner_lat[2] << "\t" << Nadirdaten[i].m_geo_nadir_corner_lon[2]
				<< "\t" << Nadirdaten[i].m_geo_nadir_corner_lat[3] << "\t" << Nadirdaten[i].m_geo_nadir_corner_lon[3]
				<< "\t" << Nadirdaten[i].m_geo_nadir_center_lat << "\t" << Nadirdaten[i].m_geo_nadir_center_lon << "\n";
		outfile << Nadirdaten[i].m_Integrationszeit << "\n";
		outfile << No_of_Pix << "\n";
		for (int j = 0; j < No_of_Pix; j++) {
			outfile << Kanal_Nr[j] << "\t" << Wellenlaenge[j] << "\t"
					<< Nadirdaten[i].m_radiance[j] << "\n";
		}
		outfile << "ERRORS\n";
		for (int j = 0; j < No_of_Pix; j++) {
			outfile << Kanal_Nr[j] << "\t" << Wellenlaenge[j] << "\t"
					<< Nadirdaten[i].m_error[j];
			if ((j < (No_of_Pix - 1)) || (i < No_of_Messungen - 1)) {
				outfile << "\n";
			}
		}

	}// ende for i
	//cout<<"Speichern fertig\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Save_Nadir_Ascii();
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//int Save_Nadir_l_mpl_binary();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Save_Nadir_n_mpl_binary(string Datei_out,
							string textheader[7], int No_of_Messungen,
							int No_of_Pix, int *Kanal_Nr, float *Wellenlaenge,
							Nadir_Datensatz *Nadirdaten)
{
	int Anzahl_Textheaderzeilen = 7;
	///////////////////////////////////////////////////////////////////////
	// n_mpl_binary Datei speichern
	//////////////////////////////////////////////////////////////////////
	ofstream outfile(Datei_out.c_str(), ios_base::binary);
	// Zuerst alles was nur einmal geschrieben werden muss hinschreiben
	for (int i = 0; i < Anzahl_Textheaderzeilen; i++) {
		outfile.write((char *) textheader[i].c_str(), sizeof(char) * 100);
	}
	outfile.write((char *) &No_of_Messungen, sizeof(int));
	outfile.write((char *) &No_of_Pix, sizeof(int));
	outfile.write((char *) Kanal_Nr, sizeof(int)*No_of_Pix);
	outfile.write((char *) Wellenlaenge, sizeof(float)*No_of_Pix);
	// Messungsspezifische Daten hinschreiben
	for (int i = 0; i < No_of_Messungen; i++) {
		outfile.write((char *) &Nadirdaten[i].m_Messung_ID, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_state_ID, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_Jahr, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_Monat, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_Tag, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_Stunde, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_Minute, sizeof(int));
		outfile.write((char *) &Nadirdaten[i].m_Sekunde, sizeof(float));
		outfile.write((char *) Nadirdaten[i].m_SZA_TOA, sizeof(float) * 3);
		outfile.write((char *) Nadirdaten[i].m_SAA_TOA, sizeof(float) * 3);
		outfile.write((char *) Nadirdaten[i].m_LOS_Zenit_Winkel, sizeof(float) * 3);
		outfile.write((char *) Nadirdaten[i].m_LOS_Azimut_Winkel, sizeof(float) * 3);
		outfile.write((char *) &Nadirdaten[i].m_Hoehe, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_Sat_Lat, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_Sat_Lon, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_Sat_Erdradius, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_orbit_phase, sizeof(float));
		outfile.write((char *) Nadirdaten[i].m_geo_nadir_corner_lat, sizeof(float) * 4);
		outfile.write((char *) Nadirdaten[i].m_geo_nadir_corner_lon, sizeof(float) * 4);
		outfile.write((char *) &Nadirdaten[i].m_geo_nadir_center_lat, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_geo_nadir_center_lon, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_Integrationszeit, sizeof(float));
		outfile.write((char *) &Nadirdaten[i].m_radiance[0], sizeof(float)*No_of_Pix);
		outfile.write((char *) &Nadirdaten[i].m_error[0], sizeof(float)*No_of_Pix);
	}
	//cout<<"BINARY Nadir speichern fertig\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Save_Nadir_l_mpl_binary();
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                          KONVERTIEREN                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////int Limb_Ascii_2_l_mpl_binary();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Limb_Ascii_2_l_mpl_binary(string Datei_in, string Datei_out)
{
	///////////////////////////////////////////////////////
	// Datenstrukturen bereitstellen
	//////////////////////////////////////////////////////
	string textheader[31];
	int no_of_alt;
	int no_of_pix;
	int Orbitstate[5];
	int Datum[6];
	float Center_Lat_Lon[10];
	float orbit_phase;
	float *Wellenlaengen = 0;
	Limb_Datensatz *Limbdaten = 0;
	/////////////////////////////////////////////////////////////
	// LADEN
	//////////////////////////////////////////////////////////////
	//cout<<"Lade Ascii\n";
	Load_Limb_Ascii(Datei_in,
					textheader, no_of_alt, no_of_pix,
					Orbitstate, Datum,
					Center_Lat_Lon, orbit_phase,
					Wellenlaengen,
					Limbdaten);
	//  cout<<no_of_alt<<"\n";
	//  cout<<no_of_pix<<"\n";
	//  cout<<Wellenlaengen[30]<<"\n";
	/////////////////////////////////////////////////////////////
	// SPEICHERN
	/////////////////////////////////////////////////////////////
	//cout<<"Speichere binär\n";
	Save_Limb_l_mpl_binary(Datei_out,
						   textheader, no_of_alt, no_of_pix,
						   Orbitstate, Datum,
						   Center_Lat_Lon, orbit_phase,
						   Wellenlaengen,
						   Limbdaten);
	/////////////////////////////////////////////////////////////
	// AUFRÄUMEN
	//////////////////////////////////////////////////////////////
	delete[] Wellenlaengen;
	delete[] Limbdaten;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Limb_Ascii_2_l_mpl_binary();
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//int Nadir_Ascii_2_n_mpl_binary();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Nadir_Ascii_2_n_mpl_binary(string Datei_in, string Datei_out)
{
	///////////////////////////////////////////////////////
	// Datenstrukturen bereitstellen
	//////////////////////////////////////////////////////
	string textheader[7];
	int No_of_Messungen;  // Das variiert bei nadir doch...standard ist 65
	int No_of_Pix;
	int *Kanal_Nr = 0;
	float *Wellenlaenge = 0;
	Nadir_Datensatz *Nadirdaten = 0;
	/////////////////////////////////////////////////////////////
	// LADEN
	//////////////////////////////////////////////////////////////
	//cout<<"Lade Nadir\n";
	//sleep(1);
	Load_Nadir_Ascii(Datei_in,
					 textheader, No_of_Messungen, No_of_Pix,
					 Kanal_Nr, Wellenlaenge, Nadirdaten);
	/////////////////////////////////////////////////////////////
	// SPEICHERN
	/////////////////////////////////////////////////////////////
	//cout<<"SAVE Nadir\n";
	//sleep(1);
	Save_Nadir_n_mpl_binary(Datei_out,
							textheader, No_of_Messungen, No_of_Pix,
							Kanal_Nr, Wellenlaenge, Nadirdaten);
	/////////////////////////////////////////////////////////////
	// AUFRÄUMEN
	//////////////////////////////////////////////////////////////
	//cout<<"Räume auf Nadir\n";
	//sleep(1);
	delete[] Kanal_Nr;
	delete[] Wellenlaenge;
	delete[] Nadirdaten;
	//cout<<"Umwandlung fertig Nadir\n";
	//sleep(1);
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Nadir_Ascii_2_n_mpl_binary();
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//int Limb_l_mpl_binary_2_Ascii();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Limb_l_mpl_binary_2_Ascii(string Datei_in, string Datei_out)
{
	///////////////////////////////////////////////////////
	// Datenstrukturen bereitstellen
	//////////////////////////////////////////////////////
	string textheader[31];
	int no_of_alt;
	int no_of_pix;
	int Orbitstate[5];
	int Datum[6];
	float Center_Lat_Lon[10];
	float orbit_phase;
	std::vector<float> Wellenlaengen;
	std::vector<Limb_Datensatz> Limbdaten;
	/////////////////////////////////////////////////////////////
	// LADEN
	//////////////////////////////////////////////////////////////
	//cout<<"Lade binär \n";
	Load_Limb_l_mpl_binary(Datei_in,
						   textheader, no_of_alt, no_of_pix,
						   Orbitstate, Datum,
						   Center_Lat_Lon,
						   orbit_phase,
						   Wellenlaengen,
						   Limbdaten);
	/////////////////////////////////////////////////////////////
	// SPEICHERN
	/////////////////////////////////////////////////////////////
	//cout<<"Speichere Ascii\n";
	Save_Limb_Ascii(Datei_out,
					textheader, no_of_alt, no_of_pix,
					Orbitstate, Datum,
					Center_Lat_Lon,
					orbit_phase,
					Wellenlaengen,
					Limbdaten);

	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Limb_l_mpl_binary_2_Ascii();
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//int Nadir_n_mpl_binary_2_Ascii();Funktionsstart
////////////////////////////////////////////////////////////////////////////////
int Nadir_n_mpl_binary_2_Ascii(string Datei_in, string Datei_out)
{
	///////////////////////////////////////////////////////
	// Datenstrukturen bereitstellen
	//////////////////////////////////////////////////////
	string textheader[7];
	int No_of_Messungen;  // Das variiert bei nadir doch...standard ist 65
	int No_of_Pix;
	std::vector<int> Kanal_Nr;
	std::vector<float> Wellenlaenge;
	std::vector<Nadir_Datensatz> Nadirdaten;
	/////////////////////////////////////////////////////////////
	// LADEN
	//////////////////////////////////////////////////////////////
	//cout<<"Lade Nadir\n";
	Load_Nadir_n_mpl_binary(Datei_in,
							textheader, No_of_Messungen, No_of_Pix,
							Kanal_Nr, Wellenlaenge, Nadirdaten);
	/////////////////////////////////////////////////////////////
	// SPEICHERN
	/////////////////////////////////////////////////////////////
	//cout<<"Speichere Nadir\n";
	Save_Nadir_Ascii(Datei_out,
					 textheader, No_of_Messungen, No_of_Pix,
					 Kanal_Nr, Wellenlaenge, Nadirdaten);
	//cout<<"Umwandlung abgeschlossen\n";
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE int Nadir_n_mpl_binary_2_Ascii();
////////////////////////////////////////////////////////////////////////////////
