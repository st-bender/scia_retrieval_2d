/*
 * Data_IO.h
 *
 *  Created on: 09.08.2010
 *      Author: martin
 */

////////////////////////////////////////////////////////////////////////////////
// UMWANDLUNGSROUTINEN VON SCIA_ASCII_DATEN IN BINÄRDATEN UND UMGEKEHRT
//
// Dieses Funktionspaket übernimmt die Implementierungen der Funktionen aus
// LN_Umwandlung.h und MPLBIN2ASCII.h und packt diese in ein einziges Paket
// zusammen...dabei werden die Teilschritte des Ladens und des Speicherns
// als selbstständige Funktionen Implementiert
////////////////////////////////////////////////////////////////////////////////

#ifndef LIMBNADIR_IO_HH_
#define LIMBNADIR_IO_HH_
//Zunächst die Struktur, in der eine Nadirmessung gespeichert werden kann
// FUNKTIONDEKLARATIONEN WEITER UNTEN

#include <string>
#include <vector>
#include <fstream>

#ifndef NADIR_DATENSATZ_HH_
#define NADIR_DATENSATZ_HH_
//struktur zur Speicherung von Daten aus Nadirdatei
// Der Header, in dem das Ausgangsprodukt steht, wird nur einmal in die Datei
// gespeichert, muss also nicht redundant in jeden Datensatz mit aufgenommen
// werden
class Nadir_Datensatz
{
public:
	//Zuweisungsoperator
	Nadir_Datensatz &operator=(const Nadir_Datensatz &RHS);
	void read_from_mpl_binary(std::ifstream *stream, int no_of_pix);
	int m_Messung_ID; //das was hinter # steht
	int m_state_ID;       //Art des Messungszustands, steht. z.b. im Sciabuch
	int m_Jahr;
	int m_Monat;
	int m_Tag;
	int m_Stunde;
	int m_Minute;
	float m_Sekunde;
	float m_SZA_TOA[3];
	float m_SAA_TOA[3];
	float m_LOS_Zenit_Winkel[3];
	float m_LOS_Azimut_Winkel[3];
	float m_Hoehe;
	float m_Sat_Lat;
	float m_Sat_Lon;
	float m_Sat_Erdradius;
	float m_orbit_phase;
	float m_geo_nadir_corner_lat[4];
	float m_geo_nadir_corner_lon[4];
	float m_geo_nadir_center_lat;
	float m_geo_nadir_center_lon;
	float m_Integrationszeit;
	int m_N_radiances;
	// float m_no_of_pix;                   wird einmal Zentral gespeichert
	// float* m_Pixelnummer             wird einmal Zentral gespeichert
	// float* m_Wellenlaengen           wird einmal Zentral gespeichert
	std::vector<float> m_radiance; // eher Teilchen/(cm^2nm)
					   // der für dieses Feld allokierte Speicher muss wieder
					   // gelöscht werden
	std::vector<float> m_error;

};
//implementation der Übergabe
inline Nadir_Datensatz &Nadir_Datensatz::operator=(const Nadir_Datensatz &RHS)
{
	if (this == &RHS) {
		return *this;
	}
	m_Messung_ID = RHS.m_Messung_ID;
	m_state_ID = RHS.m_state_ID;
	m_Jahr = RHS.m_Jahr;
	m_Monat = RHS.m_Monat;
	m_Tag = RHS.m_Tag;
	m_Stunde = RHS.m_Stunde;
	m_Minute = RHS.m_Minute;
	m_Sekunde = RHS.m_Sekunde;
	for (int i = 0; i < 3; i++) {
		m_SZA_TOA[i] = RHS.m_SZA_TOA[i];
		m_SAA_TOA[i] = RHS.m_SAA_TOA[i];
		m_LOS_Zenit_Winkel[i] = RHS.m_LOS_Zenit_Winkel[i];
		m_LOS_Azimut_Winkel[i] = RHS.m_LOS_Azimut_Winkel[i];
	}
	m_Hoehe = RHS.m_Hoehe;
	m_Sat_Lat = RHS.m_Sat_Lat;
	m_Sat_Lon = RHS.m_Sat_Lon;
	m_Sat_Erdradius = RHS.m_Sat_Erdradius;
	m_orbit_phase = RHS.m_orbit_phase;
	for (int i = 0; i < 4; i++) {
		m_geo_nadir_corner_lat[i] = RHS.m_geo_nadir_corner_lat[i];
		m_geo_nadir_corner_lon[i] = RHS.m_geo_nadir_corner_lon[i];
	}
	m_geo_nadir_center_lat = RHS.m_geo_nadir_center_lat;
	m_geo_nadir_center_lon = RHS.m_geo_nadir_center_lon;
	m_Integrationszeit = RHS.m_Integrationszeit;
	m_N_radiances = RHS.m_N_radiances;

	std::copy(RHS.m_radiance.begin(), RHS.m_radiance.end(),
			m_radiance.begin());
	std::copy(RHS.m_error.begin(), RHS.m_error.end(),
			m_error.begin());

	return *this;
}


//Struktur zur Speicherung von Limb Daten
class Limb_Datensatz
{
public:
	//Zuweisungsoperator
	Limb_Datensatz &operator=(const Limb_Datensatz &RHS);
	void read_from_mpl_binary(std::ifstream *stream, int no_of_pix);
	float m_Sub_Sat_Lat;
	float m_Sub_Sat_Lon;
	float m_TP_Lat;
	float m_TP_Lon;
	float m_Tangentenhoehe;
	float m_TP_SZA;
	float m_TP_SAA;
	float m_TP_LOS_Zenit;
	float m_TOA_SZA;
	float m_TOA_SAA;
	float m_TOA_LOS_Zenit;
	float m_Sat_SZA;
	float m_Sat_SAA;
	float m_Sat_LOS_Zenit;
	float m_Sat_Hoehe;
	float m_Erdradius;
	int m_N_radiances;
	std::vector<float> m_radiance;
	std::vector<float> m_error;
};

inline Limb_Datensatz &Limb_Datensatz::operator=(const Limb_Datensatz &RHS)
{
	if (this == &RHS) {
		return *this;
	}
	m_Sub_Sat_Lat = RHS.m_Sub_Sat_Lat;
	m_Sub_Sat_Lon = RHS.m_Sub_Sat_Lon;
	m_TP_Lat = RHS.m_TP_Lat;
	m_TP_Lon = RHS.m_TP_Lon;
	m_Tangentenhoehe = RHS.m_Tangentenhoehe;
	m_TP_SZA = RHS.m_TP_SZA;
	m_TP_SAA = RHS.m_TP_SAA;
	m_TP_LOS_Zenit = RHS.m_TP_LOS_Zenit;
	m_TOA_SZA = RHS.m_TOA_SZA;
	m_TOA_SAA = RHS.m_TOA_SAA;
	m_TOA_LOS_Zenit = RHS.m_TOA_LOS_Zenit;
	m_Sat_SZA = RHS.m_Sat_SZA;
	m_Sat_SAA = RHS.m_Sat_SAA;
	m_Sat_LOS_Zenit = RHS.m_Sat_LOS_Zenit;
	m_Sat_Hoehe = RHS.m_Sat_Hoehe;
	m_Erdradius = RHS.m_Erdradius;
	m_N_radiances = RHS.m_N_radiances;

	std::copy(RHS.m_radiance.begin(), RHS.m_radiance.end(),
			m_radiance.begin());
	std::copy(RHS.m_error.begin(), RHS.m_error.end(),
			m_error.begin());

	return *this;
}
#endif /* NADIR_DATENSATZ_HH_ */

template<typename T>
std::istream& binary_read(std::istream* stream, T& value, size_t N = 1) {
	return stream->read(reinterpret_cast<char *>(&value), N * sizeof(T));
}

// Laden
int Load_Limb_Ascii(std::string Datei_in,
					std::string textheader[31], int &no_of_alt, int &no_of_pix,
					int Orbitstate[5], int Datum[6], float Center_Lat_Lon[10],
					float &orbit_phase, float*& Wellenlaengen,
					Limb_Datensatz*& Limbdaten);

int Load_Limb_l_mpl_binary(std::string Datei_in,
						   std::string textheader[31], int &no_of_alt,
						   int &no_of_pix, int Orbitstate[5], int Datum[6],
						   float Center_Lat_Lon[10], float &orbit_phase,
						   std::vector<float> &Wellenlaengen,
						   std::vector<Limb_Datensatz> &Limbdaten);

int Load_Nadir_Ascii(std::string Datei_in,
					 std::string textheader[7], int &No_of_Messungen, int &No_of_Pix,
					 int*& Kanal_Nr, float &orbit_phase,
					 float*& Wellenlaenge, Nadir_Datensatz*& Nadirdaten);
int Load_Nadir_n_mpl_binary(std::string Datei_in,
							std::string textheader[7], int &No_of_Messungen,
							int &No_of_Pix, int*& Kanal_Nr,
							float*& Wellenlaenge, Nadir_Datensatz*& Nadirdaten);
//Speichern
int Save_Limb_Ascii(std::string Datei_out,
					std::string textheader[31], int &no_of_alt, int &no_of_pix,
					int Orbitstate[5], int Datum[6], float Center_Lat_Lon[10],
					float &orbit_phase, std::vector<float> Wellenlaengen,
					std::vector<Limb_Datensatz> &Limbdaten);

int Save_Limb_l_mpl_binary(std::string Datei_out,
						   std::string textheader[31], int &no_of_alt,
						   int &no_of_pix, int Orbitstate[5], int Datum[6],
						   float Center_Lat_Lon[10], float &orbit_phase,
						   float*& Wellenlaengen, Limb_Datensatz*& Limbdaten);

int Save_Nadir_Ascii(std::string Datei_out,
					 std::string textheader[7], int No_of_Messungen, int No_of_Pix,
					 int *Kanal_Nr, float &orbit_phase, float *Wellenlaenge,
					 Nadir_Datensatz *Nadirdaten);

int Save_Nadir_n_mpl_binary(std::string Datei_out,
							std::string textheader[7], int No_of_Messungen,
							int No_of_Pix, int *Kanal_Nr, float *Wellenlaenge,
							Nadir_Datensatz *Nadirdaten);
// Konvertieren
int Limb_Ascii_2_l_mpl_binary(std::string Datei_in, std::string Datei_out);
int Nadir_Ascii_2_n_mpl_binary(std::string Datei_in, std::string Datei_out);
int Limb_l_mpl_binary_2_Ascii(std::string Datei_in, std::string Datei_out);
int Nadir_n_mpl_binary_2_Ascii(std::string Datei_in, std::string Datei_out);

#endif /* LIMBNADIR_IO_HH_ */
