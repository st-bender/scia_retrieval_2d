/*
 * Messungs_ausschliessende_Tests.cpp
 *
 *  Created on: 13.10.2010
 *      Author: martin
 */

#include"Messungs_ausschliessende_Tests.h"
#include"Messung_Limb.h"
#include"Messung_Nadir.h"
#include<vector>
#include"MPL_Vektor.h"
#include"Koordinatentransformation.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
// Funktionsstart Test_auf_Nachtmessung_Limb
/////////////////////////////////////////////////////////////////////////
int Test_auf_Nachtmessung_Limb(Messung_Limb &Tropo, bool &ist_Nachtmessung)
{
	// Der Test überprüft im Intervall 290nm bis 295 nm, ob die mittleren Signale bei der Tangentenhöhe -1km (wie bei Marco)
	// einen kritischen Mittelwert unterschreiten (1E10)  bei Tagmessungen wird dieser Wert zumeist
	// mindestenz um einen Faktor 10 überschritten, während bei Nachtmessungen die Unterschreitung noch deutlicher ist
	// 1E7 bis 1E8
	ist_Nachtmessung = false;
	int Index1 = 608;
	int Index2 = 653; //Die Indizes könnte man auch ermitteln, aber die sind ja immer gleich
	double Troposignal = 0;
	for (int i = Index1; i <= Index2; i++) {
		Troposignal += Tropo.m_Intensitaeten[i];
	}
	Troposignal /= (Index2 - Index1 + 1);
	if (Troposignal < 1E10) {
		//cout<<"Nachtmessung detektiert\n";
		//cout<<"Troposignal:"<<Troposignal<<"\n";
		ist_Nachtmessung = true;
	}
	return 0;
}
//////////////////////////////////////////////////////////////////////////
// ENDE Test_auf_Nachtmessung_Limb
/////////////////////////////////////////////////////////////////////////

int Test_auf_Nachtmessung_Limb_meso_thermo(Messung_Limb &niedrigste_hoehe, bool &ist_Nachtmessung)
{
	//Wie im Limbfall nur leider nicht bei -1km TH, weil Messung nicht vorhanden
	// TODO Schwellenenergie heraufinden
	ist_Nachtmessung = false;
	int Index1 = 608;
	int Index2 = 653; //Die Indizes könnte man auch ermitteln, aber die sind ja immer gleich
	double Signal = 0;
	for (int i = Index1; i <= Index2; i++) {
		Signal += niedrigste_hoehe.m_Intensitaeten[i];
	}
	Signal /= (Index2 - Index1 + 1);
	//cerr<<"Signal 290nm-295nm 53km:"<<Signal<<"\n";
	double threshold = 5E10; // ist doch fast gleich
	if (Signal < threshold) {
		//cerr<<"Nachtmessung detektiert\n";
		//cerr<<"Signal 290nm-295nm 53km:"<<Signal<<"\n";
		ist_Nachtmessung = true;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////
// Funktionsstart Test_auf_NLC_Limb
/////////////////////////////////////////////////////////////////////////
int Test_auf_NLC_Limb(vector<Messung_Limb> &Rohdaten, bool &ist_NLC_Messung)
{
	//Die Intensitäten in einem Breiten Bereich zwischen etwa 1/4 des Spektrums und 3/4 des Spektrums werden
	//für die 7 Tangentenhöhen gemittelt. Das Signal fällt, falls keine Wolken in der Mesosphäre sind,
	//mit zunehmender Höhe ziemlich Stark ab.
	//Wolken Reflektieren das Sonnenlicht. Daher gibt es bei Wolken in der Mesosphäre (Polar Mesospheric Clouds
	//PMCs, oder auch Nachtleuchtende Wolken NLCs genannt) einen Peak, der ungefähr auf der Höhe der Mesopause
	//(im Sommer bei 84 km) liegt. Fällt das Signal also nicht monoton mit der Hoehe, so wurde eine PMC detektiert.
	// Das erste Viertel abzuschneiden ist eine gute Idee, da dort das Signal extrem verrauscht und stark ist.
	// 71 ist in 0;  91 ist in 6

	int Index1 = 242;
	int Index2 = 727; //Die Indizes könnte man auch ermitteln, aber die sind ja immer gleich
	ist_NLC_Messung = false;
	double mittleres_Signal[7];
	for (int Hoehenlevel = 0; Hoehenlevel < 7; Hoehenlevel++) {
		mittleres_Signal[Hoehenlevel] = 0;
		for (int i = Index1; i <= Index2; i++) {
			mittleres_Signal[Hoehenlevel] += Rohdaten[Hoehenlevel].m_Intensitaeten[i];
		}
		mittleres_Signal[Hoehenlevel] /= (Index2 - Index1 + 1);
	}
	for (int i = 0; i < 6; i++) {
		if (mittleres_Signal[i + 1] > mittleres_Signal[i]) {
			ist_NLC_Messung = true;
		}
	}
	return 0;
}
//////////////////////////////////////////////////////////////////////////
// ENDE Test_auf_NLC_Limb
/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Funktionsstart Test_auf_korrekte_geolocations_Limb
/////////////////////////////////////////////////////////////////////////
int Test_auf_korrekte_geolocations_Limb(vector<Messung_Limb> &Rohdaten, int *counter_Winkel_nicht_ok)
{
	//Die Positionen des Tangentenpunktes und des Satelliten werden in Karthesischen Koordinaten bestimmt.
	//Der Verbindungsvektor beider Punkte und der Ortsvektor des Tangentenpunkts bilden, wenn man sie vorher normiert
	//im Skalarprodukt den cosinus des Winkels a; dieser entspricht dem sinus von b=(a-90°).
	//Für kleine Winkel gilt sin x= tan x, sodass man aus der Abweichung von 90° eine Steigung erhält und mit der
	//Sehnenlänge in der Atmosphäre eine Höhenverschiebung des Tangentenpunkts
	//Derartige Abschätzungen ergeben eine Verschiebung der Tangentenhoehe um 1,13 km bei einer Abweichung von 0.02°
	//                                                                                                      und um 2,84 km bei einer Abweichung von 0.05°
	//Besser gesagt, liegt dann der Tangentenpunkt nicht am angegebenen Tangentenpunkt...das ist nicht linear und kompliziert
	//aber auch kleiner als der oben abgeschätzte Worst-Case (erst bei Messungen mit 5° Unterschied, kam eine verdächtige
	//Fehlermeldung)
	//Die meisten Messwerte liegen zwischen 0.01 und 0.02
	const double pi = 3.1415926535897;
	unsigned int i = 0; //i wie INDEX
	while (i < Rohdaten.size()) {
		//Ortsvektoren bestimmen
		MPL_Vektor Ort_Sat(3), Ort_TP(3);
		Umwandlung_Kugel_in_Karthesisch(Rohdaten[i].m_Erdradius + Rohdaten[i].m_Hoehe_Sat,
										Rohdaten[i].m_Longitude_Sat,
										Rohdaten[i].m_Lattidude_Sat,
										Ort_Sat(0), Ort_Sat(1), Ort_Sat(2));
		Umwandlung_Kugel_in_Karthesisch(Rohdaten[i].m_Erdradius + Rohdaten[i].m_Hoehe_TP,
										Rohdaten[i].m_Longitude_TP,
										Rohdaten[i].m_Lattidude_TP,
										Ort_TP(0), Ort_TP(1), Ort_TP(2));

		//Verbindungsvektor
		MPL_Vektor Verbindung(3);
		Verbindung = Ort_TP - Ort_Sat;
		//normieren
		MPL_Vektor TP_normiert(3);
		Verbindung = Verbindung / Verbindung.Betrag_ausgeben();
		TP_normiert = Ort_TP / Ort_TP.Betrag_ausgeben();
		//Skalarprodukt bilden
		double SIN_Winkel = Verbindung * TP_normiert;
		double Winkel = 180.0 / pi * asin(SIN_Winkel);
		//aussortieren und counter setzen
		if (Winkel > 0.02) {
			(*counter_Winkel_nicht_ok)++;
			Rohdaten.erase(Rohdaten.begin() + i);
		} else {
			i++;
		}
	} //ende while
	return 0;
}
//////////////////////////////////////////////////////////////////////////
// ENDE Test_auf_korrekte_geolocations_Limb
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// Funktionsstart Test_auf_Nachtmessung_Nadir
/////////////////////////////////////////////////////////////////////////
int Test_auf_Nachtmessung_Nadir(Messung_Nadir *Rohdaten, int Anzahl_Datensaetze, bool &ist_Nachtmessung)
{
	// Der Test überprüft im Intervall 290nm bis 295 nm, ob die mittleren Signale einen kritischen Mittelwert unterschreiten (5E8)
	// bei Tagmessungen wird dieser Wert zumeist mindestenz um einen Faktor 10 überschritten
	// Bei Nadir ist der Grenzwert kleiner, der Rest ist fast gleich....die erste Datei wird zuerst gemessen, ist also am Nordpol
	// die Dunkelste
	// Die letzte Datei ist am Südpol die dunkelste
	// Es wird also erst auf Nord oder Südhalbkugel geprüft
	int n;// n wie Messungnummer
	if (Rohdaten[0].m_Lattitude_Sat > 0) {
		n = 0;
	} else {
		n = Anzahl_Datensaetze - 1;
	}
	///////////////////////////////////////////////// Der Rest ist ziemlich analog zu Limb
	ist_Nachtmessung = false;
	int Index1 = 608;
	int Index2 = 653; //Die Indizes könnte man auch ermitteln, aber die sind ja immer gleich
	double Signal = 0;
	for (int i = Index1; i <= Index2; i++) {
		Signal += Rohdaten[n].m_Intensitaeten[i];
	}
	Signal /= (Index2 - Index1 + 1);
	if (Signal < 5E8) {
		//cout<<"Nachtmessung detektiert\n";
		//cout<<"Signal_Nadir:"<<Signal<<"\n";
		ist_Nachtmessung = true;
	}
	return 0;
}
//////////////////////////////////////////////////////////////////////////
// ENDE Test_auf_Nachtmessung_Nadir
/////////////////////////////////////////////////////////////////////////
