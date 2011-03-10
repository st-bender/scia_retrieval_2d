/*
 * Messungs_ausschliessende_Tests.h
 *
 *  Created on: 13.10.2010
 *      Author: martin
 */

#ifndef MESSUNGS_AUSSCHLIESSENDE_TESTS_HH_
#define MESSUNGS_AUSSCHLIESSENDE_TESTS_HH_

#include"Messung_Limb.h"
#include"Messung_Nadir.h"
#include<vector>

using namespace std;
//Falls Test_auf_Nachtmessung oder Test_auf_NLC positiv sind, so werden ganze Dateien aussortiert
//Bei den selten vorkommenden Geolocations ist zumeist nur eine Messung falsch (die einzige, die je aufgefallen ist)
int Test_auf_Nachtmessung_Limb(Messung_Limb& Tropo, bool &ist_Nachtmessung);
int Test_auf_Nachtmessung_Limb_meso_thermo(Messung_Limb& niedrigste_hoehe, bool &ist_Nachtmessung);
int Test_auf_NLC_Limb(vector<Messung_Limb> &Rohdaten, bool &ist_NLC_Messung);
int Test_auf_korrekte_geolocations_Limb(vector<Messung_Limb> &Rohdaten,int* counter_Winkel_nicht_ok);
int Test_auf_Nachtmessung_Nadir(Messung_Nadir* Rohdaten,int Anzahl_Datensaetze,bool &ist_Nachtmessung);


#endif /* MESSUNGS_AUSSCHLIESSENDE_TESTS_HH_ */
