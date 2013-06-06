/*
 * Messungs_ausschliessende_Tests.h
 *
 *  Created on: 13.10.2010
 *      Author: martin
 */

#ifndef MESSUNGS_AUSSCHLIESSENDE_TESTS_HH_
#define MESSUNGS_AUSSCHLIESSENDE_TESTS_HH_

#include <vector>

//Falls Test_auf_Nachtmessung oder Test_auf_NLC positiv sind,
//so werden ganze Dateien aussortiert
//Bei den selten vorkommenden Geolocations ist zumeist nur eine Messung falsch
//(die einzige, die je aufgefallen ist)
bool Test_auf_Nachtmessung_Limb(class Messung_Limb &Tropo, class Konfiguration &Konf);
bool Test_auf_Nachtmessung_Limb_meso_thermo(class Messung_Limb &niedrigste_hoehe,
		class Konfiguration &Konf);
bool test_auf_SAA_limb(class Messung_Limb &space, Konfiguration &konf);
bool Test_auf_NLC_Limb(std::vector<class Messung_Limb> &Rohdaten,
		class Konfiguration &Konf);
int Test_auf_korrekte_geolocations_Limb(std::vector<class Messung_Limb> &Rohdaten,
		int &counter_Winkel_nicht_ok);
bool Test_auf_Nachtmessung_Nadir(std::vector<class Messung_Nadir> &Rohdaten,
		int Anzahl_Datensaetze);


#endif /* MESSUNGS_AUSSCHLIESSENDE_TESTS_HH_ */
