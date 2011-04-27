/*
 * Glaetten_2D.h
 *
 *  Created on: 24.01.2011
 *      Author: martin
 */

#ifndef GLAETTEN_2D_HH_
#define GLAETTEN_2D_HH_

#include <string>
#include <vector>
#include "Ausgewertete_Messung_Limb.h"

int SCD_Glaettung(std::vector<Ausgewertete_Messung_Limb> &AML,
		int Anzahl_Linien,
		std::string limbmesothermo);


#endif /* GLAETTEN_2D_HH_ */
