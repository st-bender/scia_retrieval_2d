/*
 * Retrievaliteration.h
 *
 *  Created on: 14.09.2010
 *      Author: martin
 */

#ifndef RETRIEVALITERATION_HH_
#define RETRIEVALITERATION_HH_

#include "MPL_Matrix.h"
#include "Konfiguration.h"

using namespace std;

int Retrievaliteration(MPL_Matrix &Dichten,
					   MPL_Matrix &Dichten_apriori,
					   MPL_Matrix &Saeulendichten,
					   MPL_Matrix &S_apriori,
					   MPL_Matrix &S_y,
					   MPL_Matrix S_Breite,
					   MPL_Matrix S_Hoehe,
					   const double &Lambda_Breite,
					   const double &Lambda_Hoehe,
					   MPL_Matrix &AMF,
					   Konfiguration &Konf);
int Retrievaliteration_old(MPL_Matrix &Dichten,             // Zur Sicherheit behalten
						   MPL_Matrix &Dichten_apriori,         // Iteration über anpassen des apriori
						   MPL_Matrix &Saeulendichten,
						   MPL_Matrix &S_apriori,
						   MPL_Matrix &S_y,
						   MPL_Matrix S_Breite,
						   MPL_Matrix S_Hoehe,
						   MPL_Matrix S_letzte_Hoehe,
						   const double &Lambda_Breite,
						   const double &Lambda_Hoehe,
						   MPL_Matrix &AMF,
						   Konfiguration &Konf);

#endif /* RETRIEVALITERATION_HH_ */
