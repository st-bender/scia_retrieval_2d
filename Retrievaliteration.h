/*
 * Retrievaliteration.h
 *
 *  Created on: 14.09.2010
 *      Author: martin
 */

#ifndef RETRIEVALITERATION_HH_
#define RETRIEVALITERATION_HH_

#include "MPL_Matrix.h"

int Retrievaliteration(MPL_Matrix &Dichten,
					   MPL_Matrix &Dichten_apriori,
					   MPL_Matrix &Saeulendichten,
					   MPL_Matrix &S_apriori,
					   MPL_Matrix &S_y,
					   MPL_Matrix &S_Breite,
					   MPL_Matrix &S_Hoehe,
					   const double &Lambda_Breite,
					   const double &Lambda_Hoehe,
					   MPL_Matrix &AMF,
					   class Konfiguration &Konf);
int Retrievaliteration_old(MPL_Matrix &Dichten,  // Zur Sicherheit behalten
						   // Iteration über anpassen des apriori
						   MPL_Matrix &Dichten_apriori,
						   MPL_Matrix &Saeulendichten,
						   MPL_Matrix &S_apriori,
						   MPL_Matrix &S_y,
						   MPL_Matrix &S_Breite,
						   MPL_Matrix &S_Hoehe,
						   MPL_Matrix &S_letzte_Hoehe,
						   const double &Lambda_Breite,
						   const double &Lambda_Hoehe,
						   MPL_Matrix &AMF,
						   class Konfiguration &Konf);
#ifdef HAVE_EIGEN3
int Retrievaliteration_Eigen(MPL_Matrix &Dichten,
					   MPL_Matrix &Dichten_apriori,
					   MPL_Matrix &Saeulendichten,
					   MPL_Matrix &S_apriori,
					   MPL_Matrix &S_y,
					   MPL_Matrix &S_Breite,
					   MPL_Matrix &S_Hoehe,
					   const double &Lambda_Breite,
					   const double &Lambda_Hoehe,
					   MPL_Matrix &AMF,
					   class Konfiguration &Konf);
#endif /* HAVE_EIGEN3 */

#endif /* RETRIEVALITERATION_HH_ */
