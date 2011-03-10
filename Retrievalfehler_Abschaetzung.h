/*
 * Retrievalfehler_Abschaetzung.h
 *
 *  Created on: 16.09.2010
 *      Author: martin
 */

#ifndef RETRIEVALFEHLER_ABSCHAETZUNG_HH_
#define RETRIEVALFEHLER_ABSCHAETZUNG_HH_


#include "MPL_Matrix.h"
#include "Konfiguration.h"

using namespace std;


int Retrievalfehler_Abschaetzung(MPL_Matrix&             S_x,
                                                     MPL_Matrix&            Averaging_Kernel_Matrix,
                                                     const MPL_Matrix&   S_apriori,
                                                     const MPL_Matrix&   S_y,
                                                     MPL_Matrix               S_Breite,
                                                     MPL_Matrix               S_Hoehe,
                                                     MPL_Matrix               S_letzte_Hoehe,
                                                     const double&          Lambda_Breite,
                                                     const double&          Lambda_Hoehe,
                                                     MPL_Matrix               AMF,
                                                     const Konfiguration& Konf);
void Matrix_Invertieren(MPL_Matrix& M);


#endif /* RETRIEVALFEHLER_ABSCHAETZUNG_HH_ */
