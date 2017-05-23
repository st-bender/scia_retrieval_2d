/*
 * Retrievalfehler_Abschaetzung.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 16.09.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef RETRIEVALFEHLER_ABSCHAETZUNG_HH_
#define RETRIEVALFEHLER_ABSCHAETZUNG_HH_


#include "MPL_Matrix.h"

void Retrievalfehler_Abschaetzung(MPL_Matrix &S_x,
								 MPL_Matrix &S_x_meas,
								 MPL_Matrix &Averaging_Kernel_Matrix,
								 const MPL_Matrix &S_apriori,
								 const MPL_Matrix &S_y,
								 MPL_Matrix &S_Breite,
								 MPL_Matrix &S_Hoehe,
								 const double &Lambda_Breite,
								 const double &Lambda_Hoehe,
								 MPL_Matrix &AMF,
								 const class Konfiguration &Konf);
void Matrix_Invertieren(MPL_Matrix &M);


#endif /* RETRIEVALFEHLER_ABSCHAETZUNG_HH_ */
