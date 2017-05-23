/*
 * Glaetten_2D.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2011 Martin Langowski
 *
 * Initial version created on: 24.01.2011
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef GLAETTEN_2D_HH_
#define GLAETTEN_2D_HH_

#include <string>
#include <vector>

int SCD_Glaettung(std::vector<class Ausgewertete_Messung_Limb> &AML,
		int Anzahl_Linien,
		std::string limbmesothermo);


#endif /* GLAETTEN_2D_HH_ */
