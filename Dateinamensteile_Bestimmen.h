/*
 * Dateinamensteile_Bestimmen.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 23.09.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef DATEINAMENSTEILE_BESTIMMEN_HH_
#define DATEINAMENSTEILE_BESTIMMEN_HH_


// Sinn dieser Funktion ist es im Batchprozess die richtigen Dateinamen
// zu bestimmen die erzeugt werden sollen
// xxxxx Orbit_ID
// yyyymmdd_hhmm Das Datum der ersten Limbdatei im Orbit

#include<string>

std::string sb_basename(std::string filename);
std::string xxxxx_Bestimmen(std::string Orbitlistenpfad);
std::string yyyymmdd_hhmm_Bestimmen(std::string Name_erste_Limbdatei);


#endif /* DATEINAMENSTEILE_BESTIMMEN_H_ */
