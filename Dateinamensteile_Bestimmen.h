/*
 * Dateinamensteile_Bestimmen.h
 *
 *  Created on: 23.09.2010
 *      Author: martin
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
