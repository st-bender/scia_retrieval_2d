/*
 * Koordinatentransformation.h
 *
 *  Created on: 08.09.2010
 *      Author: martin
 */

#ifndef KOORDINATENTRANSFORMATION_HH_
#define KOORDINATENTRANSFORMATION_HH_

// Funktionen zur Umwandlung von karthesichen Koordinaten in andere Systeme und umgekehrt

// Kugelkoordinaten
void  Umwandlung_Kugel_in_Karthesisch(double r, double phi, double theta, double &x, double &y, double &z);
void  Umwandlung_Karthesisch_in_Kugel(double x, double y, double z, double &r, double &phi, double &theta);
//

#endif /* KOORDINATENTRANSFORMATION_HH_ */
