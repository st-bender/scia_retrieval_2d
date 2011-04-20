/*
 * constants.h
 *
 *  Created on: 20-Apr-2011
 *      Author: bender-s
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#define PHYS_h 6.62606896e-34 // Planck constant h [Js]
#define PHYS_h_eV 4.13566733e-15 // Planck constant h [eV*s]
#define PHYS_c 299792458 // speed of light in vacuum [m/s]
#define PHYS_k 1.3806504e-23 // Boltzmann constant k [J/K]
#define PHYS_k_eV 8.617343e-5 // Boltzmann constant k [eV/K]
#define PHYS_k2 69.50356 // k/hc [1/(m*K)]
#define PHYS_flux 8.829e-21 // pi*e^2/mc^2, [cm^2/Å]

namespace phys {
	const double h = 6.62606896e-34; // Planck constant h [Js]
	const double h_eV = 4.13566733e-15; // Planck constant h [eV*s]
	const double c = 299792458.; // speed of light in vacuum [m/s]
	const double k = 1.3806504e-23; // Boltzmann constant k [J/K]
	const double k_eV = 8.617343e-5; // Boltzmann constant k [eV/K]
	const double k2 = 69.50356; // k/hc [1/(m*K)]
	const double flux = 8.829e-21; // pi*e^2/mc^2, [cm^2/Å]
}

#endif /* CONSTANTS_H_ */
