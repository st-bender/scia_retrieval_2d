/*
 * simple NO regression model interface
 *
 * Copyright (c) 2015-2017 Stefan Bender
 *
 * Initial version created on: 20.01.2015
 *      Author: Stefan Bender
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 *
 */
#ifndef NO_REGRESS_MODEL_HH_
#define NO_REGRESS_MODEL_HH_

#include <vector>

/* helper function to run the python model code with the right parameters */
std::vector<double> NO_regress_model_python(class Ausgewertete_Messung_Limb &aml,
		class Konfiguration &Konf, std::vector<double> alts, double lat);

#endif /* NO_REGRESS_MODEL_HH_ */
