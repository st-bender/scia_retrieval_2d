/*
 * simple NO regression model interface
 *
 */
#ifndef NO_REGRESS_MODEL_HH_
#define NO_REGRESS_MODEL_HH_

#include <vector>

/* helper function to run the python model code with the right parameters */
std::vector<double> NO_regress_model_python(class Ausgewertete_Messung_Limb &aml,
		class Konfiguration &Konf, std::vector<double> alts, double lat);

#endif /* NO_REGRESS_MODEL_HH_ */
