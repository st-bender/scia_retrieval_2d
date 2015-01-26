/*
 * simple NO regression model
 *
 */

#include <ext/stdio_filebuf.h>
#include <iterator>
#include <istream>
#include <sstream>

#include "Ausgewertete_Messung_Limb.h"
#include "Konfiguration.h"

/* run the python model code with the right parameters */
std::vector<double> NO_regress_model_python(Ausgewertete_Messung_Limb &aml,
		Konfiguration &Konf, std::vector<double> alts, double lat)
{
	/* default script and data files */
	std::string script{"DATA/NO_regress_model.py"};
	std::string config{"DATA/NO_regress_model.json"};
	/* prepare the date string from the measurement date */
	std::string date{std::to_string(aml.m_Jahr) + "-"
			+ std::to_string(aml.m_Monat) + "-" + std::to_string(aml.m_Tag)};
	/* prepare the string with the altitudes
	 * for passing it to the python script */
	std::stringstream ss;
	std::string alt_string;
	ss << "\"";
	std::copy(alts.begin(), alts.end(), std::ostream_iterator<double>(ss, ","));
	ss << "\"";
	ss >> alt_string;
	ss.str(std::string());
	ss.clear();

	/* prepare the command line to run */
	std::string command{"python " + script + " " + config + " "
			+ date + " " + alt_string + " " + std::to_string(lat)};

	std::vector<double> NO_apriori;
	double NO_val;
	std::string line;
	int posix_handle{fileno(::popen(command.c_str(), "r"))};
	__gnu_cxx::stdio_filebuf<char> filebuf(posix_handle, std::ios::in);
	std::istream is(&filebuf);
	std::getline(is, line);
	ss << line;
	while (ss >> NO_val)
		NO_apriori.push_back(NO_val);
	//std::copy(NO_apriori.begin(), NO_apriori.end(),
	//		std::ostream_iterator<double>(std::cerr, ","));
	return NO_apriori;
}

