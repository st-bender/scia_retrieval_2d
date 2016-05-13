/*****************************************
 * Messung.cpp
 *
 *  Created on: 27.01.2015
 *      Author: Stefan Bender
 *****************************************/

#include "Messung.h"
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <numeric>

#include <fstream>  //für Ausgabe
#include <iostream>//für Ausgabe
#include <sstream>
#include <cstdlib>  //für Ausgabe
#include <cstdio>   //Filekram
#include <iomanip>
#include <sys/stat.h>
#include "Konfiguration.h"
#include "Ausdrucke.h"
#include "Fit_Polynom.h"
#include "Glaetten.h"
#include "Speziesfenster.h"
#include "NO_emiss.h"
#include "Sonnenspektrum.h"
#include "Dateinamensteile_Bestimmen.h"

extern "C" {
#include "nrlmsise-00.h"
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
} //Lapackroutine zur Lösung eines linearen GLeichungssystems


using std::cout;
using std::string;
using std::stringstream;
using std::vector;

Messung::Messung(std::string filename) :
	m_Dateiname_L1C(filename)
{
	//initialisierung
	m_Zeilendichte = 0;
	m_Fehler_Zeilendichten = 0;
	m_Number_of_Wavelength = 0;
	m_Jahr = 0;
	m_Monat = 0;
	m_Tag = 0;
	m_Stunde = 0;
	m_Minute = 0;
	m_Sekunde = 0;
	m_Deklinationswinkel = 0;
	m_Sonnen_Longitude = 0;
	m_Latitude_Sat = 0;
	m_Longitude_Sat = 0;
	m_Hoehe_Sat = 0;
	m_Erdradius = 0;
	m_orbit_phase = 0.;
	// limb initialisierung
	total_number_density = 0.;
	m_Latitude_TP = 0;
	m_Longitude_TP = 0;
	m_Hoehe_TP = 0;
	m_TP_SZA = 0.;
	m_TP_rel_SAA = 0.;
	center_lat = 0.;
	center_lon = 0.;
	// nadir initialisierung
	// Herkunftsmerkmale
	m_Messung_ID = -1;
	// Geolokationen für Raytrace
	m_Latitude_Ground = 0;
	m_Longitude_Ground = 0;
	//statische Felder werden erstmal nicht 0 gesetzt
}
//========================================
//
//copyconstructor
//
//========================================
Messung::Messung(const Messung &rhs)
{
	*this = rhs;
}//copyconstructor ende

//========================================

//========================================
//
// Assignmentoperator Overload
//
//========================================
Messung &Messung::operator =(const Messung &rhs)
{
	//TODO das nochmal anpassen
	// Prevent self assignment. We say two Strings
	// are equal if their memory addresses are equal.
	if (this == &rhs)
		return *this;
	// Ergebnisse
	m_Zeilendichte = rhs.m_Zeilendichte;
	m_Fehler_Zeilendichten = rhs.m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	m_Deklinationswinkel = rhs.m_Deklinationswinkel;
	m_Sonnen_Longitude = rhs.m_Sonnen_Longitude;
	// Dateiname
	m_Dateiname_L1C = rhs.m_Dateiname_L1C;
	// Datum
	m_Jahr = rhs.m_Jahr;
	m_Monat = rhs.m_Monat;
	m_Tag = rhs.m_Tag;
	m_Stunde = rhs.m_Stunde;
	m_Minute = rhs.m_Minute;
	m_Sekunde = rhs.m_Sekunde;
	// Geolocation
	m_Latitude_Sat = rhs.m_Latitude_Sat;
	m_Longitude_Sat = rhs.m_Longitude_Sat;
	m_Hoehe_Sat = rhs.m_Hoehe_Sat;
	m_Erdradius = rhs.m_Erdradius;
	m_orbit_phase = rhs.m_orbit_phase;

	// limb only variables
	total_number_density = rhs.total_number_density;
	m_Latitude_TP = rhs.m_Latitude_TP;
	m_Longitude_TP = rhs.m_Longitude_TP;
	m_Hoehe_TP = rhs.m_Hoehe_TP;
	m_TP_SZA = rhs.m_TP_SZA;
	m_TP_rel_SAA = rhs.m_TP_rel_SAA;
	center_lat = rhs.center_lat;
	center_lon = rhs.center_lon;

	// nadir only variables
	// Herkunftsmerkmale
	m_Messung_ID = rhs.m_Messung_ID;
	//Geolocations
	m_Latitude_Ground = rhs.m_Latitude_Ground;
	m_Longitude_Ground = rhs.m_Longitude_Ground;

	m_Number_of_Wavelength = rhs.m_Number_of_Wavelength;
	// copy vectors
	m_Wellenlaengen = rhs.m_Wellenlaengen;
	m_Intensitaeten = rhs.m_Intensitaeten;
	m_Intensitaeten_relativer_Fehler = rhs.m_Intensitaeten_relativer_Fehler;
	m_Sonne = rhs.m_Sonne;
	m_Intensitaeten_durch_piF = rhs.m_Intensitaeten_durch_piF;
	m_Intensitaeten_durch_piF_Gamma = rhs.m_Intensitaeten_durch_piF_Gamma;
	m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand
		= rhs.m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand;

	// Return a reference to *this object.
	return *this;
}// Assignmentoperator Overload ende
//========================================
//========================================
//Methoden
//========================================

void Messung::Intensitaeten_normieren(vector<double> &Sonnen_Intensitaet)
{
	//Teiler wurde vorher interpoliert
	//todo prüfen
	// der Teiler ist das interpolierte Sonnenspektrum
	std::transform(m_Intensitaeten.begin(), m_Intensitaeten.end(),
			Sonnen_Intensitaet.begin(), m_Intensitaeten_durch_piF.begin(),
			std::divides<double>());
	m_Sonne = Sonnen_Intensitaet;
}
//========================================
//========================================
void Messung::Intensitaeten_durch_piF_Gamma_berechnen(Speziesfenster Spezfenst, double wl_gamma)
{

	//Auf dem ganzen Fenster...Verschwendung !!!!!...
	std::transform(m_Intensitaeten_durch_piF.begin(), m_Intensitaeten_durch_piF.end(),
			m_Intensitaeten_durch_piF_Gamma.begin(),
			std::bind2nd(std::divides<double>(), wl_gamma));
}
void Messung::Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen(Speziesfenster Spezfenst)
{

	//Auf dem ganzen Fenster...Verschwendung !!!!!

	// Wir berechnen den Gitterabstand nur einmal
	// Am besten gleich bei der Wellenlänge des Übergangs....
	// eigentlich reicht 0,11nm, falls es mal schneller gehn soll
	// die Gitterabstände sind aber über große Bereiche doch schon nicht linear

	// rough default to prevent it from being uninitialised.
	double Delta_WL = 0.11;
	// Nun alles damit multiplizieren....wie gesagt..das ist etwas langsam,
	// da es sich um nen konstanten Faktor handelt
	for (int i = 0; i < m_Number_of_Wavelength; i++) { //langsam, optimierbar
		//m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[i]=m_Intensitaeten
		//_durch_piF_Gamma[i]*Delta_WL;
		// Delta_Wl ist in nm gegeben...
		// dann muss beim Peakfit nicht in nm umgerechnet werden
		// wenn integriert wird
		if (i + 1 < m_Number_of_Wavelength)
			Delta_WL = m_Wellenlaengen[i + 1] - m_Wellenlaengen[i];

		m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[i]
			= m_Intensitaeten_durch_piF_Gamma[i] / Delta_WL;
		// glaub man muss dividieren
	}
}
//========================================
//========================================
void Messung::Deklinationswinkel_bestimmen()
{
	const double pi = M_PI;
	// Formel nach der englischen Wikipedia
	// https://en.wikipedia.org/wiki/Position_of_the_Sun#Declination_of_the_Sun_as_seen_from_Earth
	// theta = -23,44*cos(2*pi/365 * (N+10));
	// dieser Winkel ändert sich nicht sehr stark von Tag zu Tag
	int Monatstage[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	if (m_Jahr % 4 == 0 &&
			!(m_Jahr % 100 == 0 && m_Jahr % 400 != 0))
		Monatstage[1] = 29;
	// reicht auch auf Tagesgenauigkeit
	double Tage = -1.;
	for (int i = 0; i < (this->m_Monat - 1); i++) {
		Tage += Monatstage[i];
	}
	Tage += this->m_Tag;
	Tage += (double) this->m_Stunde / 24.0;

	this->m_Deklinationswinkel = -180.0 / pi * std::asin(
			0.39779 * std::cos(2. * pi / 365.24 * (Tage + 10.)
				+ 0.0334 * std::sin(2. * pi / 365.24 * (Tage - 2.))));
}// Deklinationswinkel_bestimmen() ende

//========================================
//========================================
// Funktionsstart  Sonnen_Longitude_bestimmen
void Messung::Sonnen_Longitude_bestimmen()
{
	// 12 Uhr Mittags (UTC) ist die Sonne bei Phi=0
	// (glaub Greenwich, oder zumindest in etwa) im Zenit
	double Stunden = 0.0;
	Stunden += this->m_Stunde;
	Stunden += (double) this->m_Minute / 60.0;

	this->m_Sonnen_Longitude = 180.0 - 360.0 * (Stunden / 24.0);
}
//ENDE Sonnen_Longitude_bestimmen
//========================================
//========================================
void Messung::calc_SunEarthDistance()
{
	const double a = 149598022.96; // semi-major axis
	const double e = 0.0167086342; // excentricity
	int Monatstage[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	if (m_Jahr % 4 == 0 &&
			!(m_Jahr % 100 == 0 && m_Jahr % 400 != 0))
		Monatstage[1] = 29;
	double Tage = -3.; // perihelion on 3 Jan.
	for (int i = 0; i < (this->m_Monat - 1); i++) {
		Tage += Monatstage[i];
	}
	Tage += this->m_Tag;
	Tage += (double) this->m_Stunde / 24.0;

	double theta = Tage / 365.24 * 2. * M_PI;
	m_SunEarthDistance = a * (1. - e * e) / (1. + e * std::cos(theta));
}

class lin_x {
	public:
	lin_x(double a0) : a(a0) {}
	double operator()(double x, double y) { return x + a*y; }
	private:
	double a;
};
double Messung::fit_NO_spec(NO_emiss &NO,
		std::vector<double> &x, std::vector<double> &y,
		double &rms_err)
{
	int l0 = sb_Get_Index(x.at(0)) + 1;

	double sum_gy = std::inner_product(y.begin(), y.end(),
			NO.spec_scia_res.begin() + l0, 0.0);
	double sum_gg = std::inner_product(NO.spec_scia_res.begin() + l0,
			NO.spec_scia_res.begin() + l0 + x.size(),
			NO.spec_scia_res.begin() + l0, 0.0);
	double A = sum_gy / sum_gg;

	std::vector<double> diffs;
	std::transform(y.begin(), y.end(), NO.spec_scia_res.begin() + l0,
			std::back_inserter(diffs), lin_x(-A));

	double err = std::inner_product(diffs.begin(), diffs.end(),
			diffs.begin(), 0.0) / (diffs.size() - 1) / sum_gg;
	rms_err = std::sqrt(err);

	return A;
}
double Messung::fit_NO_spec_weighted(NO_emiss &NO,
		std::vector<double> &x, std::vector<double> &y,
		std::vector<double> &yerr, double &rms_err)
{
	int l0 = sb_Get_Index(x.at(0)) + 1;

	/* prepare weighted vectors, 'yerr' should be something like sigma of y */
	std::vector<double> yw, NO_specw;
	std::transform(y.begin(), y.end(), yerr.begin(),
			std::back_inserter(yw), std::divides<double>());
	std::transform(NO.spec_scia_res.begin() + l0,
			NO.spec_scia_res.begin() + l0 + x.size(),
			yerr.begin(), std::back_inserter(NO_specw),
			std::divides<double>());

	/* In matrix terms we now calculate X^T*W*y = X'^T*y' and
	 * X^T*W*X = X'^T*X' with X' = w*X and y' = w*y, where
	 * w = sqrt(W) = 1 / yerr. */
	double sum_gwy = std::inner_product(yw.begin(), yw.end(),
			NO_specw.begin(), 0.0);
	double sum_gwg = std::inner_product(NO_specw.begin(),
			NO_specw.end(), NO_specw.begin(), 0.0);

	double A = sum_gwy / sum_gwg;

	rms_err = std::sqrt(1. / sum_gwg);

	return A;
}

class rayl {
	public:
	rayl(double f) : f_sol(f) {}
	double operator()(double x, double sol) {
		return f_sol * sigma_rayleigh(x) * sol;
	}
	private:
	double f_sol;
};
//========================================
/* Finds the Rayleigh fit factor for the Rayleigh background by fitting
 * (solar spectrum x Rayleigh cross section) to the limb spectrum
 * in the selected wavelength range.
 * Linearly interpolates large peak values, note that this alters
 * the values of m_Intensitaeten, which should be kept in mind whenever
 * m_Intensitaeten is used.
 */
double Messung::fit_rayleigh_and_interp_peaks(Sonnenspektrum &sol_spec,
		double wl_min, double wl_max, bool debug)
{
	// threshold for peak detection in the NO wavelength range
	// starting at 6*10^10 at 247 nm (NO(0, 2)) and increasing ~ lambda^4
	// because of Rayleigh scattering (see below)
	const double peak_threshold = 6.e10;
	int i_min = sb_Get_closest_index(wl_min);
	int i_max = sb_Get_closest_index(wl_max);
	std::vector<double> fit_spec, ones(i_max - i_min + 1, 1.);

	for (int i = i_min; i < i_max + 1; ++i) {
		double sol_i = sol_spec.m_Int_interpoliert.at(i);
		double rad_i = m_Intensitaeten.at(i);
		double wl = m_Wellenlaengen.at(i);
		// peak detection: unusual high radiance
		// threshold is 6*10^10 (see above) at 247 nm (NO(0, 2))
		// and scales ~ lambda^4 like Rayleigh scattering
		if (rad_i > peak_threshold * std::pow(wl / 247.0, 4)
				&& i > i_min + 2 && i < i_max - 2
				// make sure that the surrounding points are smaller, i.e.,
				// that we have a single large spike in the spectrum
				&& m_Intensitaeten.at(i - 1) < rad_i
				&& m_Intensitaeten.at(i + 1) < rad_i) {
			// exclude the previous, the current, and the next point.
			// That means we pop the last one and don't include the current
			// one, and interpolate the next point of the fit spectrum.
			if (!fit_spec.empty())
				fit_spec.pop_back();
			// interpolate three points of the peak linearly
			double y0 = m_Intensitaeten.at(i - 2);
			double yN = m_Intensitaeten.at(i + 2);
			double a = 0.25 * (yN - y0);
			for (int k = 0; k < 3; k++)
				m_Intensitaeten.at(i - 1 + k) = (k + 1)*a + y0;
			// done interpolating
			++i;
		} else {
			fit_spec.push_back(rad_i / (sigma_rayleigh(wl) * sol_i));
		}
	}
	ones.resize(fit_spec.size(), 1.);
	double f_sol_fit = fit_spectra(ones, fit_spec);
	if (debug == true)
		std::cout << "# solar fit factor = " << f_sol_fit << std::endl;
	if (f_sol_fit < 0.)
		f_sol_fit = 0.;
	return f_sol_fit;
}
//========================================
void Messung::slant_column_NO(NO_emiss &NO, string mache_Fit_Plots,
		Sonnenspektrum &sol_spec, int index,
		Speziesfenster &Spezfenst, std::string Arbeitsverzeichnis,
		bool debug)
{
	// I/(piFGamma)=integral(AMF n ds) mit AMF = s exp(-tau) ...aber zu der
	// Formel später nochmal zurück Das spätere Retrieval ermittelt dann die
	// Dichte n aus der rechten Seite

	//Zunächst Indizes der Wellenlaengen der Basisfenster bestimmen
	int i, j;
	int NO_NJ = NO.get_NJ();
	double wl;
	double f_sol_fit;
	double min_lambda_NO = 1000., max_lambda_NO = 0.;
	double gamma_threshold = 0.25 * NO.get_spec_scia_max();
	for (i = 0; i < NO_NJ; i++) {
		for (j = 0; j < 12; j++) {
			wl = NO.get_lambda_K(j, i);
			if (NO.get_gamma_j(j, i) > gamma_threshold) {
				if (wl > 0. && wl < min_lambda_NO) min_lambda_NO = wl;
				if (wl > max_lambda_NO) max_lambda_NO = wl;
			}
		}
	}
	// inner and outer baseline and peak window offset
	// the defaults (from M.L.) are base_offset_o = 3. and base_offset_i = 1.
	double base_offset_o = 1.5, base_offset_i = 0.3;
	int i_basewin_l_min = sb_Get_closest_index(min_lambda_NO - base_offset_o);
	int i_basewin_l_max = sb_Get_closest_index(min_lambda_NO - base_offset_i) - 1;
	int i_basewin_r_min = sb_Get_closest_index(max_lambda_NO + base_offset_i) + 1;
	int i_basewin_r_max = sb_Get_closest_index(max_lambda_NO + base_offset_o);
	int i_peakwin_min = sb_Get_closest_index(min_lambda_NO - base_offset_i);
	int i_peakwin_max = sb_Get_closest_index(max_lambda_NO + base_offset_i);
	// Speicherplatzbedarf für die Fenster ermitteln
	int base_l = (i_basewin_l_max - i_basewin_l_min + 1);
	int base_r = (i_basewin_r_max - i_basewin_r_min + 1);
	int N_fit_tot = i_basewin_r_max - i_basewin_l_min + 1;
	int N_base = base_l + base_r;
	int N_peak = i_peakwin_max - i_peakwin_min + 1;
	// Speicher anfordern
	std::vector<double> basewin_rad(N_base);
	std::vector<double> peakwin_wl(N_peak);
	std::vector<double> peakwin_rad(N_peak);
	std::vector<double> rad = m_Intensitaeten;
	std::vector<double> sol_rad = sol_spec.m_Int_interpoliert;

	/* prints the geolocation of the tangent point for later inspection */
	if (debug == true) {
		std::cout << "# TP: lat = " << m_Latitude_TP;
		std::cout << ", lon = " << m_Longitude_TP;
		std::cout << ", height = " << m_Hoehe_TP << std::endl;
		std::cout << "# orbit_phase = " << m_orbit_phase << std::endl;
		std::cout << "# NO band emission = " << NO.get_scia_band_emiss()
			<< ", NO rotational band emission = " << NO.get_band_emiss()
			<< std::endl;
	}

	f_sol_fit = fit_rayleigh_and_interp_peaks(sol_spec,
			min_lambda_NO - base_offset_o,
			max_lambda_NO + base_offset_o, debug);
	// reassign because m_Intensitaeten may have changed
	rad = m_Intensitaeten;

	// prepare baseline and rayleigh data
	std::vector<double> baseline_wl, baseline_rad, rayleigh_rad;
	std::vector<double> y, y_weights(N_fit_tot, 1.);
	std::copy(m_Wellenlaengen.begin() + i_basewin_l_min,
			m_Wellenlaengen.begin() + i_basewin_l_min + N_fit_tot,
			std::back_inserter(baseline_wl));
	std::transform(m_Wellenlaengen.begin() + i_basewin_l_min,
			m_Wellenlaengen.begin() + i_basewin_l_min + N_fit_tot,
			sol_rad.begin() + i_basewin_l_min,
			std::back_inserter(rayleigh_rad),
			rayl(f_sol_fit));
	std::transform(rad.begin() + i_basewin_l_min,
			rad.begin() + i_basewin_l_min + N_fit_tot,
			rayleigh_rad.begin(),
			std::back_inserter(y),
			std::minus<double>());

	// Basisfenster WL und I auffüllen
	std::copy(y.begin(), y.begin() + base_l,
			basewin_rad.begin());
	std::copy(y.begin() + base_l + N_peak, y.end(),
			basewin_rad.begin() + base_l);
	/* construct new baseline vectors by removing outliers
	 * This currently discards 20% (10% left and 10% right)
	 * of the baseline points. */
	std::vector<double> rad_sort(basewin_rad);
	std::sort(rad_sort.begin(), rad_sort.end());
	size_t ten_perc = rad_sort.size() / 10;
	double rad0 = rad_sort.at(ten_perc);
	double rad1 = rad_sort.at(rad_sort.size() - ten_perc - 1);

	for (i = 0; i < N_fit_tot; i++) {
		int idx = i_basewin_l_min + i;

		// prepare radiances and weights for the Whittaker smoother
		double radi = y.at(i);
		// exclude the peak window and outliers by zeroing the weights
		if ((idx >= i_peakwin_min && idx <= i_peakwin_max)
			|| radi < rad0 || radi > rad1)
			y_weights.at(i) = 0.;
	}
	// reset N_base
	N_base = std::accumulate(y_weights.begin(), y_weights.end(), 0);

	// replace the linear baseline by the Whittaker smoothed radiances
	// excluding the peak window and outliers as in the linear case.
	// the original (linear) baseline behaviour can be obtained by commenting
	// this line or by setting lambda (the 4th argument) to something large,
	// e.g. ~ 1.e9. (quick test showed that 3.e5 is quite close)
	double rms_err_base;
	baseline_rad = my_whittaker_smooth(y, y_weights, 2, 1.e4, rms_err_base);

	//Peakfenster WL und I auffüllen
	// lineare Funktion von Intensitäten des Peakfenster abziehen
	std::copy(m_Wellenlaengen.begin() + i_peakwin_min,
			m_Wellenlaengen.begin() + i_peakwin_min + N_peak,
			peakwin_wl.begin());
	std::transform(y.begin() + i_peakwin_min - i_basewin_l_min,
			y.begin() + i_peakwin_min - i_basewin_l_min + N_peak,
			baseline_rad.begin() + i_peakwin_min - i_basewin_l_min,
			peakwin_rad.begin(), std::minus<double>());
	double rms_err_peak;
	m_Zeilendichte = fit_NO_spec(NO, peakwin_wl, peakwin_rad,
			rms_err_peak);
	m_Fehler_Zeilendichten = std::sqrt(rms_err_base * rms_err_base
			+ rms_err_peak * rms_err_peak);

	if (mache_Fit_Plots == "ja" && Spezfenst.plot_fit) {
		// prepare data to plot
		std::vector<double> wavelengths, spec_wo_rayleigh = y, NO_fit;
		for (i = 0; i < base_l; i++) {
			wavelengths.push_back(m_Wellenlaengen.at(i_basewin_l_min + i));
			NO_fit.push_back(m_Zeilendichte *
					NO.get_spec_scia_res(i_basewin_l_min + i)
					+ baseline_rad.at(i));
		}
		for (size_t k = 0; k < peakwin_wl.size(); k++) {
			wavelengths.push_back(m_Wellenlaengen.at(i_peakwin_min + k));
			NO_fit.push_back(m_Zeilendichte *
					NO.get_spec_scia_res(i_peakwin_min + k)
					+ baseline_rad.at(i_peakwin_min - i_basewin_l_min + k));
		}
		for (i = 0; i < base_r; i++) {
			wavelengths.push_back(m_Wellenlaengen.at(i_basewin_r_min + i));
			NO_fit.push_back(m_Zeilendichte *
					NO.get_spec_scia_res(i_basewin_r_min + i)
					+ baseline_rad.at(i_basewin_r_min - i_basewin_l_min + i));
		}

		std::transform(spec_wo_rayleigh.begin(), spec_wo_rayleigh.end(),
				spec_wo_rayleigh.begin(),
				std::bind1st(std::multiplies<double>(), 1.e-9));
		std::transform(NO_fit.begin(), NO_fit.end(), NO_fit.begin(),
				std::bind1st(std::multiplies<double>(), 1.e-9));

		// plot the data to postscript files
		std::string s_OrbNum;
		std::stringstream buf;
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...
		// falls / im Namen ist das schlecht
		std::string Datnam = sb_basename(m_Dateiname_L1C);

		//TODO Pfad anpassen
		std::string plot_dir = Arbeitsverzeichnis + "/Plots";
		mkdir(plot_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		buf << Datnam.c_str() << "_" << Spezfenst.m_Spezies_Name.c_str()
			<< "_" << index << "_"
			<< std::setw(3) << std::setfill('0') << std::setprecision(0)
			<< std::fixed << m_Hoehe_TP << "km_"
			<< std::setw(3) << std::setfill('0') << std::setprecision(0)
			<< std::internal << std::fixed << m_Latitude_TP << "deg.pdf";
		if (m_Messung_ID != -1) {
		buf.str(string());
		buf << Datnam.c_str() << "_" << Spezfenst.m_Spezies_Name.c_str()
			<< "_" << m_Messung_ID << "_" << index << ".pdf";
		}
		std::string new_datnam(buf.str());
		std::string s1(plot_dir + "/" + new_datnam);
		// s1 ist der Volle Pfad der Datei... diesen wegspeichern,
		// um später die .pdf files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		// Orbitnummer ermitteln
		// die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			std::cerr << " kein .dat in Limbdateiname... "
					  << "Orbitnummer nicht findbar" << std::endl;
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		buf.str(std::string());
		buf << "Orbit " << s_OrbNum.c_str() << ", NO ("
			<< NO.get_vu() << ", " << NO.get_vl() << "), "
			<< std::resetiosflags(std::ios::fixed)
			<< std::setprecision(3) << m_Latitude_TP << u8"°N, "
			<< std::setprecision(3) << m_Longitude_TP << u8"°E, "
			<< std::setprecision(3) << m_Hoehe_TP << " km";
		if (m_Messung_ID != -1) {
		buf.str(string());
		buf << "Orbit " << s_OrbNum.c_str() << ", Nadir GP:"
			<< " lat: " << m_Latitude_Ground << " deg N,"
			<< " lon: " << m_Longitude_Ground << " deg E; Sat:"
			<< " lat: " << m_Latitude_Sat << " deg N,"
			<< " lon: " << m_Longitude_Sat << " deg E";
		}
		std::string s2(buf.str());

		Plot_2xy(Arbeitsverzeichnis.c_str(), s1.c_str(), s2.c_str(),
				 "wavelength [nm]",
				 "residual radiance [10^9 ph cm^{-2} s^{-1} nm^{-1}]",
				 wavelengths, spec_wo_rayleigh, wavelengths, NO_fit,
				 0, wavelengths.size() - 1,
				 m_Zeilendichte, m_Fehler_Zeilendichten);
	}

	if (debug == true) {
		std::cout << "# slant column = " << m_Zeilendichte;
		std::cout << ", error = " << m_Fehler_Zeilendichten << std::endl;
		std::cout << "# emissivity = "
			<< std::accumulate(peakwin_rad.begin(), peakwin_rad.end(), 0.)*0.11
			<< std::endl;
	}
}
//========================================
//Methoden ende
//========================================

////////////////////////////////////////////////////////////////////////////////
//
//Hilfsfunktionen
//
////////////////////////////////////////////////////////////////////////////////
int Messung::Get_Index(double WL)
{
	// Die Wellenlängen sind monoton steigend sortiert
	// Quicksearch algorithmus kann angewendet werden
	// Startelement-> 825/2 =412
	// Wellenlängen sind etwa 0,1 nm auseinander...also gebe Toleranz von 0,08nm
	// Maximale Schrittzahl ist ceil(log_2 N) also 10 Schritte
	bool gefunden = false;
	int unterer_Fensterindex = 0;
	int oberer_Fensterindex = 825;
	int startindex = 412;
	int aktueller_Index = startindex;

	if (abs(this->m_Wellenlaengen[oberer_Fensterindex] - WL) < 0.08) {
		return oberer_Fensterindex;
	}
	while (!gefunden) {
		//eventuell durch forschleife ersetzen,
		//weil Programm sich hier festhängen könnte
		// Wellenlänge gefunden
		if (abs(this->m_Wellenlaengen[aktueller_Index] - WL) < 0.08) {
			return aktueller_Index;// Schleifenabbruch und sofortige Rückgabe
		}
		if (m_Wellenlaengen[aktueller_Index] < WL) {
			unterer_Fensterindex = aktueller_Index;
		} else {
			oberer_Fensterindex = aktueller_Index;
		}//if m_WL<WL
		aktueller_Index = (oberer_Fensterindex + unterer_Fensterindex) / 2;
		//achtung hier wird abgerundet
	}//while
	return -1;
	// soweit sollte es eigentlich nicht kommen...
	// aber damit die Warnung verschwindet geben wir mal was zurück
}//ende int Messung::Get_Index(double WL)

int Messung::sb_Get_Index(double WL)
{
	vector<double>::iterator low;
	low = lower_bound(m_Wellenlaengen.begin(), m_Wellenlaengen.end(), WL);

	// catch edge cases
	if (low == m_Wellenlaengen.begin()) return 0;
	if (low == m_Wellenlaengen.end()) --low;

	return distance(m_Wellenlaengen.begin(), low) - 1;
}

int Messung::sb_Get_closest_index(double WL)
{
	int i = sb_Get_Index(WL);
	return (WL - m_Wellenlaengen[i]) < (m_Wellenlaengen[i + 1] - WL) ? i : i + 1;
}

void Messung::Fit_Linear(double *x, double *y, double &a0, double &a1,
		double &rms_err, int Anfangsindex, int Endindex)
{
	//fit der Funktion y=a0+a1x;
	//Bestimmung von a und b im Intervall zwischen Anfangs und endindex
	a0 = 0;
	a1 = 0;
	int i;
	// benötigt werden die Mittelwerte von x,y,x*y,und x^2 =====================
	double xsum = 0.;
	double ysum = 0.;
	double xysum = 0.;
	double xxsum = 0.;
	for (i = Anfangsindex; i <= Endindex; i++) {
		xsum += x[i];
		ysum += y[i];
		xysum += x[i] * y[i];
		xxsum += x[i] * x[i];
	}
	//Mittelwerte
	double N = Endindex - Anfangsindex + 1.;
	double x_m = xsum / N;
	double y_m = ysum / N;
	double xy_m = xysum / N;
	double xx_m = xxsum / N;
	//==========================================================================
	// Parameter b
	a1 = (xy_m - y_m * x_m) / (xx_m - x_m * x_m);
	//Parameter a
	a0 = y_m - a1 * x_m;

	double err = 0.;
	for (i = Anfangsindex; i <= Endindex; i++) {
		double diff = y[i] - a0 - a1 * x[i];
		err += diff * diff;
	}
	err /= N;
	rms_err = std::sqrt(err);
}
void Messung::Fit_Linear(vector<double> &x, vector<double> &y,
		double &a0, double &a1, double &rms_err,
		int Anfangsindex, int Endindex)
{
	//fit der Funktion y=a0+a1x;
	//Bestimmung von a und b im Intervall zwischen Anfangs und endindex
	a0 = 0.;
	a1 = 0.;
	int i;
	// benötigt werden die Mittelwerte von x,y,x*y,und x^2 =====================
	double xsum = 0.;
	double ysum = 0.;
	double xysum = 0.;
	double xxsum = 0.;
	for (i = Anfangsindex; i <= Endindex; i++) {
		xsum += x[i];
		ysum += y[i];
		xysum += x[i] * y[i];
		xxsum += x[i] * x[i];
	}
	//Mittelwerte
	double N = Endindex - Anfangsindex + 1.;
	double x_m = xsum / N;
	double y_m = ysum / N;
	double xy_m = xysum / N;
	double xx_m = xxsum / N;
	//==========================================================================
	// Parameter b
	a1 = (xy_m - y_m * x_m) / (xx_m - x_m * x_m);
	//Parameter a
	a0 = y_m - a1 * x_m;

	double err = 0.;
	for (i = Anfangsindex; i <= Endindex; i++) {
		double diff = y[i] - a0 - a1 * x[i];
		err += diff * diff;
	}
	err /= N;
	rms_err = std::sqrt(err);
}//Ende Fit_linear

void Messung::Fit_Polynom_4ten_Grades(double *x, double *y, double x0,
		double *Par_a0, double *Par_a1, double *Par_a2, double *Par_a3,
		double *Par_a4, int Anfangsindex, int Endindex)
{
	// Das geht auch analytisch, aber das ist eine laaaaaaaaaaaaaaaaange Formel,
	// deren Ableitung zwar trivial, aber Fehleranfällig ist(so vieeeeel zu
	// schreiben) Das Diagonalisieren des Linearen Gleichungssystems könnte man
	// auslagern,
	//sodass nur das Rückeinsetzen benutzt werden muss...  bei einer 5*6 Matrix
	//ist ist das aber vermutlich noch zu verschmerzen...wir werden sehn

	//Für jede Messung müssen 30 konstanten bestimmt werden, von denen aber
	//einige doppelt vorkommen

	double a0 = 0;
	double b0 = 0;
	double c0 = 0;
	double d0 = 0;
	double e0 = 0;
	double f0 = 0;
	/*b0*/              /*c0*/             /*d0*/             /*e0*/
	double e1 = 0;
	double f1 = 0;
	/*c0*/             /*d0*/             /*e0*/              /*e1*/
	double e2 = 0;
	double f2 = 0;
	/*d0*/             /*e0*/             /*e1*/              /*e2*/
	double e3 = 0;
	double f3 = 0;
	/*e0*/             /*e1*/            /*e2*/               /*e3*/
	double e4 = 0;
	double f4 = 0;
	// Nur 14 Parameter müssen echt bestimmt werden  //sieh 5x5 Matrix oben
	double h; //h wie hilfsvariable (soll nur ein Buchstabe sein)
	for (int i = Anfangsindex; i <= Endindex; i++) {
		h = x[i] - x0;
		a0 += 1;
		b0 += h;
		c0 += h * h;
		d0 += h * h * h;
		e0 += h * h * h * h;
		e1 += h * h * h * h * h;
		e2 += h * h * h * h * h * h;
		e3 += h * h * h * h * h * h * h;
		e4 += h * h * h * h * h * h * h * h;

		f0 += y[i];
		f1 += y[i] * h;
		f2 += y[i] * h * h;
		f3 += y[i] * h * h * h;
		f4 += y[i] * h * h * h * h;
	}
	// Zur Lösung des Gleichungssystems wird eine Lapack routine benutzt da
	// diese in Fortran geschrieben sind, muss die Transponierte Matrix
	// übergeben werden da die RHS nur aus einem Vektor besteht, gibt es keine
	// Verwirrung mit transponierten Matrix und Vektor aufbauen
	double LHS[25];
	double RHS[5];
	// LHS Matrix Spaltenweise eingeben
	LHS[0] = a0;
	LHS[5] = b0;
	LHS[10] = c0;
	LHS[15] = d0;
	LHS[20] = e0;
	LHS[1] = b0;
	LHS[6] = c0;
	LHS[11] = d0;
	LHS[16] = e0;
	LHS[21] = e1;
	LHS[2] = c0;
	LHS[7] = d0;
	LHS[12] = e0;
	LHS[17] = e1;
	LHS[22] = e2;
	LHS[3] = d0;
	LHS[8] = e0;
	LHS[13] = e1;
	LHS[18] = e2;
	LHS[23] = e3;
	LHS[4] = e0;
	LHS[9] = e1;
	LHS[14] = e2;
	LHS[19] = e3;
	LHS[24] = e4;

	RHS[0] = f0;
	RHS[1] = f1;
	RHS[2] = f2;
	RHS[3] = f3;
	RHS[4] = f4;
	// Restliche Vorbereitungen für Lapackroutine
	int N = 5;  //<-- Feldgröße Speed propto N^3 , LHS ist quadratisch, N ist Anzahl der Gitterpunkte
	int *IPIV;  //array mit der Pivotisierungsmatrix sollte so groß wie N sein, alle Elemente 0
	IPIV = new int[N];
	// ------ RHS oben definiert
	int NRHS = 1; //Spalten von RHS 1 nehmen, um keine c/Fortran Verwirrungen zu provozieren
	int LDA = N;
	int LDB = N;
	int INFO;
	//int Anzahl=N*N;  Anzahl sollte die Integer grenzen nicht überschreiten,
	//aber danbn sollte der Aufbau von LHS schon stören
	// AUFRUF   A ist LHS.transponiert und B ist RHS
	dgesv_(&N, &NRHS, LHS, &LDA, IPIV, RHS, &LDB, &INFO);
	// ENDE LU ZERLEGUNG
	delete[] IPIV;
	IPIV = 0;
	// Ergebnis steckt nun in RHS
	*Par_a0 = RHS[0];
	*Par_a1 = RHS[1];
	*Par_a2 = RHS[2];
	*Par_a3 = RHS[3];
	*Par_a4 = RHS[4];
	return;
}

void Messung::Fit_Peak_hyperbolic(double *x, double *y, double x0,
		double FWHM, double &A, int Anfangsindex, int Endindex)
{
	//Folgende Funktion ist fürs Integral über alles ordentlich auf 1 normiert
	//Spaltfunktion
	//cnorm           = 4.*PI*sqrt(2.) / FWHM**3      ! Normierung stimmt MPL
	//SlitFuncSPEC    = 1./( ( (.5*FWHM)**4 + X**4 ) * cnorm )
	//Im folgenden nenne ich SlitFuncSPEC==g
	//
	//Für einen Linearen Parameterfit der Funktion A*g gilt:
	// cih^2=sum(y-Ag)^2=sum(y^2-2Agy+A^2g^2)
	//dchi^2/dA=-2sum(gy)+2 A sum(g^2) das soll 0 sein
	//-> A=sum(gy)/sum(g^2)
	//
	// FWHM muss gegeben werden und wir werten die Funktion um den Mittelwert
	// x0 aus also statt X-> X-x0

	// In dieser Funktion wird die Fläche A der Spaltfunktion bestimmt, da die
	// Funktionwerte y=I/(piFGamma) sind so ist A dann die Säulendichte

	// Zahl der Messwertpaare
	// double, damit später keine Probleme beim weiterrechnen
	double sum_gy = 0.;
	double sum_gg = 0.;
	double g;

	for (int i = Anfangsindex; i <= Endindex; i++) {
		//g berechnen
		g = slit_func(FWHM, x0, x[i]);
		//eine Rechnung...nicht  Zeitkritisch
		// sum_gy erhöhen
		sum_gy += g * y[i];
		// sum_gg erhöhen
		sum_gg += g * g;
	}
	A = sum_gy / sum_gg;
}

void Messung::Fit_Peak_hyperbolic(vector<double> &x, vector<double> &y,
		double x0, double FWHM, double &A, int Anfangsindex, int Endindex)
{
	//Folgende Funktion ist fürs Integral über alles ordentlich auf 1 normiert
	//Spaltfunktion
	//cnorm           = 4.*PI*sqrt(2.) / FWHM**3      ! Normierung stimmt MPL
	//SlitFuncSPEC    = 1./( ( (.5*FWHM)**4 + X**4 ) * cnorm )
	//Im folgenden nenne ich SlitFuncSPEC==g
	//
	//Für einen Linearen Parameterfit der Funktion A*g gilt:
	// cih^2=sum(y-Ag)^2=sum(y^2-2Agy+A^2g^2)
	//dchi^2/dA=-2sum(gy)+2 A sum(g^2) das soll 0 sein
	//-> A=sum(gy)/sum(g^2)
	//
	// FWHM muss gegeben werden und wir werten die Funktion um den Mittelwert
	// x0 aus also statt X-> X-x0

	// In dieser Funktion wird die Fläche A der Spaltfunktion bestimmt, da die
	// Funktionwerte y=I/(piFGamma) sind so ist A dann die Säulendichte

	// Zahl der Messwertpaare
	// double, damit später keine Probleme beim weiterrechnen
	double sum_gy = 0.;
	double sum_gg = 0.;
	double g;

	for (int i = Anfangsindex; i <= Endindex; i++) {
		//g berechnen
		g = slit_func(FWHM, x0, x[i]);
		//eine Rechnung...nicht  Zeitkritisch
		// sum_gy erhöhen
		sum_gy += g * y[i];
		// sum_gg erhöhen
		sum_gg += g * g;
	}
	A = sum_gy / sum_gg;
}

double Messung::Evaluate_Error_primitive(double *x, double *y, double a0,
		double a1, double A, double FWHM, double x0, int Anfangsindex,
		int Endindex)
{
	/***************************************************************************
	Wie der Name schon sagt, ist dies eine eher einfache Berechnung des Fehlers.
	Summe der Quadratischen Abweichungen-> Chi^2 hmm nicht gut... aber als
	Wichtungsfaktor noch akzeptabel
	 **************************************************************************/
	double Error = 0;
	//double y_quadrat=0;
	for (int i = Anfangsindex; i < Endindex + 1; i++) {
		//Funktionswert Bestimmen
		double Basis = a0 + a1 * x[i];
		double Peak = A * slit_func(FWHM, x0, x[i]);
		double Funktionswert = Peak + Basis;
		// Quadratische Abweichung des Funktionswerts zum Messwert Bestimmen
		// und aufaddieren
		Error += (Funktionswert - y[i]) * (Funktionswert - y[i]);
		//y_quadrat+=y[i]*y[i];
	}
	Error /= (Endindex - Anfangsindex + 1);
	Error = sqrt(Error);
	//Das ist nach Numerical Recipes der Fehlerbalken der Messpunkte Es ist
	//vermutlich anschaulicher diesen Fehler noch durch den Mittelwert zu
	//teilen;
	//y_quadrat/=(Endindex-Anfangsindex+1);
	//double quot=Error/sqrt(y_quadrat);
	//cout<<quot<<"\n";
	return Error;
}
double Messung::Evaluate_Error_primitive(vector<double> &x,
		vector<double> &y, double a0, double a1, double A, double FWHM,
		double x0, int Anfangsindex, int Endindex)
{
	/***************************************************************************
	Wie der Name schon sagt, ist dies eine eher einfache Berechnung des Fehlers.
	Summe der Quadratischen Abweichungen-> Chi^2 hmm nicht gut... aber als
	Wichtungsfaktor noch akzeptabel
	 **************************************************************************/
	double Error = 0;
	//double y_quadrat=0;
	for (int i = Anfangsindex; i < Endindex + 1; i++) {
		//Funktionswert Bestimmen
		double Basis = a0 + a1 * x[i];
		double Peak = A * slit_func(FWHM, x0, x[i]);
		double Funktionswert = Peak + Basis;
		// Quadratische Abweichung des Funktionswerts zum Messwert Bestimmen
		// und aufaddieren
		Error += (Funktionswert - y[i]) * (Funktionswert - y[i]);
		//y_quadrat+=y[i]*y[i];
	}
	Error /= (Endindex - Anfangsindex + 1);
	Error = sqrt(Error);
	//Das ist nach Numerical Recipes der Fehlerbalken der Messpunkte Es ist
	//vermutlich anschaulicher diesen Fehler noch durch den Mittelwert zu
	//teilen;
	//y_quadrat/=(Endindex-Anfangsindex+1);
	//double quot=Error/sqrt(y_quadrat);
	//cout<<quot<<"\n";
	return Error;
}
////////////////////////////////////////////////////////////////////////////////
//
//Hilfsfunktionen ENDE
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen
//
////////////////////////////////////////////////////////////////////////////////
void Messung::Ausgabe_in_Datei(string Dateiname)
{
	//TODO hier kann sich später auch noch was verändern
	/************************************************************
	Diese Funktion gibt den Inhalt des Objekts in eine Datei aus.
	Damit kann überprüft werden:
	a) ob das einlesen aller parameter ordentlich geklappt hat
	b) ob die Unterfunktionen das richtige errechnet haben
	c) Die Felder können geplottet werden, zusammen mit den Fitfunktionen->
	   Überprüfung, ob Fit sinnvoll (z.b. mit Matlab oder Gnuplot)
	************************************************************/
	//Formatierte Ausgabe
	FILE *outfile;
	//Datei öffnen
	outfile = fopen(Dateiname.c_str(), "w");
	//checken, ob Datei auch offen fehlt...aber ok Funktion wird eh beim
	//debuggen eingesetzt...da kriegt man das schon raus..hoffentlich Datei
	//schreiben
	///////////////////////////////////////////////////////////
	//zunächst die randdaten in den header
	//fprintf(outfile,"blblblbllblblb");
	fprintf(outfile, "%20s %1.5E\n", "m_Zeilendichte[cm^-2]: ", m_Zeilendichte);
	fprintf(outfile, "%20s %1.5E\n", "Sigma für alle Messwerte[cm^-2]: ",
			m_Fehler_Zeilendichten);
	fprintf(outfile, "%20s %12s\n", "m_Dateiname_L1C: ",
			m_Dateiname_L1C.c_str()); //hoffentlich passt das auch
	fprintf(outfile, "%20s %5i\n", "m_Number_of_Wavelength: ",
			m_Number_of_Wavelength);
	fprintf(outfile, "%20s %5i\n", "m_Jahr: ", m_Jahr);
	fprintf(outfile, "%20s %5i\n", "m_Monat: ", m_Monat);
	fprintf(outfile, "%20s %5i\n", "m_Tag: ", m_Tag);
	fprintf(outfile, "%20s %5i\n", "m_Stunde: ", m_Stunde);
	fprintf(outfile, "%20s %5i\n", "m_Minute: ", m_Minute);
	fprintf(outfile, "%20s %1.5E\n", "m_Deklinationswinkel: ",
			m_Deklinationswinkel);
	fprintf(outfile, "%20s %1.5E\n", "m_Sonnen_Longitude: ",
			m_Sonnen_Longitude);
	//Geolokationen
	fprintf(outfile, "%20s %1.5E\n", "m_Latitude_Sat: ", m_Latitude_Sat);
	fprintf(outfile, "%20s %1.5E\n", "m_Longitude_Sat: ", m_Longitude_Sat);
	fprintf(outfile, "%20s %1.5E\n", "m_Hoehe_Sat: ", m_Hoehe_Sat);
	fprintf(outfile, "%20s %1.5E\n", "m_Erdradius: ", m_Erdradius);
	//Überschrift
	fprintf(outfile, "%12s %12s %12s %12s %12s\n",
			"m_Wellenlaengen", "m_Intensitaeten", "m_Intensitaeten_durch_piF",
			"m_Intensitaeten_durch_piF_Gamma",
			"m_Intensitaeten_durch_piF_Gamma*ds");
	//for(int i=0;i<m_Number_of_Wavelength;i++)
	for (int i = 0; i < m_Number_of_Wavelength; i++) {
		//die letzte Zeile der Datei ist leer, da \n in der Vorletzten steht
		fprintf(outfile, "%1.5E %1.5E %1.5E %1.5E %1.5E\n",
				m_Wellenlaengen[i], m_Intensitaeten[i],
				m_Intensitaeten_durch_piF[i],
				m_Intensitaeten_durch_piF_Gamma[i],
				m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[i]);
	}
	///////////////////////////////////////////////////////////
	// Datei schließen
	fclose(outfile);
}//Ausgabe_in_Datei ENDE
////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen ENDE
//
////////////////////////////////////////////////////////////////////////////////

