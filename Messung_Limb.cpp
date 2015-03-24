/*****************************************
 * Messung_Limb.cpp
 *
 *  Created on: 20.04.2010
 *      Author: martin
 *****************************************/

#include "Messung_Limb.h"
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

Messung_Limb::Messung_Limb(std::string filename) :
	Messung(filename)
{
	//initialisierung
	total_number_density = 0.;
	m_Latitude_TP = 0;
	m_Longitude_TP = 0;
	m_Hoehe_TP = 0;
	m_TP_SZA = 0.;
	m_TP_rel_SAA = 0.;
	center_lat = 0.;
	center_lon = 0.;
	//statische Felder werden erstmal nicht 0 gesetzt
}

//========================================
//
// Assignmentoperator Overload
//
//========================================
Messung_Limb &Messung_Limb::operator =(const Messung_Limb &rhs)
{
	//TODO das nochmal anpassen
	// Prevent self assignment. We say two Strings
	// are equal if their memory addresses are equal.
	if (this == &rhs)
		return *this;
	Messung::operator=(rhs);
	total_number_density = rhs.total_number_density;
	m_Latitude_TP = rhs.m_Latitude_TP;
	m_Longitude_TP = rhs.m_Longitude_TP;
	m_Hoehe_TP = rhs.m_Hoehe_TP;
	m_TP_SZA = rhs.m_TP_SZA;
	m_TP_rel_SAA = rhs.m_TP_rel_SAA;
	center_lat = rhs.center_lat;
	center_lon = rhs.center_lon;
	m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand
		= rhs.m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand;

	// Return a reference to *this object.
	return *this;
}// Assignmentoperator Overload ende
//========================================
//========================================
//Methoden
//========================================
//========================================
void Messung_Limb::Zeilendichte_Bestimmen(Speziesfenster &Spezfenst, int Index,
		string Arbeitsverzeichnis, string mache_Fit_Plots)
{
	//kurz:
	//Diese Fitroutine ermittelt die Fläche von I/(piF*Gamma) bei der
	//Übergangswellenlänge Da nur wenige Punkte zu einer Linie beitragen,
	//entspricht der Linie der Auflösungsfunktion
	//Diese ist ein Hyperboloide 4ter Ordnung

	//lang:
	// Die Umgebung jeder Linie wird in 3 Bereiche unterteilt. 2 sogenannte
	// Basisfenster links und rechts von der Linie und das Peakfenster selbst,
	// welches in der Umsetzung auch Bereiche aus den beiden Basisfenstern
	// beinhalten darf.  Bis Hierhin haben wir I/(piFGamma) berechnet. Nun
	// wollen wir die "normierte" Intensität für einen Peak berechnen Da der
	// Peak gerade im Kanal zumeist auf einem relativ großen Untergrundsignal
	// sitzt, wird dieses zunächst abgezogen.  Dafür wird im Bereich um die
	// Linie herum, wo keine signifikanten anderen Linien liegen(Basisfenster)
	// ein linearer Fit der Grundlinie durchgeführt, welcher vom der
	// "normierten" Intensität abgezogen wird.  Danach wird im Bereich des
	// Peakfensters die Fläche des Peaks bestimmt.  Da die Linie zumeist nur
	// aus 3 Punkten besteht,  hat das Profil im wesentlichen die Form der
	// Auflösungsfunktion, welche eine hyperboloide 4ten Grades ist. Diese hat
	// im wesentlichen 3 Parameter: Wellenlänge des Peaks, Breite und Fläche.
	// Die robusteste Variante ist es, die Wellenlänge und die Halbwertsbreite
	// der Linie vorzugeben und nur die Fläche anzufitten.  Diese Variante kann
	// durch einen linearen Parameterfit erreicht werden und ist zum einen
	// schnell und was hier viel wichtiger ist robust. Andere nichtlineare
	// Verfahren sind entweder nicht robust(Gauss Newton, insbesondere bei dem
	// schlechten Signal/Noise Verhältnis) oder deutlich langsamer (simulated
	// annealing bzw. ist das aufwändiger dies noch schnell zu implementieren...
	// das wäre Schritt 2)
	// Die ermittelte Fläche ist dann unsere gesuchte Säulenzeilendichte.

	// I/(piFGamma)=integral(AMF n ds) mit AMF = s exp(-tau) ...aber zu der
	// Formel später nochmal zurück Das spätere Retrieval ermittelt dann die
	// Dichte n aus der rechten Seite

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = sb_Get_closest_index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = sb_Get_closest_index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = sb_Get_closest_index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = sb_Get_closest_index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
	int Index_Peakfenster_min = sb_Get_closest_index(Spezfenst.m_Peakfenster_WLmin[Index]);
	int Index_Peakfenster_max = sb_Get_closest_index(Spezfenst.m_Peakfenster_WLmax[Index]);
	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int N_Basis = Bas_l + Bas_r;
	int N_Peak = Index_Peakfenster_max - Index_Peakfenster_min + 1;
	// Speicher anfordern
	vector<double> Basisfenster_WL(N_Basis);
	vector<double> Basisfenster_Intensitaet(N_Basis);
	vector<double> Peakfenster_WL(N_Peak);
	vector<double> Peakfenster_Intensitaet(N_Peak);
	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i] = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] = this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_rechts_min + i];
	}
	//Peakfenster WL und I auffüllen
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_WL[i] = m_Wellenlaengen[Index_Peakfenster_min + i];
		Peakfenster_Intensitaet[i] =
			m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Peakfenster_min + i];
	}
	// linearen Fit des Basisfensters durchführen
	// Proto: Fit_Linear(double* x,double* y, double& a0, double& a1,int
	// Anfangsindex, int Endindex)
	double a0, a1, rms_err_base;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, rms_err_base,
			0, N_Basis - 1);
	// lineare Funktion von Intensitäten des Peakfenster abziehen
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_Intensitaet[i] -= a0 + a1 * Peakfenster_WL[i];
	}
	// Hyperboloiden an Peakfenster anfitten
	// Proto: Fit_Peak_hyperbolic(double* x,double* y,double x0, double FWHM,
	// double& A, int Anfangsindex, int Endindex)
	double Flaeche;
	Fit_Peak_hyperbolic(Peakfenster_WL, Peakfenster_Intensitaet,
						Spezfenst.m_Wellenlaengen[Index],
						Spezfenst.m_FWHM, Flaeche, 0, N_Peak - 1);
	// Hier Wellenlängen in nm verwendet..das hebt sich mit dem Gitterabstand
	// raus
	//Fehler des Fits bestimmen... da Peakfenster_Intensitaet nicht mehr
	//gebraucht wird die Basislinie für die Fehlerberechnung wieder aufaddiert
	m_Zeilendichte = Flaeche;
	//cout<<m_Zeilendichte<<"\n";
	//cout<<Spezfenst.m_Liniendaten[Index].m_Gamma<<"\n";
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_Intensitaet[i] += a0 + a1 * Peakfenster_WL[i];
	}
	// Funktion double Messung_Limb::Evaluate_Error_primitive(double* x,
	// double* y, double a0,double a1, double A, double FWHM, double x0,
	// int Anfangsindex, int Endindex)
	m_Fehler_Zeilendichten = Evaluate_Error_primitive(Peakfenster_WL,
			Peakfenster_Intensitaet, a0, a1, m_Zeilendichte, Spezfenst.m_FWHM,
			Spezfenst.m_Wellenlaengen[Index], 0, N_Peak - 1);

	////////////////////////////////////////////////////////////////////////////
	// Hier kann man zur Testzwecken noch einen Plot machen  ///////////////////
	if (mache_Fit_Plots == "ja" && Spezfenst.plot_fit) {
		//TODO das als Funktion implementieren
		vector<double> Funktion(N_Peak);

		for (int i = 0; i < N_Peak; i++) {
			double Basis = a0 + a1 * Peakfenster_WL[i];
			double Peak = m_Zeilendichte *
				slit_func(Spezfenst.m_FWHM,
						Spezfenst.m_Wellenlaengen[Index], Peakfenster_WL[i]);
			//cout<<Peak<<"\n";
			//cout<<m_Zeilendichte<<"\n";
			Funktion[i] = Peak + Basis;
		}

		string s_OrbNum;
		stringstream buf;
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...
		// falls / im Namen ist das schlecht
		string Datnam = sb_basename(m_Dateiname_L1C);

		//TODO Pfad anpassen
		string plot_dir = Arbeitsverzeichnis + "/Plots";
		mkdir(plot_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		buf << Datnam.c_str() << "_" << Spezfenst.m_Spezies_Name.c_str()
			<< "_" << Index << "_" << m_Hoehe_TP << "km.pdf";
		string new_datnam(buf.str());
		string s1(plot_dir + "/" + new_datnam);
		//s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern, um
		//später die .pdf files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal wie die Datei heißt, die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Limbdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		buf.str(string());
		buf << "Orbit " << s_OrbNum.c_str() << ", Limb TP:"
			<< " Lat: " << m_Latitude_TP << " deg,"
			<< " Lon: " << m_Longitude_TP << " deg,"
			<< " Hoehe: " << m_Hoehe_TP << " km.";
		string s2(buf.str());
		//cout<<s1<<"\n";
		//int Plot_2xy(string Dateiname,string title, string xlabel,
		//string ylabel,double* x1,double*y1, double* x2,double* y2,
		//int Startindex,int Endindex);
		//Plot_2xy(s1.c_str(),s1.substr(s1.size()-50,50).c_str(),
		//"$\\lambda$ in nm","$\\frac{I}{\\piF\\gamma}$",Peakfenster_WL,
		//Peakfenster_Intensitaet,Peakfenster_WL,Funktion,0,N_Peak-1);
		//-> Fit geht
		Plot_2xy(Arbeitsverzeichnis.c_str(), s1.c_str(), s2.c_str(),
				 "Wellenlaenge in nm",
				 "Schraege Saeule bei Peakposition in cm^{-2}/nm",
				 Peakfenster_WL, Peakfenster_Intensitaet, Peakfenster_WL,
				 Funktion, 0, N_Peak - 1,
				 m_Zeilendichte, m_Fehler_Zeilendichten);
	}
	// Ende Plot ///////////////////////////////////////////////////////////////
}//int Zeilendichte_Bestimmen() ende
//========================================
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
void Messung_Limb::slant_column_NO(NO_emiss &NO, string mache_Fit_Plots,
		Sonnenspektrum &sol_spec, int index,
		Speziesfenster &Spezfenst, std::string Arbeitsverzeichnis,
		bool debug)
{
	// I/(piFGamma)=integral(AMF n ds) mit AMF = s exp(-tau) ...aber zu der
	// Formel später nochmal zurück Das spätere Retrieval ermittelt dann die
	// Dichte n aus der rechten Seite

	// threshold for peak detection in the NO wavelength range
	// starting at 6*10^10 at 247 nm (NO(0, 2)) and increasing ~ lambda^4
	// because of Rayleigh scattering (see below)
	const double peak_threshold = 6.e10;

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
	std::vector<double> fit_spec, ones(N_base + N_peak, 1.);

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

	for (i = 0; i < N_base + N_peak; i++) {
		int idx = i_basewin_l_min + i;
		double sol_i = sol_rad.at(idx);
		double rad_i = rad.at(idx);
		wl = m_Wellenlaengen.at(idx);
		// peak detection: unusual high radiance
		// threshold is 6*10^10 (see above) at 247 nm (NO(0, 2))
		// and scales ~ lambda^4 like Rayleigh scattering
		if (rad_i > peak_threshold * std::pow(wl / 247.0, 4)
				&& i > 2 && i < N_base + N_peak - 2
				// make sure that the surrounding points are smaller, i.e.,
				// that we have a single large spike in the spectrum
				&& rad.at(idx - 1) < rad_i
				&& rad.at(idx + 1) < rad_i) {
			// exclude the previous, the current, and the next point.
			// That means we pop the last one and don't include the current
			// one, and interpolate the next point of the fit spectrum.
			if (!fit_spec.empty()) fit_spec.pop_back();
			// interpolate three points of the peak linearly
			double y0 = rad.at(idx - 2);
			double yN = rad.at(idx + 2);
			double a = 0.25 * (yN - y0);
			for (int k = 0; k < 3; k++)
				rad.at(idx - 1 + k) = (k + 1)*a + y0;
			// done interpolating
			i++;
		} else
			fit_spec.push_back(rad_i / (sigma_rayleigh(wl) * sol_i));
	}
	ones.resize(fit_spec.size());
	f_sol_fit = fit_spectra(ones, fit_spec);
	if (debug == true)
		std::cout << "# solar fit factor = " << f_sol_fit << std::endl;
	if (f_sol_fit < 0.) f_sol_fit = 0.;

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
	size_t offset = rad_sort.size() / 10;
	double rad0 = rad_sort.at(offset);
	double rad1 = rad_sort.at(rad_sort.size() - offset - 1);

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
	double rms_err_peak, rms_err_tot;
	m_Zeilendichte = fit_NO_spec(NO, peakwin_wl, peakwin_rad,
			rms_err_peak);
	rms_err_tot = std::sqrt((N_base * rms_err_base * rms_err_base
		+ N_peak * rms_err_peak * rms_err_peak) / (N_base + N_peak));
	m_Fehler_Zeilendichten = rms_err_tot / NO.get_spec_scia_max();

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
			<< std::fixed << m_Latitude_TP << "deg.pdf";
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
		buf << "Orbit " << s_OrbNum.c_str() << ", "
			<< NO.get_vu() << NO.get_vl() << ", "
			<< std::resetiosflags(std::ios::fixed)
			<< " Lat: " << std::setprecision(3) << m_Latitude_TP << " deg,"
			<< " Lon: " << std::setprecision(3) << m_Longitude_TP << " deg,"
			<< " Alt: " << std::setprecision(3) << m_Hoehe_TP << " km.";
		std::string s2(buf.str());

		Plot_2xy(Arbeitsverzeichnis.c_str(), s1.c_str(), s2.c_str(),
				 "wavelength [nm]",
				 "residual radiance [10^9 ph/cm^2/s/nm]",
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

////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Saeulendichte_Bestimmen_MgI285nm
////////////////////////////////////////////////////////////////////////////////
void Messung_Limb::Saeulendichte_Bestimmen_MgI285nm(Speziesfenster &Spezfenst,
		int Index, string Arbeitsverzeichnis, string mache_Fit_Plots,
		double *mean_10_20)
{
	//kurz:
	//alternative Fitroutine für die Bestimmung  der MgI 285.21275 nm Linie
	//lang:
	// Das Signal zu Rausch Verhältnis für die Limbmessungen scheint ziemlich
	// mies zu sein, sodass die Befürchtung besteht, dass man da Rauschen als
	// Messwerte interpretiert
	// Ein Erster Rettungsversuch ist es, statt den Qoutienten I/piF direkt aus
	// den Messwerten zu bilden, die breiten Peaks im Limb und Sonnenspektrum
	// auszunutzen
	//
	// I/(Fgamma * oder / Integrations wegstück in nm) ist schon vorhanden
	// in dem entsprechendem Fenster wird die Basislinie diese Quotienten der
	// Messwerte gebildet Limb und Sonnenspektrum werden einzeln als Polynome
	// angefittet
	// Da das Minimum eines nicht geraden Polynoms nicht bei der vorgegebenen
	// Wellenlänge liegt wird das Limbspektrum so geshifted, dass die Minima
	// übereinander liegen
	// Der Quotient beider Polynome wird gebildet
	// Der konstante Faktor aus gamma und infinitesimalen Wegstück(ca. 0.11
	// also WL(2)-WL(1)) wird anmultipliziert, sodass Quotient aus Messwerten
	// und Quotient aus Polynomen die gleiche Größenordnung haben
	// Die Basislinie wird nun von beiden abgezogen
	// Wenn überhaupt, sollte man erst jetzt den Quotienten so shiften, dass
	// das Minimum bei der Wellenlänge des Übergangs liegt.  Abschließend
	// werden für jede Messung mehrere Plots erstellt

	double *Basisfenster_WL;
	double *Basisfenster_Intensitaet;
	double *Vollfenster_WL;
	double *Vollfenster_Limb;
	double *Vollfenster_Sonne;

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
	int Index_Uebergangs_Wellenlaenge = Get_Index(Spezfenst.m_Wellenlaengen[Index]) - Index_Basisfenster_links_min;

	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int N_Basis = Bas_l + Bas_r;
	//int N_Peak=Index_Peakfenster_max-Index_Peakfenster_min+1;
	int N_Vollfenster = Index_Basisfenster_rechts_max - Index_Basisfenster_links_min + 1;
	// Speicher anfordern
	Basisfenster_WL = new double[N_Basis];
	Basisfenster_Intensitaet = new double[N_Basis];
	Vollfenster_WL = new double[N_Vollfenster];
	Vollfenster_Limb = new double[N_Vollfenster];
	Vollfenster_Sonne = new double[N_Vollfenster];

	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i] = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] = this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_rechts_min + i];
	}

	//Vollfenster Limb und Sonne auffüllen
	for (int i = 0; i < N_Vollfenster; i++) {
		Vollfenster_WL[i] = m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Vollfenster_Limb[i] = m_Intensitaeten[Index_Basisfenster_links_min + i];
		Vollfenster_Sonne[i] = m_Sonne[Index_Basisfenster_links_min + i];
	}

	double *Vollfenster_Limb_mittlere_atmo; //zwischen 40 und 60km
	Vollfenster_Limb_mittlere_atmo = new double[N_Vollfenster];
	for (int i = 0; i < N_Vollfenster; i++) {
		Vollfenster_Limb_mittlere_atmo[i] = mean_10_20[Index_Basisfenster_links_min + i];
	}
	//Minima vergleichen um Linie herum (5 nachbarpunkte)
	// und shift durchführen
	// Suche mit Glattem Fenster, Verschiebung im unglatten
	int Suchbereich = 5;
	int min_Index_Limb = Index_Uebergangs_Wellenlaenge - Suchbereich;
	int min_Index_Sun = Index_Uebergangs_Wellenlaenge - Suchbereich;
	for (int i = Index_Uebergangs_Wellenlaenge - Suchbereich;
			i <= Index_Uebergangs_Wellenlaenge + Suchbereich; i++) {
		if (Vollfenster_Limb_mittlere_atmo[i] < Vollfenster_Limb_mittlere_atmo[min_Index_Limb]) {
			min_Index_Limb = i;
		}
		if (Vollfenster_Sonne[i] < Vollfenster_Sonne[min_Index_Sun]) {
			min_Index_Sun = i;
		}
	}
	int Verschiebung = min_Index_Sun - min_Index_Limb;
	if ((Verschiebung >= 1) || (Verschiebung <= -1)) {
		cout << "Verschiebung: " << Verschiebung << "\n";
		cout << "Minima auseinander\n";
	}
	if ((Verschiebung >= 3) || (Verschiebung <= -3)) {
		cout << "Verschiebung: " << Verschiebung << "\n";
		cout << "Minima weit auseinander\n";
	}
	double *vor_Verschiebung;
	vor_Verschiebung = new double[N_Vollfenster];
	for (int i = 0; i < N_Vollfenster; i++) {
		vor_Verschiebung[i] = Vollfenster_Limb[i];
	}
	for (int i = 0; i < N_Vollfenster; i++) {
		if (!(((i + Verschiebung) < 0)
			|| (i + Verschiebung >= N_Vollfenster))) {
			Vollfenster_Limb[i + Verschiebung] = vor_Verschiebung[i];
		}
	}
	delete[] vor_Verschiebung;
	delete[] Vollfenster_Limb_mittlere_atmo;
	// Spektrum glätten
	int smooth_Nachbarn = 0; //2;
	int smooth_Iterationen = 0; //4;
	smooth_data(N_Vollfenster, Vollfenster_Limb, smooth_Nachbarn, smooth_Iterationen);
	smooth_data(N_Vollfenster, Vollfenster_Sonne, smooth_Nachbarn, smooth_Iterationen);

	// linearen Fit des Basisfensters durchführen
	// Proto:
	// Fit_Linear(double* x,double* y, double& a0, double& a1,
	//   int Anfangsindex, int Endindex)
	double a0, a1, rms_err_base;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, rms_err_base, 0,
			N_Basis - 1);

	int Polynomgrad = 4;
	//Get_Index(Spezfenst.m_Wellenlaengen[Index]) darf hier nicht benutzt werden
	int Fitpunktezahl_links = 4;
	int Polyfit_Startindex = Index_Uebergangs_Wellenlaenge - Fitpunktezahl_links;
	int Polyfit_Endindex = Index_Uebergangs_Wellenlaenge + Fitpunktezahl_links;
	double *Limbfit_Parameter;
	double *Sonnefit_Parameter;
	Limbfit_Parameter = new double[Polynomgrad + 1];
	Sonnefit_Parameter = new double[Polynomgrad + 1];
	// Limb-Spektrum mit einem Polynom anfitten
	// später evlt mehrmals mit shifts in beide Richtungen
	// evtl runs-test machen
	Fit_Polynom(Vollfenster_WL, Vollfenster_Limb, Polyfit_Startindex,
				Polyfit_Endindex, Spezfenst.m_Wellenlaengen[Index], Polynomgrad,
				Limbfit_Parameter);
	// Fit sieht ok aus
	// Sonnen-Spektrum mit gleichem Polynom anfitten
	// später evlt mehrmals mit shifts in beide Richtungen
	Fit_Polynom(Vollfenster_WL, Vollfenster_Sonne, Polyfit_Startindex,
				Polyfit_Endindex, Spezfenst.m_Wellenlaengen[Index], Polynomgrad,
				Sonnefit_Parameter);
	// Beide Polynome im gesamten Fenster diskretisieren...
	// dabei 5 mal so viele Punkte nutzen, wie die ursprünglichen Messwerte
	double *Vollfenster_fein_WL;
	double *Vollfenster_fein_Limb;
	double *Vollfenster_fein_Sonne;
	int N_Vollfenster_fein = 1 + (N_Vollfenster - 1) * 5;
	Vollfenster_fein_WL = new double[N_Vollfenster_fein];
	Vollfenster_fein_Limb = new double[N_Vollfenster_fein];
	Vollfenster_fein_Sonne = new double[N_Vollfenster_fein];
	for (int i = 0; i < (N_Vollfenster - 1); i++) {
		Vollfenster_fein_WL[5 * i] = Vollfenster_WL[i];
		Vollfenster_fein_WL[5 * i + 1] = Vollfenster_WL[i] + 0.2 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
		Vollfenster_fein_WL[5 * i + 2] = Vollfenster_WL[i] + 0.4 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
		Vollfenster_fein_WL[5 * i + 3] = Vollfenster_WL[i] + 0.6 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
		Vollfenster_fein_WL[5 * i + 4] = Vollfenster_WL[i] + 0.8 * (Vollfenster_WL[i + 1] - Vollfenster_WL[i]);
	}
	//letzter Punkt, nicht verfünffachen
	Vollfenster_fein_WL[N_Vollfenster_fein - 1] = Vollfenster_WL[N_Vollfenster - 1];
	//Limb und Sonnenspektrum für diskrete Punkte ausrechnen
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		Vollfenster_fein_Limb[i] = 0;
		Vollfenster_fein_Sonne[i] = 0;
	}
	for (int i = 5 * (Polyfit_Startindex - 3); i <= 5 * (Polyfit_Endindex + 3); i++) {
		double h = Vollfenster_fein_WL[i] - Spezfenst.m_Wellenlaengen[Index]; //x-x0
		double Faktor = 1.0;
		Vollfenster_fein_Limb[i] = 0;
		Vollfenster_fein_Sonne[i] = 0;
		for (int j = 0; j <= Polynomgrad; j++) {
			Vollfenster_fein_Limb[i] += Limbfit_Parameter[j] * Faktor;
			Vollfenster_fein_Sonne[i] += Sonnefit_Parameter[j] * Faktor;
			Faktor *= h;
		} //ende for j
	}// ende for i
	//Minimimum beider Polynome bestimmen und Limbspektrum so shiften, dass das
	//Minimum auch auf dem Minimum des Sonnenspektrums liegt
	double Limb_WL_min, Limb_y_min;
	int Limb_Indexmin;
	double Sonne_WL_min, Sonne_y_min;
	int Sonne_Indexmin;
	x_zu_Minimum_von_y_in_Intervall(Vollfenster_fein_WL, Vollfenster_fein_Limb,
									Polyfit_Startindex * 5, Polyfit_Endindex * 5,
									Limb_WL_min, Limb_y_min, Limb_Indexmin);
	x_zu_Minimum_von_y_in_Intervall(Vollfenster_fein_WL, Vollfenster_fein_Sonne,
									Polyfit_Startindex * 5, Polyfit_Endindex * 5,
									Sonne_WL_min, Sonne_y_min, Sonne_Indexmin);
	int shift = Sonne_Indexmin - Limb_Indexmin; //z.b. Sonne 43 Limb 50 shift=-7
	// shiften des Limbspektrums
	// Die Randpunkte sind jetzt falsch, aber die Linie sollte nicht am Rand
	// liegen
	double *vor_shift;
//    cout<<"Sonne_WL_min: "<<Sonne_WL_min<<"\n";
//    cout<<"Limb_WL_min: "<<Limb_WL_min<<"\n";
//    cout<<"Limb_Indexmin: "<<Limb_Indexmin<<"\n";
//    cout<<"Sonne_Indexmin: "<<Sonne_Indexmin<<"\n";
//    cout<<"shift: "<<shift<<"\n";

	vor_shift = new double[N_Vollfenster_fein];
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		vor_shift[i] = Vollfenster_fein_Limb[i];
	}
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		if (!(((i + shift) < 0) || (i + shift >= N_Vollfenster_fein))) {
			Vollfenster_fein_Limb[i + shift] = vor_shift[i];
		}
	}
	delete[] vor_shift;
	//Die beiden gefitteten Spektren dividieren Limb/Sonne
	double *Fit_Quotient;
	double *Messwerte_Quotient;
	Messwerte_Quotient = new double[N_Vollfenster];
	Fit_Quotient = new double[N_Vollfenster_fein];
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] = Vollfenster_Limb[i] / Vollfenster_Sonne[i];
	}
	int Ind = Get_Index(Spezfenst.m_Wellenlaengen[Index]);
	double Delta_WL = (m_Wellenlaengen[Ind + 1] - m_Wellenlaengen[Ind]);
	double Umrechnung = 1 / (Delta_WL * Spezfenst.m_Liniendaten[Index].m_Gamma);
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] *= Umrechnung;
	}
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		Fit_Quotient[i] = Umrechnung * Vollfenster_fein_Limb[i]
						  / Vollfenster_fein_Sonne[i];
		//cout<<Vollfenster_fein_WL[i]<<"\t"<<Fit_Quotient[i]<<"\n";

	}

	//Basislinie abziehen
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] -= a0 + a1 * Vollfenster_WL[i];
	}
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		Fit_Quotient[i]             -= a0 + a1 * Vollfenster_fein_WL[i];
	}
	// Falls an der Stelle des Peaks der Wert jetzt negativ ist, alle positiven
	// Werte abschneiden, sonst negative
	if (Fit_Quotient[Sonne_Indexmin] < 0) {
		for (int i = 0; i < N_Vollfenster_fein; i++) {
			if (Fit_Quotient[i] > 0) {
				Fit_Quotient[i] = 0;
			}
		}
	} else {
		for (int i = 0; i < N_Vollfenster_fein; i++) {
			if (Fit_Quotient[i] < 0) {
				Fit_Quotient[i] = 0;
			}
		}
	}
	//Alles Abschneiden, was nicht zum peak beiträgt
	for (int i = 0; i < N_Vollfenster_fein; i++) {
		if ((Vollfenster_fein_WL[i] < 284.8) || (Vollfenster_fein_WL[i] > 285.6)) {
			Fit_Quotient[i] = 0;
		}
	}

	// Noch ebend eine Fitfunktion in den Quotienten legen
	double Flaeche;
	// Hier Wellenlängen in nm verwendet..
	// as hebt sich mit dem Gitterabstand raus
	Fit_Peak_hyperbolic(Vollfenster_fein_WL, Fit_Quotient, Sonne_WL_min,
						Spezfenst.m_FWHM, Flaeche, 0,
						N_Vollfenster_fein - 1);

	//Plotroutine aufrufen
	if (mache_Fit_Plots == "ja") {
		// Dateinamenschnickschnack
		string s_OrbNum;
		stringstream buf;
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...
		//falls / im Namen ist das schlecht
		string Datnam = sb_basename(m_Dateiname_L1C);
		//TODO Pfad anpassen
		string plot_dir = Arbeitsverzeichnis + "/Plots";
		mkdir(plot_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		buf << Datnam.c_str() << "_" << Spezfenst.m_Spezies_Name.c_str()
			<< "_" << Index << "_" << m_Hoehe_TP << "km.ps";
		string new_datnam(buf.str());
		string s1(plot_dir + "/" + new_datnam);
		//s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern, um
		//später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal, wie die Datei heißt die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Limbdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		buf.str(string());
		buf << "Orbit " << s_OrbNum.c_str() << ", Limb TP:"
			<< " Lat: " << m_Latitude_TP << " deg,"
			<< " Lon: " << m_Longitude_TP << " deg,"
			<< " Hoehe: " << m_Hoehe_TP << " km.";
		string s2(buf.str());
		//cout<<s1<<"\n";
		Plot_Slantcoloumns_polyfit_MgI(Arbeitsverzeichnis.c_str(), s1.c_str(),
									   s2.c_str(),
									   Vollfenster_WL,
									   0.35 * (N_Vollfenster - 1),
									   0.8 * (N_Vollfenster - 1),
									   Vollfenster_fein_WL,
									   0.35 * (N_Vollfenster_fein - 1),
									   0.8 * (N_Vollfenster_fein - 1),
									   Vollfenster_Limb, Vollfenster_fein_Limb,
									   Vollfenster_Sonne, Vollfenster_fein_Sonne,
									   Messwerte_Quotient, Fit_Quotient);
	}
	//Speicher freimachen
	delete[] Basisfenster_WL;
	delete[] Basisfenster_Intensitaet;
	delete[] Vollfenster_WL;
	delete[] Vollfenster_Limb;
	delete[] Vollfenster_Sonne;
	delete[] Vollfenster_fein_WL;
	delete[] Vollfenster_fein_Limb;
	delete[] Vollfenster_fein_Sonne;
	delete[] Limbfit_Parameter;
	delete[] Sonnefit_Parameter;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Saeulendichte_Bestimmen_MgI285nm
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plots_der_Spektren_erzeugen
////////////////////////////////////////////////////////////////////////////////
void Messung_Limb::Plots_der_Spektren_erzeugen(Speziesfenster &Spezfenst,
		int Index, string Arbeitsverzeichnis, string mache_Fit_Plots,
		double *mean_10_20)
{
	// Plot der Spektren Sonne und Limb und Quotient mit und ohne Rauschen
	double *Basisfenster_WL;
	double *Basisfenster_Intensitaet;
	double *Vollfenster_WL;
	double *Vollfenster_Limb;
	double *Vollfenster_Sonne;

	double *Vollfenster_Limb_abs_error;

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
//    int Index_Uebergangs_Wellenlaenge=Get_Index(Spezfenst.m_Wellenlaengen[Index])-Index_Basisfenster_links_min;
	// Speicherplatzbedarf für die Fenster ermitteln
	int Bas_l = (Index_Basisfenster_links_max - Index_Basisfenster_links_min + 1);
	int Bas_r = (Index_Basisfenster_rechts_max - Index_Basisfenster_rechts_min + 1);
	int N_Basis = Bas_l + Bas_r;
	//int N_Peak=Index_Peakfenster_max-Index_Peakfenster_min+1;
	int N_Vollfenster = Index_Basisfenster_rechts_max - Index_Basisfenster_links_min + 1;
	// Speicher anfordern
	Basisfenster_WL = new double[N_Basis];
	Basisfenster_Intensitaet = new double[N_Basis];
	Vollfenster_WL = new double[N_Vollfenster];
	Vollfenster_Limb = new double[N_Vollfenster];
	Vollfenster_Sonne = new double[N_Vollfenster];

	Vollfenster_Limb_abs_error = new double[N_Vollfenster];

	// Basisfenster WL und I auffüllen
	for (int i = 0; i < Bas_l; i++) {
		Basisfenster_WL[i] = this->m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Basisfenster_Intensitaet[i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] =
			this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma_mal_Gitterabstand[Index_Basisfenster_rechts_min + i];
	}

	//Vollfenster Limb und Sonne auffüllen
	for (int i = 0; i < N_Vollfenster; i++) {
		Vollfenster_WL[i] = m_Wellenlaengen[Index_Basisfenster_links_min + i];
		Vollfenster_Limb[i] = m_Intensitaeten[Index_Basisfenster_links_min + i];
		Vollfenster_Sonne[i] = m_Sonne[Index_Basisfenster_links_min + i];
		Vollfenster_Limb_abs_error[i] = Vollfenster_Limb[i]
			* m_Intensitaeten_relativer_Fehler[Index_Basisfenster_links_min + i];
	}

	// TODO braucht man das hier überhaupt noch für irgendwas...wenn nicht weg
	// damit..das verwirrt nur
	//double* Vollfenster_Limb_mittlere_atmo; //zwischen 40 und 60km
	//Vollfenster_Limb_mittlere_atmo=new double[N_Vollfenster];
	//for(int i=0;i<N_Vollfenster;i++)
	//{
	//  Vollfenster_Limb_mittlere_atmo[i]
	//    = mean_10_20[Index_Basisfenster_links_min+i];
	//}

	// Spektrum glätten bei 0 wird nichts geglättet
	// ...könnte man auch auskommentieren
	int smooth_Nachbarn = 0; //1;//2;
	int smooth_Iterationen = 0; //6;//4;
	smooth_data(N_Vollfenster, Vollfenster_Limb, smooth_Nachbarn, smooth_Iterationen);
	smooth_data(N_Vollfenster, Vollfenster_Sonne, smooth_Nachbarn, smooth_Iterationen);

	// linearen Fit des Basisfensters durchführen
	// TODO Basisfenster neu berechnen
	// Proto:
	// Fit_Linear(double* x,double* y, double& a0, double& a1,int Anfangsindex,
	// int Endindex)
	double a0, a1, rms_err_base;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, rms_err_base, 0,
			N_Basis - 1);
	//Die beiden gefitteten Spektren dividieren Limb/Sonne
	double *Messwerte_Quotient;
	double *Messwerte_Quotient_error;  // Es wird angenommen, der Fehler liegt nur in Limb vor
	Messwerte_Quotient = new double[N_Vollfenster];
	Messwerte_Quotient_error = new double[N_Vollfenster];
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] = Vollfenster_Limb[i] / Vollfenster_Sonne[i];
		// Die relativen Fehler von MW_Q und Vf_L sind gleich...Limb ist linear,
		// also auch Fehler linear skalieren
		Messwerte_Quotient_error[i] = Vollfenster_Limb_abs_error[i]
									  / Vollfenster_Sonne[i];

	}
	int Ind = Get_Index(Spezfenst.m_Wellenlaengen[Index]);
	double Delta_WL = (m_Wellenlaengen[Ind + 1] - m_Wellenlaengen[Ind]);
	double Umrechnung = 1 / (Delta_WL * Spezfenst.m_Liniendaten[Index].m_Gamma);
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] *= Umrechnung;
		Messwerte_Quotient_error[i] *= Umrechnung;
	}

	//Basislinie abziehen
	for (int i = 0; i < N_Vollfenster; i++) {
		Messwerte_Quotient[i] -= a0 + a1 * Vollfenster_WL[i];
		// Im besten Fall verändert sich der absolute Fehler nicht...aber der
		// relative, wenn die Baseline im Vergleich zum peak hoch liegt
	}
	// Fehler aus Residuuen abschätzen
	/////////////////////(und nicht den gegebenen Fehler nehmen
	// Zunächst nochmal den Mittelwert bilden
	double *Messwerte_Quotient_stabw;
	Messwerte_Quotient_stabw = new double[N_Vollfenster];
	double Mean = 0;
	for (int i = 0; i < N_Vollfenster; i++)  {
		Mean += Messwerte_Quotient[i];
	}
	Mean /= N_Vollfenster;
	double standardabweichung = 0;
	for (int i = 0; i < N_Vollfenster; i++) {
		standardabweichung += (Messwerte_Quotient[i] - Mean)
							  * (Messwerte_Quotient[i] - Mean);
	}
	standardabweichung /= N_Vollfenster - 1;
	standardabweichung = sqrt(standardabweichung);
	for (int i = 0; i < N_Vollfenster; i++)    {
		Messwerte_Quotient_stabw[i] = standardabweichung;
	}
	// Ende Fehler aus Residuum abschätzen //////////////

	//Plotroutine aufrufen
	if (mache_Fit_Plots == "ja") {
		// Dateinamenschnickschnack
		string s_OrbNum;
		stringstream buf;
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt
		// ...falls / im Namen ist das schlecht
		string Datnam = sb_basename(m_Dateiname_L1C);
		//TODO Pfad anpassen
		string plot_dir = Arbeitsverzeichnis + "/Plots";
		mkdir(plot_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		buf << Datnam.c_str() << "_" << Spezfenst.m_Spezies_Name.c_str()
			<< "_" << Index << "_" << m_Hoehe_TP << "km.ps";
		string new_datnam(buf.str());
		string s1(plot_dir + "/" + new_datnam);
		//s1 ist der Volle Pfad der Datei...diesen kann man wegspeichern,
		//um später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal, wie die Datei heißt die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Limbdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		buf.str(string());
		buf << "Orbit " << s_OrbNum.c_str() << ", Limb TP:"
			<< " Lat: " << m_Latitude_TP << " deg,"
			<< " Lon: " << m_Longitude_TP << " deg,"
			<< " Hoehe: " << m_Hoehe_TP << " km.";
		string s2(buf.str());
		//cout<<s1<<"\n";
		Plot_Spektren_und_Quotient(Arbeitsverzeichnis.c_str(),
								   s1.c_str(), s2.c_str(),
								   Vollfenster_WL, 0,
								   N_Vollfenster - 1,
								   Vollfenster_Limb, Vollfenster_Limb_abs_error,
								   Vollfenster_Sonne, Messwerte_Quotient,
								   Messwerte_Quotient_error);
		/*Plot_Quotient_mit_Fehler(Arbeitsverzeichnis.c_str(),s1.c_str(), s2.c_str(),
		                   Vollfenster_WL,0 ,N_Vollfenster-1,
		                   Messwerte_Quotient,Messwerte_Quotient_error,
		                   Messwerte_Quotient_stabw);*/
	}
	///////////////////////////////////////////////
	// Fit der Säulendichte  /////////////
	double Flaeche;
	// Hier Wellenlängen in nm verwendet..
	// das hebt sich mit dem Gitterabstand raus
	Fit_Peak_hyperbolic(Vollfenster_WL, Messwerte_Quotient,
						Spezfenst.m_Wellenlaengen[Index],
						Spezfenst.m_FWHM, Flaeche, 0,
						N_Vollfenster - 1);
	//Fehler des Fits bestimmen... da Peakfenster_Intensitaet nicht mehr
	//gebraucht wird die Basislinie für die Fehlerberechnung wieder aufaddiert
	m_Zeilendichte = Flaeche;
	// Funktion double Messung_Limb::Evaluate_Error_primitive(double* x,
	// double* y, double a0,double a1, double A, double FWHM, double x0,
	// int Anfangsindex, int Endindex)
	m_Fehler_Zeilendichten =
		Evaluate_Error_primitive(Vollfenster_WL,
								 Messwerte_Quotient, a0, a1, m_Zeilendichte,
								 Spezfenst.m_FWHM,
								 Spezfenst.m_Wellenlaengen[Index], 0,
								 N_Vollfenster - 1);
	//////////////////////////////////////////////
	/////////////////////////////////////////////


	//Speicher freimachen
	delete[] Basisfenster_WL;
	delete[] Basisfenster_Intensitaet;
	delete[] Vollfenster_WL;
	delete[] Vollfenster_Limb;
	delete[] Vollfenster_Sonne;
	delete[] Vollfenster_Limb_abs_error;
	delete[] Messwerte_Quotient_error;
	delete[] Messwerte_Quotient_stabw;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Plots_der_Spektren_erzeugen
////////////////////////////////////////////////////////////////////////////////
void Messung_Limb::moving_average(int window_size)
{
	my_moving_average(m_Intensitaeten, window_size);
}
void Messung_Limb::savitzky_golay(int window_size)
{
	my_savitzky_golay(m_Intensitaeten, window_size);
}

double Messung_Limb::msise_temperature(Konfiguration &Konf)
{
	struct nrlmsise_output output;
	struct nrlmsise_input input;
	struct nrlmsise_flags flags;

	// set the flags
	flags.switches[0] = 0;
	for (int i = 1; i < 24; i++)
		flags.switches[i] = 1;

	// construct the input for the temperature calculation
	input.year = m_Jahr; // year, but ignored

	// get the day of the year
	int days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	if (m_Jahr % 4 == 0 && !(m_Jahr % 100 == 0 && m_Jahr % 400 != 0))
		days[1] = 29;

	input.doy = 0;
	for (int i = 0; i < (m_Monat - 1); i++) {
		input.doy += days[i];
	}
	input.doy += m_Tag;

	// ut seconds in day
	input.sec = m_Stunde * 3600. + m_Minute * 60. + m_Sekunde;

	// geo data
	input.alt = m_Hoehe_TP;
	input.g_lat = m_Latitude_TP;
	input.g_long = m_Longitude_TP;
	// local apparent solar time (quick default)
	input.lst = input.sec / 3600. + input.g_long / 15.;

	// solar data from spidr data files
	double f107 = spidr_value_from_file(m_Jahr, m_Monat, m_Tag,
			Konf.m_Pfad_f107_index);
	double ap = spidr_value_from_file(m_Jahr, m_Monat, m_Tag,
			Konf.m_Pfad_Ap_index);
	std::cout << "# msis parameters: f10.7 = " << f107
		<< ", ap = " << ap << std::endl;
	input.f107A = f107;
	input.f107 = f107;
	input.ap = ap;

	gtd7(&input, &flags, &output);

	total_number_density = output.d[0] + output.d[1] + output.d[2]
		+ output.d[3] + output.d[4] + output.d[6] + output.d[7];

	/*
	// http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html
	// refraction at pressure p [Pa] and T [Celsius] is related to this:
	// n' − 1 = (n − 1)*p[1 + p*(60.1 - 0.972*T)*10^-10]/[96095.43(1 + 0.003661*T)]
	double T_K = output.t[1];
	double T_C = T_K - 273.15;
	double p_Pa = 1.e6 * total_number_density * 1.3806504e-23 * T_K;
	double n_f = p_Pa*(1.+p_Pa*(60.1-0.972*T_C)*1.e-10)/(96095.43*(1.+0.003661*T_C));
	// */

	return output.t[1];
}

Ausgewertete_Messung_Limb Messung_Limb::Ergebnis_Zusammenfassen()
{
	Ausgewertete_Messung_Limb aus;
	//Ergebnisse
	aus.m_Zeilendichte = this->m_Zeilendichte;
	aus.m_Fehler_Zeilendichten = this->m_Fehler_Zeilendichten;
	// total number density
	aus.total_number_density = this->total_number_density;
	//Zwischenergebnisse
	aus.m_Deklination = this->m_Deklinationswinkel;
	aus.m_Sonnen_Longitude = this->m_Sonnen_Longitude;
	aus.m_Wellenlaenge = 0;
	// Nullinitialisierung...
	// die Wellenlänge des Übergangs steckt im Speziesfenster
	//Datum
	aus.m_Jahr = this->m_Jahr;
	aus.m_Monat = this->m_Monat;
	aus.m_Tag = this->m_Tag;
	aus.m_Stunde = this->m_Stunde;
	aus.m_Minute = this->m_Minute;
	aus.m_Sekunde = this->m_Sekunde;
	// Geolocation
	aus.m_Latitude_Sat = this->m_Latitude_Sat;
	aus.m_Longitude_Sat = this->m_Longitude_Sat;
	aus.m_Hoehe_Sat = this->m_Hoehe_Sat;
	aus.m_Latitude_TP = this->m_Latitude_TP;
	aus.m_Longitude_TP = this->m_Longitude_TP;
	aus.m_Hoehe_TP = this->m_Hoehe_TP;
	aus.m_Erdradius = this->m_Erdradius;
	// phase of orbit (0...1)
	aus.m_orbit_phase = this->m_orbit_phase;
	aus.center_lat = this->center_lat;
	aus.center_lon = this->center_lon;
	return aus;
}//Ausgewertete_Messung_Limb Ergebnis_Zusammenfassen() ende
//========================================

//========================================
//Methoden ende
//========================================

////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen
//
////////////////////////////////////////////////////////////////////////////////
void Messung_Limb::Ausgabe_in_Datei(string Dateiname)
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
	fprintf(outfile, "%20s %1.5E\n", "m_Latitude_TP: ", m_Latitude_TP);
	fprintf(outfile, "%20s %1.5E\n", "m_Longitude_TP: ", m_Longitude_TP);
	fprintf(outfile, "%20s %1.5E\n", "m_Hoehe_TP: ", m_Hoehe_TP);
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

