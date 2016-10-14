/*
 * Glaetten.cpp
 *
 *  Created on: 05.11.2010
 *      Author: martin
 */

#include <iostream>
#include <vector>
#include <deque>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "MPL_Matrix.h"

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B,
			int *LDB, int *INFO);
}

using std::vector;

void smooth_data(int Datasize, double *Data, int Anzahl_Nachbarn_eine_Seite,
		int Zahl_der_Iterationen)
{
	double *Data_old;
	Data_old = new double[Datasize];
	// Altes Feld übergeben
	for (int i = 0; i < Datasize; i++) {
		Data_old[i] = Data[i];
	}
	// randpunkte werden nicht verändert
	for (int j = 0; j < Zahl_der_Iterationen; j++) {
		for (int i = Anzahl_Nachbarn_eine_Seite;
				i < (Datasize - Anzahl_Nachbarn_eine_Seite - 1); i++) {
			//Für jeden Punkt innerhalb des Glättungsintervalls über
			//nachbarpunkte mitteln
			for (int k = i - Anzahl_Nachbarn_eine_Seite;
					k <= i + Anzahl_Nachbarn_eine_Seite; k++) {
				Data[i] += Data_old[k];
			}
			Data[i] /= 2 * Anzahl_Nachbarn_eine_Seite + 1;
		}
		//Altes Datenfeld anpassen für nächsten iterationsschrittt
		for (int i = 0; i < Datasize; i++) {
			Data_old[i] = Data[i];
		}
	}
	// Speicher freigeben
	delete[] Data_old;
}

double my_clamp(double d, double min, double max)
{
	const double t = d < min ? min : d;
	return t > max ? max : t;
}
/*
 * smooth base function (for the partition of unity)
 * theta(t) = 0, t <= 0
 * theta(t) = exp(-1/t), t > 0
 */
double my_theta(double t)
{
	if (t <= 0.) return 0.;
	return std::exp(-1. / t);
}
/*
 * smooth transition function
 * phi(x) = 0, x < 0
 * phi(x) = 1, x >= 1
 */
double my_phi(double x)
{
	double tx = my_theta(x);
	double t1mx = my_theta(1. - x);
	return tx / (tx + t1mx);
}
/*
 * hard step function at 1.0
 * phi(x) = 0, x < 1
 * phi(x) = 1, x >= 1
 */
double hardstep(double x)
{
	return x >= 1.0;
}
/*
 * C1 step function like my_phi(x) from wikipedia
 * https://en.wikipedia.org/wiki/Smoothstep
 */
double smoothstep(double x)
{
	x = my_clamp(x, 0., 1.);
	// Evaluate polynomial
	return x*x*(3 - 2*x);
}
/*
 * C2 step function from wikipedia
 * https://en.wikipedia.org/wiki/Smoothstep
 */
double smootherstep(double x)
{
	x = my_clamp(x, 0., 1.);
	// Evaluate polynomial
	return x*x*x*(x*(x*6 - 15) + 10);
}
/*
 * flexible transition function
 * using the first argument as transition from zero to one
 *
 * Phi(x) = 0, x < a - w
 * Phi(x) = 1, a <= x <= b
 * Phi(x) = 0, x > b + w
 */
double Phi_func(double (*trans01)(double),
				double a, double b, double w, double x)
{
	return trans01((x - a + w) / w) * trans01((b - x + w) / w);
}

void my_moving_average(vector<double> &y, int ws)
{
	vector<double> y_neu;
	vector<double>::iterator y_it;
	const int wsh = ws / 2;

	double sum = accumulate(y.begin(), y.begin() + wsh, 0.);
	int pts = wsh;

	for (y_it = y.begin(); y_it != y.end(); ++y_it) {
		double a = 0., b = 0.;

		if (y_it <= y.begin() + wsh) {
			a = *(y_it + wsh);
			++pts;
		} else if (y_it >= y.end() - wsh) {
			b = *(y_it - wsh - 1);
			--pts;
		} else {
			a = *(y_it + wsh);
			b = *(y_it - wsh - 1);
		}

		sum += a - b;
		y_neu.push_back(sum / pts);
	}

	y = y_neu;
}

void my_convolution_1d(vector<double> &y, vector<double> &weights)
{
	size_t i;
	const size_t ws = weights.size(), wsh = ws / 2;

	vector<double> y_neu;
	vector<double>::iterator y_it;

	for (y_it = y.begin(); y_it != y.end(); ++y_it) {
		double avg = 0.;
		double wnorm = 0.;
		size_t start_shift = 0, end_shift = 0;

		if (y_it < y.begin() + wsh)
			start_shift = distance(y_it, y.begin() + wsh);

		if (y_it >= y.end() - wsh)
			end_shift = distance(y.end() - wsh, y_it);

		for (i = start_shift; i < ws - end_shift; i++) {
			avg += (*(y_it - wsh + i)) * weights.at(i);
			wnorm += weights.at(i);
		}

		avg /= wnorm;
		y_neu.push_back(avg);
	}

	y = y_neu;
}

int my_savitzky_golay(vector<double> &y, int ws)
{
	const double weights5[5] = { -3., 12., 17., 12., -3. };
	const double weights7[7] = { -2., 3., 6., 7., 6., 3., -2. };
	const double weights9[9] = { -21., 14., 39., 54., 59., 54., 39., 14., -21 };

	vector<double> wgts5(weights5, weights5 + 5);
	vector<double> wgts7(weights7, weights7 + 7);
	vector<double> wgts9(weights9, weights9 + 9);

	switch (ws) {
	case 5:
		my_convolution_1d(y, wgts5);
		break;
	case 7:
		my_convolution_1d(y, wgts7);
		break;
	case 9:
		my_convolution_1d(y, wgts9);
		break;
	default:
		std::cerr << "unsupported window size for Savitzky-Golay." << std::endl;
		return -1;
	}

	return 0;
}

void my_gauss_blur_1d(vector<double> &y)
{
	const double wgts[7] = { 0.0044, 0.054, 0.242, 0.399, 0.242, 0.054, 0.0044 };

	vector<double> weights(wgts, wgts + 7);

	my_convolution_1d(y, weights);
}

/* transforms the sao solar reference (0.01 nm resolution)
 * to the sciamachy resolution (0.11 nm), FWHM = 0.22 nm */
void my_sciamachy_blur(vector<double> &y)
{
	const double wgts[99] = { .01036709959893280509, .01125594763314678180,
		.01224195115079379920, .01333810263675415044, .01455946003010629420,
		.01592352876117980314, .01745072461985826728, .01916493674234588559,
		.02109421512879353307, .02327161372799369268, .02573622872401204464,
		.02853448287139796698, .03172172140873178286, .03536420439181033864,
		.03954160579083501440, .04435016350408290365, .04990666941003913911,
		.05635354856088212371, .06386535677234359197, .07265713303153179253,
		.08299518623931035337, .09521108608535462417, .10971987877325275129,
		.12704387424595680166, .14784376394825968360, .17295932601533818760,
		.20346252348277792508, .24072628774254172983, .28651242260620790224,
		.34308124003100880633, .41332246271612842636, .50089912671489344275,
		.61037927956088312110, .74729724035889854875, .91802854227707242631,
		1.12927701981369864455, 1.38688141382934607104, 1.69364419934577293193,
		2.04617344581160470660, 2.43155922406783040600, 2.82596221301082016019,
		3.19773981108263911061, 3.51578751556480512812, 3.75955643096287940250,
		3.92480353990930230846, 4.02202126872896617699, 4.06983092244636659800,
		4.08787956882413925119, 4.09206739791390581960, 4.09234689162320942107,
		4.09206739791390581960, 4.08787956882413925119, 4.06983092244636659800,
		4.02202126872896617699, 3.92480353990930230846, 3.75955643096287940250,
		3.51578751556480512812, 3.19773981108263911061, 2.82596221301082016019,
		2.43155922406783040600, 2.04617344581160470660, 1.69364419934577293193,
		1.38688141382934607104, 1.12927701981369864455, .91802854227707242631,
		.74729724035889854875, .61037927956088312110, .50089912671489344275,
		.41332246271612842636, .34308124003100880633, .28651242260620790224,
		.24072628774254172983, .20346252348277792508, .17295932601533818760,
		.14784376394825968360, .12704387424595680166, .10971987877325275129,
		.09521108608535462417, .08299518623931035337, .07265713303153179253,
		.06386535677234359197, .05635354856088212371, .04990666941003913911,
		.04435016350408290365, .03954160579083501440, .03536420439181033864,
		.03172172140873178286, .02853448287139796698, .02573622872401204464,
		.02327161372799369268, .02109421512879353307, .01916493674234588559,
		.01745072461985826728, .01592352876117980314, .01455946003010629420,
		.01333810263675415044, .01224195115079379920, .01125594763314678180,
		.01036709959893280509 };

	vector<double> weights(wgts, wgts + 99);

	my_convolution_1d(y, weights);
}

/* lowess smoothing, code inspired by the biopython module found in
 * <biopython>/Bio/Statistics/lowess.py
 *
 * For more information, see
 *
 * William S. Cleveland: "Robust locally weighted regression and smoothing
 * scatterplots", Journal of the American Statistical Association, Dec 1979,
 * volume 74, number 368, pp. 829-836.

 * William S. Cleveland and Susan J. Devlin: "Locally weighted regression: An
 * approach to regression analysis by local fitting", Journal of the American
 * Statistical Association, Sep 1988, volume 83, number 403, pp. 596-610.
 */
void my_lowess(vector<double> &x, vector<double> &y, double f)
{
	size_t n = x.size();
	size_t r = ceil(f * n);
	vector<double> y_neu;

	for (size_t i = 0; i < n; i++) {
		vector<double> dist, dist_sort, wgts;
		// calculate the distances
		for (size_t j = 0; j < n; j++) {
			double d = abs(x.at(i) - x.at(j));
			dist.push_back(d);
			dist_sort.push_back(d);
		}
		// sort the distances
		sort(dist_sort.begin(), dist_sort.end());

		// calculate the weights
		for (size_t j = 0; j < n; j++) {
			double w = dist.at(j) / dist_sort.at(r);
			if (w >= 0. && w < 1.) {
				w = 1. - w * w * w;
				w = w * w * w;
			} else w = 0.;
			wgts.push_back(w);
		}

		// sums for the weighted linear regression
		double w_x_sum = 0., w_y_sum = 0.;
		double w_xx_sum = 0., w_xy_sum = 0.;
		double w_sum = 0.;

		for (size_t j = 0; j < n; j++) {
			double w = wgts.at(j), a = x.at(j), b = y.at(j);
			w_sum += w;
			w_x_sum += w * a;
			w_y_sum += w * b;
			w_xx_sum += w * a * a;
			w_xy_sum += w * a * b;
		}

		double det = w_sum * w_xx_sum - w_x_sum * w_x_sum;
		double beta1 = (w_xx_sum * w_y_sum - w_x_sum * w_xy_sum) / det;
		double beta2 = (w_sum * w_xy_sum - w_x_sum * w_y_sum) / det;
		double yval = beta1 + beta2 * x.at(i);

		y_neu.push_back(yval);
	}
	y = y_neu;
}

/* linear equation solver helper function
 * LHS = Ax = b = RHS
 * the original RHS is replaced by the solution x */
int my_solve(MPL_Matrix &LHS, MPL_Matrix &RHS)
{
	// Fortran Matrizen sind zu C++ Matrizen transponiert
	MPL_Matrix A = LHS.transponiert_full();
	// N ist Anzahl der Gitterpunkte
	int N = LHS.m_Zeilenzahl;
	int M = RHS.m_Spaltenzahl;
	// array mit der Pivotisierungsmatrix sollte so groß wie N sein,
	int *IPIV = new int[N];
	// Spalten von RHS 1 nehmen, um keine C/Fortran Verwirrungen zu provozieren
	int NRHS = M;
	int LDA = N;
	int LDB = N;
	int INFO;

	// AUFRUF A ist LHS.transponiert und B ist RHS
	dgesv_(&N, &NRHS, A.m_Elemente, &LDA, IPIV, RHS.m_Elemente, &LDB, &INFO);

	delete[] IPIV;

	return INFO;
}

/* Whittaker smoothing for background subtraction
 * method described in Anal. Chem. 75, 3631--3636 (2003)
 * returns the smoothed vector with the same length as the input vector (y),
 * and takes a weight vector (w) (0 for points to be skipped, 1 otherwise). */
std::vector<double> my_whittaker_smooth(std::vector<double> &y,
		std::vector<double> &w, int order, double lambda, double &err)
{
	int m = y.size();
	MPL_Matrix dummy(m, m);
	MPL_Matrix E = dummy.unity();
	MPL_Matrix D = E.row_diff();
	MPL_Matrix Y(m, 1), Z(m, 1), W(m, m);

	while (--order > 0)
		D = D.row_diff();

	// prepare RHS and W
	for (int i = 0; i < m; i++) {
		Z(i) = Y(i) = w.at(i) * y.at(i);
		W(i, i) = w.at(i);
	}

	// prepare LHS
	MPL_Matrix A = W + lambda * D.transponiert() * D;
	my_solve(A, Z);

	/* According to the paper, the Whittaker smoothing matrix H
	 * is given by H = (W + lambda*D^T*D)^-1 * W = A^-1 * W. */
	MPL_Matrix H(E);
	my_solve(A, H);
	H = H * W;

	// calculate the rms error of the smoothed points
	double N = std::accumulate(w.begin(), w.end(), 0.);
	MPL_Matrix R = Y - Z;
	MPL_Matrix Res = R.transponiert() * W * R;
	err = std::sqrt(Res(0, 0) / (N - 1.25 * H.trace() + 0.5));

	// generate return vector
	std::vector<double> z(Z.m_Elemente, Z.m_Elemente + Z.m_Elementanzahl);

	return z;
}

/* not really smoothing functions but this file seems to be the best place for now */
double interpolate(std::vector<double> &x, std::vector<double> &y, double x0)
{
	long i;
	std::vector<double>::iterator x_it;
	x_it = std::upper_bound(x.begin(), x.end(), x0);

	// linear extrapolation outside the range
	if (x_it == x.begin())
		i = 0;
	else if (x_it == x.end())
		i = x.size() - 2;
	else
		i = distance(x.begin(), x_it) - 1;

	return y.at(i)
		+ (x0 - x.at(i)) * (y.at(i) - y.at(i + 1)) / (x.at(i) - x.at(i + 1));
}

double fit_spectra(std::vector<double> &x, std::vector<double> &y)
{
	double sum_gy = std::inner_product(x.begin(), x.end(), y.begin(), 0.);
	double sum_gg = std::inner_product(x.begin(), x.end(), x.begin(), 0.);

	return sum_gy / sum_gg;
}

double n_air(double wl)
{
	double sigma = 1.e6 / (wl * wl);
	/* Peck and Reeder, J. Optic. Soc. Am., vol. 62, no. 8, pp. 958--962, 1972
	 * doi: 10.1364/JOSA.62.000958 */
	return 1. + (8060.51 + 2480990. / (132.274 - sigma)
					+ 17455.7 / (39.32957 - sigma)) * 1.e-8;
}
/* the rough rayleigh cross section for wl in [nm] in [cm^2] */
double sigma_rayleigh(double wl)
{
	/* ref.: planet. space sci., vol. 32, no. 6, pp. 785--790, 1984
	 * http://dx.doi.org/10.1016/0032-0633(84)90102-8 */
	std::vector<double> wl_Fk{ 200, 205, 210, 215, 220, 225, 230, 240, 250, 260,
		270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400,
		450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 };
	std::vector<double> F_ks{ 1.080, 1.077, 1.074, 1.072, 1.070, 1.068, 1.066, 1.064,
		1.062, 1.060, 1.059, 1.057, 1.056, 1.055, 1.055, 1.054, 1.053, 1.053, 1.052,
		1.052, 1.052, 1.051, 1.051, 1.051, 1.050, 1.049, 1.049, 1.048, 1.048,
		1.048, 1.048, 1.047, 1.047, 1.047, 1.047, 1.047 };
	double n = n_air(wl);
	double nm1_div_NA = (n - 1.) / 2.687e-2;
	double F_k = interpolate(wl_Fk, F_ks, wl);

	return 32. * M_PI*M_PI*M_PI / 3. * nm1_div_NA*nm1_div_NA
		* F_k / (wl*wl*wl*wl) * 1.e-14;
}
double shift_wavelength(double wl)
{
	return wl / n_air(wl);
}
double spidr_value_from_file(int year, int month, int day,
		std::string filename, double defvalue, unsigned lag)
{
	double ret = defvalue;
	std::deque<double> lagged_vals;
	std::string line, date;
	std::stringstream ss;
	size_t pos;

	// construct the date string from the variables
	ss << year
		<< "-" << std::setw(2) << std::setfill('0') << month
		<< "-" << std::setw(2) << std::setfill('0') << day;
	ss >> date;

	std::ifstream f(filename.c_str());
	if (!f.is_open()) {
		std::cerr << "Error opening `" << filename << "'." << std::endl;
		return ret;
	}
	while (std::getline(f, line)) {
		pos = line.find(date);
		if (pos == std::string::npos && lag > 0) {
			std::istringstream iss(line);
			std::string dummy1, dummy2;
			double val;
			// skip the first two items (date and time)
			iss >> dummy1 >> dummy2;
			// the third is what we need
			iss >> val;
			lagged_vals.push_front(val);
			if (lagged_vals.size() > lag)
				lagged_vals.pop_back();
		}
		if (pos != std::string::npos) {
			if (lag > 0)
				ret = lagged_vals.at(lag - 1);
			else {
				std::istringstream iss(line);
				std::string dummy1, dummy2;
				// skip the first two items (date and time)
				iss >> dummy1 >> dummy2;
				// the third is what we need
				iss >> ret;
			}
			break;
		}
	}

	return ret;
}
