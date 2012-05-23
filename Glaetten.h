/*
 * Glaetten.h
 *
 *  Created on: 05.11.2010
 *      Author: martin
 */

#ifndef GLAETTEN_HH_
#define GLAETTEN_HH_

void smooth_data(int Datasize, double *Data, int Anzahl_Nachbarn_eine_Seite,
		int Zahl_der_Iterationen);

int my_moving_average(std::vector<double> &y, int ws);
int my_savitzky_golay(std::vector<double> &y, int ws);
int my_gauss_blur_1d(std::vector<double> &y);
int my_sciamachy_blur(std::vector<double> &y);
int my_lowess(std::vector<double> &x, std::vector<double> &y, double f);
std::vector<double> my_whittaker_smooth(std::vector<double> &y,
		std::vector<double> &w, int order, double lambda, double &err);

double interpolate(std::vector<double> &x, std::vector<double> &y, double x0);
double n_air(double wl);
double shift_wavelength(double wl);
double spidr_value_from_file(int year, int month, int day,
		std::string filename);

#endif /* GLAETTEN_HH_ */
