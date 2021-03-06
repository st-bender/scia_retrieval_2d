/*
 * Glaetten.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 05.11.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#ifndef GLAETTEN_HH_
#define GLAETTEN_HH_

void smooth_data(int Datasize, double *Data, int Anzahl_Nachbarn_eine_Seite,
		int Zahl_der_Iterationen);

double my_clamp(double d, double min, double max);
double my_phi(double x);
double hardstep(double x);
double smoothstep(double x);
double smootherstep(double x);
double Phi_func(double (*f)(double), double a, double b, double w, double x);
std::vector<double> my_moving_average(std::vector<double> &y, int ws);
std::vector<double> my_savitzky_golay(std::vector<double> &y, int ws);
std::vector<double> my_gauss_blur_1d(std::vector<double> &y);
std::vector<double> my_sciamachy_blur(std::vector<double> &y);
std::vector<double> my_lowess(std::vector<double> &x, std::vector<double> &y, double f);
std::vector<double> my_whittaker_smooth(std::vector<double> &y,
		std::vector<double> &w, int order, double lambda, double &err);

double interpolate(std::vector<double> &x, std::vector<double> &y, double x0);
double fit_spectra(std::vector<double> &x, std::vector<double> &y);
double n_air(double wl);
double sigma_rayleigh(double wl);
double shift_wavelength(double wl);
double spidr_value_from_file(int year, int month, int day,
		std::string filename, double defvalue = 0.0, unsigned lag = 0);

#endif /* GLAETTEN_HH_ */
