/*
 * Glaetten.cpp
 *
 *  Created on: 05.11.2010
 *      Author: martin
 */

#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>

using namespace std;

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

int my_moving_average(vector<double> &y, int ws)
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

	return 0;
}

int my_convolution_1d(vector<double> &y, vector<double> &weights)
{
	int i;
	const int ws = weights.size(), wsh = ws / 2;

	vector<double> y_neu;
	vector<double>::iterator y_it;

	for (y_it = y.begin(); y_it != y.end(); ++y_it) {
		double avg = 0.;
		double wnorm = 0.;
		int start_shift = 0, end_shift = 0;

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

	return 0;
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
		return my_convolution_1d(y, wgts5);
	case 7:
		return my_convolution_1d(y, wgts7);
	case 9:
		return my_convolution_1d(y, wgts9);
	default:
		cerr << "unsupported window size for Savitzky-Golay." << endl;
		return -1;
	}

	return 0;
}

int my_gauss_blur_1d(vector<double> &y)
{
	const double wgts[7] = { 0.0044, 0.054, 0.242, 0.399, 0.242, 0.054, 0.0044 };

	vector<double> weights(wgts, wgts + 7);

	return my_convolution_1d(y, weights);
}

/* not really a smoothing function but this file seems to be the best place for now */
double shift_wavelength(double wl)
{
	double sigma = 1.e6 / (wl * wl);
	double n_air = 1. + 0.000064328 + 0.0294981 / (146. - sigma)
		+ 0.0002554 / (41. - sigma);

	return wl / n_air;
}
