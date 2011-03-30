/*
 * Glaetten.cpp
 *
 *  Created on: 05.11.2010
 *      Author: martin
 */

#include <iostream>
#include <vector>
#include <iterator>

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
	int i, wsh = ws / 2;

	for (y_it = y.begin(); y_it != y.end(); ++y_it) {
		double avg = 0.;

		if (y_it < y.begin() + wsh || y_it >= y.end() - wsh) {
			y_neu.push_back(*y_it);
			continue;
		}

		for (i = 0; i < ws; i++)
			avg += *(y_it - wsh + i);

		avg /= ws;
		y_neu.push_back(avg);
	}

	y = y_neu;

	return 0;
}

int my_savitzky_golay(vector<double> &y, int ws)
{
	int i, wsh = ws / 2, w_ind = (ws - 5) / 2;
	const double weights5[5] = { -3., 12., 17., 12., -3. };
	const double weights7[7] = { -2., 3., 6., 7., 6., 3., -2. };
	const double weights9[9] = { -21., 14., 39., 54., 59., 54., 39., 14., -21 };
	const double weightss[3] = { 35., 21., 231. };

	vector<double> y_neu;
	vector<double>::iterator y_it;
	vector<vector<double> > weights;
	vector<double> wgts5(weights5, weights5 + 5);
	vector<double> wgts7(weights7, weights7 + 7);
	vector<double> wgts9(weights9, weights9 + 9);

	weights.push_back(wgts5);
	weights.push_back(wgts7);
	weights.push_back(wgts9);

	if (ws != 5 && ws != 7 && ws != 9) {
		cerr << "unsupported window size for Savitzky-Golay." << endl;
		return -1;
	}

	for (y_it = y.begin(); y_it != y.end(); ++y_it) {
		double avg = 0.;

		if (y_it < y.begin() + wsh || y_it >= y.end() - wsh) {
			y_neu.push_back(*y_it);
			continue;
		}

		for (i = 0; i < ws; i++) {
			avg += (*(y_it - wsh + i)) * weights[w_ind][i];
		}

		avg /= weightss[w_ind];
		y_neu.push_back(avg);
	}

	y = y_neu;

	return 0;
}
