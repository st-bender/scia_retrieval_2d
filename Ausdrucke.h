/*
 * Ausdrucke.h
 *
 *  Created on: 07.05.2010
 *      Author: martin
 */
/******************************************************************************
 * Dieses Paket enthält ein par Funktionen mit den Plots erzeugt werden
 * können. Dafür wird gnuplot benötigt.
 *
 * Zunächst wird eine Gnuplotskriptdatei geschrieben. "Skript.xxx"
 * Diese wird dann ausgeführt (spawn).
 * Wenn der Plot abgeschlossen ist, so wird die Skriptdatei wieder gelöscht.
 *
 ******************************************************************************/

#ifndef AUSDRUCKE_HH_
#define AUSDRUCKE_HH_

#include<string>
#include<vector>

using namespace std;

int Plot_2xy(string Arbeitsverzeichnis, string Dateiname,
			 string title, string xlabel, string ylabel,
			 vector<double> &x1, vector<double> &y1,
			 vector<double> &x2, vector<double> &y2,
			 int Startindex, int Endindex,
			 double Mittelwert, double Fehler);
int Plot_Slantcoloumns_polyfit_MgI(string Arbeitsverzeichnis, string Dateiname,
								   string title, double *WL, int Startind,
								   int Endind, double *WL_fein,
								   int Startind_fein, int Endind_fein,
								   double *Limb, double *Fitlimb, double *Sonne,
								   double *Fitsonne, double *Quotient,
								   double *Fitquotient);
int Plot_Spektren_und_Quotient(string Arbeitsverzeichnis, string Dateiname,
							   string title, double *WL, int Startind,
							   int Endind, double *Limb, double *Limb_error,
							   double *Sonne, double *Quotient,
							   double *Quotient_error);
int Plot_Quotient_mit_Fehler(string Arbeitsverzeichnis, string Dateiname,
							 string title, double *WL, int Startind, int Endind,
							 double *Quotient, double *Quotient_error,
							 double *Quotient_stabw);


//int Plot_Slantcoloumns_smooth
int Plots_Zusammenfassen(string Pfad_multips2pdf, string Pfad_multips2ps,
						 string Name_pdf_Datei,
						 vector<string> Liste_der_ps_Dateinamen);


#endif /* AUSDRUCKE_HH_ */
