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

int Plot_2xy(std::string Arbeitsverzeichnis, std::string Dateiname,
			 std::string title, std::string xlabel, std::string ylabel,
			 std::vector<double> &x1, std::vector<double> &y1,
			 std::vector<double> &x2, std::vector<double> &y2,
			 int Startindex, int Endindex,
			 double Mittelwert, double Fehler);
int Plot_Slantcoloumns_polyfit_MgI(std::string Arbeitsverzeichnis, std::string Dateiname,
								   std::string title, double *WL, int Startind,
								   int Endind, double *WL_fein,
								   int Startind_fein, int Endind_fein,
								   double *Limb, double *Fitlimb, double *Sonne,
								   double *Fitsonne, double *Quotient,
								   double *Fitquotient);
int Plot_Spektren_und_Quotient(std::string Arbeitsverzeichnis, std::string Dateiname,
							   std::string title, double *WL, int Startind,
							   int Endind, double *Limb, double *Limb_error,
							   double *Sonne, double *Quotient,
							   double *Quotient_error);
int Plot_Quotient_mit_Fehler(std::string Arbeitsverzeichnis, std::string Dateiname,
							 std::string title, double *WL, int Startind, int Endind,
							 double *Quotient, double *Quotient_error,
							 double *Quotient_stabw);


//int Plot_Slantcoloumns_smooth
int Plots_Zusammenfassen(std::string Pfad_multips2pdf, std::string Pfad_multips2ps,
						 std::string Name_pdf_Datei,
						 std::vector<std::string> Liste_der_ps_Dateinamen);


#endif /* AUSDRUCKE_HH_ */
