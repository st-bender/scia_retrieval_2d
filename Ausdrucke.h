/*
 * Ausdrucke.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 07.05.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
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

void Plot_2xy(std::string Arbeitsverzeichnis, std::string Dateiname,
			 std::string title, std::string xlabel, std::string ylabel,
			 std::vector<double> &x1, std::vector<double> &y1,
			 std::vector<double> &x2, std::vector<double> &y2,
			 int Startindex, int Endindex,
			 double Mittelwert, double Fehler, bool keep = false, bool run = true);
void Plot_Slantcoloumns_polyfit_MgI(std::string Arbeitsverzeichnis, std::string Dateiname,
								   std::string title, double *WL, int Startind,
								   int Endind, double *WL_fein,
								   int Startind_fein, int Endind_fein,
								   double *Limb, double *Fitlimb, double *Sonne,
								   double *Fitsonne, double *Quotient,
								   double *Fitquotient);
void Plot_Spektren_und_Quotient(std::string Arbeitsverzeichnis, std::string Dateiname,
							   std::string title, double *WL, int Startind,
							   int Endind, double *Limb, double *Limb_error,
							   double *Sonne, double *Quotient,
							   double *Quotient_error);
void Plot_Quotient_mit_Fehler(std::string Arbeitsverzeichnis, std::string Dateiname,
							 std::string title, double *WL, int Startind, int Endind,
							 double *Quotient, double *Quotient_error,
							 double *Quotient_stabw);


//int Plot_Slantcoloumns_smooth
void Plots_Zusammenfassen(std::string Pfad_multips2pdf, std::string Pfad_multips2ps,
						 std::string Name_pdf_Datei,
						 std::vector<std::string> Liste_der_ps_Dateinamen);


#endif /* AUSDRUCKE_HH_ */
