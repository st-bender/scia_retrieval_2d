/*
 * Ausdrucke.cpp
 *
 *  Created on: 07.05.2010
 *      Author: martin
 */
#include "Ausdrucke.h"
#include <string>
#include<vector>
#include<fstream>
#include<iostream>
#include <cstdlib>
#include <sstream>
#include <algorithm>

using namespace std;
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plot_2xy
////////////////////////////////////////////////////////////////////////////////
int Plot_2xy(string Arbeitsverzeichnis, string Dateiname,
			 string title, string xlabel, string ylabel,
			 vector<double> &x1, vector<double> &y1,
			 vector<double> &x2, vector<double> &y2,
			 int Startindex, int Endindex, double Mittelwert, double Fehler,
			 bool keep)
{
	string Rohdaten_Name = Dateiname + ".raw.dat";
	string Temp_Skript_Name = Dateiname + ".plt";

	ofstream outfile1, outfile2;
	outfile1.open(Rohdaten_Name.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben
	////////////////////////////////////////////////////////////////////////////
	for (int i = Startindex; i < Endindex; i++) {
		outfile1 << x1[i] << "\t" << y1[i] << "\t" << x2[i] << "\t"
				 << y2[i] << "\n";
	}
	//letzte Zeile
	outfile1 << x1[Endindex] << "\t" << y1[Endindex] << "\t" << x2[Endindex]
			 << "\t" << y2[Endindex];
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben ENDE
	////////////////////////////////////////////////////////////////////////////
	outfile1.close();
	// x und y Grenzen des Datensatzes finden //////////////////////////////
	double x_min = std::min(*std::min_element(x1.begin(), x1.end()),
			*std::min_element(x2.begin(), x2.end()));
	double x_max = std::max(*std::max_element(x1.begin(), x1.end()),
			*std::max_element(x2.begin(), x2.end()));
	double y_min = std::min(*std::min_element(y1.begin(), y1.end()),
			*std::min_element(y2.begin(), y2.end()));
	double y_max = std::max(*std::max_element(y1.begin(), y1.end()),
			*std::max_element(y2.begin(), y2.end()));

	stringstream buf;
	buf.precision(4);
	buf << "scd: " << Mittelwert << " cm^{-2}";
	string text_messwert(buf.str());
	buf.str(string()); // clear the stream
	buf << "error: " << Fehler << " cm^{-2}";
	string text_Fehler(buf.str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile2.open(Temp_Skript_Name.c_str());
	//outfile2<<"set terminal x11\n";   // Betriebssystemabhängig, standard auf
	//was ungefährliches setzen
	outfile2 << "#!/usr/bin/env gnuplot" << endl;
	outfile2 << "set terminal postscript landscape enhanced color "
			 << "font \"Helvetica\" 24 solid linewidth 2\n";
	outfile2 << "set size ratio 0.5\n";
	outfile2 << "set format y \"%.1t{/Symbol \\327}10^{%T}\"\n";
	outfile2 << "set style line 1 lc 1 lt 1 lw 3 pt 7 ps 4\n";
	outfile2 << "set style line 2 lc 3 lt 2 lw 3 pt 5 ps 4\n";
	outfile2 << "set output '" << Dateiname.c_str() << "'\n";
	outfile2 << "set title \'" << title.c_str() << "\'\n";
	outfile2 << "set xlabel \'" << xlabel.c_str() << "\'\n";
	outfile2 << "set ylabel \'" << ylabel.c_str() << "\'\n";
	outfile2 << "set xtics 1\n";
	outfile2 << "set nokey \n"; //keine Legende
	outfile2 << "set label \"" << text_messwert.c_str() << "\" at "
			 << x_min + 0.4 * (x_max - x_min) << ","
			 << y_min + 0.98 * (y_max - y_min) << "\n";
	outfile2 << "set label \"" << text_Fehler.c_str() << "\" at "
			 << x_min + 0.4 * (x_max - x_min) << ","
			 << y_min + 0.9 * (y_max - y_min) << "\n";
	// nun beide Datenreihen mit  Linien Plotten
	outfile2 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2 with lines ls 1, '" << Rohdaten_Name.c_str()
			 << "' using 3:4 with lines ls 2\n";

//    outfile2<<"set terminal postscript landscape enhanced color "
//            <<"\"NimbusSans-Regu\" 28\n";
	//outfile2<<"set terminal epslatex \"NimbusSans-Regu\" 28\n";
	//-> funzt nicht, wies soll
//   outfile2<<"set output '"<<Dateiname.c_str()<<"'\n";
//    outfile2<<"replot\n";
//    outfile2<<"set terminal x11\n";
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript schreiben ENDE
	////////////////////////////////////////////////////////////////////////////
	outfile2.close();
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript ausführen
	////////////////////////////////////////////////////////////////////////////
	string befehl;
	befehl = "gnuplot " + Temp_Skript_Name;
	system(befehl.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript löschen
	////////////////////////////////////////////////////////////////////////////
	// hmm Wartet der, bis Gnuplot fertig ist?---sollte er, in system steckt ja
	// waitpid drin
	remove(Temp_Skript_Name.c_str());
	// remove the plot data file only if requested (the default)
	if (!keep) remove(Rohdaten_Name.c_str());
	return 0;
};
////////////////////////////////////////////////////////////////////////////////
// ENDE Plot_2xy
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plot_Slantcoloumns_polyfit_MgI
////////////////////////////////////////////////////////////////////////////////
int Plot_Slantcoloumns_polyfit_MgI(string Arbeitsverzeichnis, string Dateiname,
								   string title, double *WL, int Startind,
								   int Endind, double *WL_fein,
								   int Startind_fein, int Endind_fein,
								   double *Limb, double *Fitlimb, double *Sonne,
								   double *Fitsonne, double *Quotient,
								   double *Fitquotient)
{
	string Rohdaten_Name = Arbeitsverzeichnis
		+ "/Rohdaten_extrem_unwahrscheinliche_Endung_g"; //geht das parallel???
	string Rohdaten_Name_fein = Arbeitsverzeichnis
		+ "/Rohdaten_fein_extrem_unwahrscheinliche_Endung_g";
	string Temp_Skript_Name = Arbeitsverzeichnis
		+ "/Skript_extrem_unwahrscheinliche_Endung_g";
	ofstream outfile1, outfile2, outfile3;         //3 ist später dazugekommen
	outfile1.open(Rohdaten_Name.c_str());
	outfile2.open(Rohdaten_Name_fein.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben
	////////////////////////////////////////////////////////////////////////////
	for (int i = Startind; i < Endind; i++) {
		outfile1 << WL[i] << "\t" << Limb[i] << "\t" << Sonne[i]
				 << "\t" << Quotient[i] << "\n";
	}
	//letzte Zeile
	outfile1 << WL[Endind] << "\t" << Limb[Endind] << "\t" << Sonne[Endind]
			 << "\t" << Quotient[Endind];
	// Datei mit feinerer Wellenlängenauflösung
	for (int i = Startind_fein; i < Endind_fein; i++) {
		outfile2 << WL_fein[i] << "\t" << Fitlimb[i] << "\t" << Fitsonne[i]
				 << "\t" << Fitquotient[i] << "\n";
	}
	//letzte Zeile
	outfile2 << WL_fein[Endind_fein] << "\t" << Fitlimb[Endind_fein] << "\t"
			 << Fitsonne[Endind_fein] << "\t" << Fitquotient[Endind_fein];
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben ENDE
	////////////////////////////////////////////////////////////////////////////
	outfile1.close();
	outfile2.close();
	////////////////////////////////////////////////////////////////////////////
	// Achsengrenzen finden
	////////////////////////////////////////////////////////////////////////////
	//optional später auch machbar
	////////////////////////////////////////////////////////////////////////////
	// ENDE Achsengrenzen finden
	////////////////////////////////////////////////////////////////////////////
	outfile3.open(Temp_Skript_Name.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile3 << "set style line 1 lc 3 lt 1 lw 3 pt 7 ps 4\n";
	outfile3 << "set style line 2 lc 2 lt 1 lw 3 pt 5 ps 4\n";
	outfile3 << "set style line 3 lc 1 lt 1 lw 3 pt 3 ps 4\n";

	outfile3 << "set terminal postscript landscape enhanced color "
			 << "\"NimbusSans-Regu\" 12\n";
	outfile3 << "set output '" << Dateiname.c_str() << "'\n";
	outfile3 << "set multiplot\n";
	outfile3 << "set size 0.5, 0.33\n"; //# Größe der Einzelnen Plotfenster
	outfile3 << "set pointsize 1\n";
	outfile3 << "set xtics 1.0 \n";
	//########
	//#Plot 1
	//########
	outfile3 << "set origin 0.0, 0.66\n"; //# oben links
	outfile3 << "set title '" << title.c_str()
			 << "' font \"NimbusSans-Regu,10\" \n"; // todo input!!!!
	//set title 'Orbit 12485 Limb TP: Lat 68.298 {/Symbol \260}N Lon 201.956
	// {/Symbol \260}E TH 84.04 km' font \"NimbusSans-Regu,10\"
	outfile3 << "set xlabel 'wavelength[nm]'\n";
	outfile3 << "set ylabel 'l:limb, r:sun'\n";
	outfile3 << "set nokey\n";
	outfile3 << "set mxtics 10\n";
	outfile3 << "set autoscale  y\n";
	outfile3 << "set autoscale y2\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2 with lines ls 1 axes x1y1, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:3 with lines ls 2 axes x1y2\n";
	// ########
	//#Plot 2
	//########
	outfile3 << "set origin 0.5, 0.66;\n"; // # oben rechts
	outfile3 << "set title 'wie links, mit y-Achsen'\n";
	outfile3 << "set xlabel 'wavelength[nm]'\n";
	outfile3 << "set ylabel  \"Limb[phot/(s cm^2 nm)]\" textcolor lt 3\n";
	outfile3 << "set y2label \"Sun[phot/(s cm^2 nm)]\" textcolor lt 2\n";
	outfile3 << "set nokey\n";
	outfile3 << "set mxtics 10\n";
	outfile3 << "set ytics nomirror\n";
	outfile3 << "set y2tics\n";
	outfile3 << "set tics out\n";
	outfile3 << "set autoscale  y\n";
	outfile3 << "set autoscale y2\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2 with lines ls 1 axes x1y1, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:3 with lines ls 2 axes x1y2\n";
	outfile3 << "set nokey\n";
	outfile3 << "set mxtics 10\n";
	outfile3 << "unset y2tics\n";
	outfile3 << "unset y2label\n";
	//########
	//#Plot 3
	//########
	outfile3 << "set origin 0.0, 0.33\n"; //# mitte links
	outfile3 << "set title \"Quotient Messwerte\"\n";
	outfile3 << "set ylabel \"SCD [1/(nm cm^2)]\" textcolor lt 3\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:4 axes x1y1 with lines ls 1\n";
	//########
	//#Plot 4
	//########
	outfile3 << "set title \"Polynomfit Limb\"\n";
	outfile3 << "set ylabel  \"Limb[phot/(s cm^2 nm)]\" textcolor lt 3\n";
	outfile3 << "set origin 0.5, 0.33;\n"; //# mitte rechts
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2  axes x1y1 with points 3 7,'"
			 << Rohdaten_Name_fein.c_str()
			 << "' using 1:2 axes x1y1 with lines ls 3\n";
	//########
	//#Plot 5
	//########
	outfile3 << "set origin 0.0, 0.0;\n"; //# unten links
	outfile3 << "set title \"Quotient Fit\"\n";
	outfile3 << "set ylabel \"SCD [1/(nm cm^2)]\"\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:4 axes x1y1 with lines ls 1, '"
			 << Rohdaten_Name_fein.c_str()
			 << "' using 1:4 axes x1y1 with lines ls 2\n";
	//########
	//#Plot 6
	//########
	outfile3 << "set origin 0.5, 0.0;\n"; //# unten rechts
	outfile3 << "set title \"Polynomfit Sun\"\n";
	outfile3 << "set ylabel  \"Sun[phot/(s cm^2 nm)]\" textcolor lt 3\n";
	//outfile3<<"set yrange [0:5000]\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:3  axes x1y1 with points 3 7,'"
			 << Rohdaten_Name_fein.c_str()
			 << "' using 1:3 axes x1y1 with lines ls 3\n";
	outfile3 << "unset multiplot\n";
	outfile3 << "set output\n";
	outfile3 << "reset\n";
	////////////////////////////////////////////////////////////////////////////
	//ENDE Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile3.close();
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript ausführen
	////////////////////////////////////////////////////////////////////////////
	string befehl;
	befehl = "gnuplot " + Temp_Skript_Name;
	system(befehl.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript löschen
	////////////////////////////////////////////////////////////////////////////
	remove(Temp_Skript_Name.c_str());
	remove(Rohdaten_Name.c_str());
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Plot_Slantcoloumns_polyfit
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plot_Spektren_und_Quotient
////////////////////////////////////////////////////////////////////////////////
int Plot_Spektren_und_Quotient(string Arbeitsverzeichnis, string Dateiname,
							   string title, double *WL, int Startind,
							   int Endind, double *Limb, double *Limb_error,
							   double *Sonne, double *Quotient,
							   double *Quotient_error)
{
	string Rohdaten_Name = Arbeitsverzeichnis
		+ "/Rohdaten_extrem_unwahrscheinliche_Endung_g"; //geht das parallel???
	string Temp_Skript_Name = Arbeitsverzeichnis
		+ "/Skript_extrem_unwahrscheinliche_Endung_g";
	ofstream outfile1, outfile3;          //3 ist später dazugekommen
	outfile1.open(Rohdaten_Name.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben
	////////////////////////////////////////////////////////////////////////////
	for (int i = Startind; i < Endind; i++) {
		outfile1 << WL[i] << "\t" << Limb[i] << "\t" << Limb_error[i] << "\t"
				 << Sonne[i] << "\t" << Quotient[i] << "\t" << Quotient_error[i]
				 << "\n";
	}
	//letzte Zeile
	outfile1 << WL[Endind] << "\t" << Limb[Endind] << "\t" << Limb_error[Endind]
			 << "\t" << Sonne[Endind]
			 << "\t" << Quotient[Endind] << "\t" << Quotient_error[Endind];
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben ENDE
	////////////////////////////////////////////////////////////////////////////
	outfile1.close();
	outfile3.open(Temp_Skript_Name.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile3 << "set style line 1 lc 3 lt 1 lw 3 pt 7 ps 4\n";
	outfile3 << "set style line 2 lc 2 lt 1 lw 3 pt 5 ps 4\n";
	outfile3 << "set style line 3 lc 1 lt 1 lw 3 pt 3 ps 4\n";
	outfile3 << "set style line 4 lc 0 lt 1 lw 3 pt 1 ps 1\n";

	outfile3 << "set terminal postscript landscape enhanced color "
			 << "\"NimbusSans-Regu\" 12\n";
	outfile3 << "set output '" << Dateiname.c_str() << "'\n";
	outfile3 << "set multiplot\n";
	outfile3 << "set size 1, 0.5\n"; //# Größe der Einzelnen Plotfenster
	outfile3 << "set pointsize 1\n";
	outfile3 << "set xtics 2.0 \n";
	//########
	//#Plot 1
	//########
	outfile3 << "set origin 0.0, 0.5\n"; //# oben links
	outfile3 << "set title '" << title.c_str()
			 << "' font \"NimbusSans-Regu,10\" \n"; // todo input!!!!
	//set title 'Orbit 12485 Limb TP: Lat 68.298 {/Symbol \260}N Lon 201.956
	//{/Symbol \260}E TH 84.04 km' font \"NimbusSans-Regu,10\"
	outfile3 << "set xlabel 'wavelength[nm]'\n";
	outfile3 << "set ylabel 'l:limb, r:sun'\n";
	outfile3 << "set nokey\n";
	outfile3 << "set mxtics 10\n";
	outfile3 << "set autoscale  y\n";
	outfile3 << "set autoscale y2\n";
	//outfile3<<"set y2label \"Sun[phot/(s cm^2 nm)]\" textcolor lt 2\n";
	//outfile3<<"set ytics nomirror\n";
	//outfile3<<"set y2tics\n";
	//outfile3<<"set tics out\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2 axes x1y1 with lines ls 1, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:4 axes x1y2 with lines ls 2, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:2:3 axes x1y1 with yerrorbars ls 4\n";
	//outfile3<<"unset y2tics\n";
	//outfile3<<"unset y2label\n";
	//########
	//#Plot 2
	//########
	outfile3 << "set origin 0.0, 0.0\n"; //# mitte links
	outfile3 << "set title \"Quotient Messwerte\"\n";
	outfile3 << "set ylabel \"SCD [1/(nm cm^2)]\" textcolor lt 3\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:5 axes x1y1 with lines ls 1, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:5:6 axes x1y1 with yerrorbars ls 4\n";

	outfile3 << "unset multiplot\n";
	outfile3 << "set output\n";
	outfile3 << "reset\n";
	////////////////////////////////////////////////////////////////////////////
	//ENDE Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile3.close();
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript ausführen
	////////////////////////////////////////////////////////////////////////////
	string befehl;
	befehl = "gnuplot " + Temp_Skript_Name;
	system(befehl.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript löschen
	////////////////////////////////////////////////////////////////////////////
	remove(Temp_Skript_Name.c_str());
	remove(Rohdaten_Name.c_str());
	return 0;

}
////////////////////////////////////////////////////////////////////////////////
// ENDE Plot_Spektren_und_Quotient
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plot_Quotient_mit_Fehler
////////////////////////////////////////////////////////////////////////////////
int Plot_Quotient_mit_Fehler(string Arbeitsverzeichnis, string Dateiname,
							 string title, double *WL, int Startind, int Endind,
							 double *Quotient, double *Quotient_error,
							 double *Quotient_stabw)
{
	string Rohdaten_Name = Arbeitsverzeichnis
		+ "/Rohdaten_extrem_unwahrscheinliche_Endung_g"; //geht das parallel.???
	string Temp_Skript_Name = Arbeitsverzeichnis
		+ "/Skript_extrem_unwahrscheinliche_Endung_g";
	ofstream outfile1, outfile3;          //3 ist später dazugekommen
	outfile1.open(Rohdaten_Name.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben
	////////////////////////////////////////////////////////////////////////////
	for (int i = Startind; i < Endind; i++) {
		outfile1 << WL[i] << "\t" << Quotient[i] << "\t" << Quotient_error[i]
				 << "\t" << Quotient_stabw[i] << "\n";
	}
	//letzte Zeile
	outfile1 << WL[Endind] << "\t" << Quotient[Endind] << "\t"
			<< Quotient_error[Endind] << "\t" << Quotient_stabw[Endind] << "\n";
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript Rohdatenfile schreiben ENDE
	////////////////////////////////////////////////////////////////////////////
	outfile1.close();
	outfile3.open(Temp_Skript_Name.c_str());

	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile3 << "set style line 1 lc 3 lt 1 lw 3 pt 7 ps 4\n";
	outfile3 << "set style line 2 lc 2 lt 1 lw 3 pt 5 ps 4\n";
	outfile3 << "set style line 3 lc 1 lt 1 lw 3 pt 3 ps 4\n";
	outfile3 << "set style line 4 lc 0 lt 1 lw 3 pt 1 ps 1\n";

	outfile3 << "set terminal postscript landscape enhanced color "
			 << "\"NimbusSans-Regu\" 12\n";
	outfile3 << "set output '" << Dateiname.c_str() << "'\n";
	outfile3 << "set multiplot\n";
	outfile3 << "set size 1, 0.5\n"; //# Größe der Einzelnen Plotfenster
	outfile3 << "set pointsize 1\n";
	outfile3 << "set xtics 2.0 \n";
	//########
	//#Plot 1
	//########
	outfile3 << "set origin 0.0, 0.5\n"; //# oben links
	//outfile3<<"set title '"<<title.c_str()<<"' font \"NimbusSans-Regu,10\" \n";
	// todo input!!!!
	outfile3 << "set title \"" << title.c_str() << "\" \n"; // todo input!!!!
	outfile3 << "set x2label \"Quotient Messwerte mit Fehlern\"\n";
	outfile3 << "set xlabel 'wavelength[nm]'\n";
	outfile3 << "set nokey\n";
	outfile3 << "set mxtics 10\n";
	outfile3 << "set autoscale  y\n";

	outfile3 << "set ylabel \"SCD [1/(nm cm^2)]\" textcolor lt 3\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2 axes x1y1 with lines ls 1, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:2:3 axes x1y1 with yerrorbars ls 4\n";

	//########
	//#Plot 2
	//########
	outfile3 << "unset x2label \n";
	outfile3 << "set origin 0.0, 0.0\n"; //# mitte links
	outfile3 << "set title \"Quotient Messwerte mit Standardabweichung der Residuen\"\n";
	outfile3 << "set ylabel \"SCD [1/(nm cm^2)]\" textcolor lt 3\n";
	outfile3 << "plot '" << Rohdaten_Name.c_str()
			 << "' using 1:2 axes x1y1 with lines ls 1, '"
			 << Rohdaten_Name.c_str()
			 << "' using 1:2:4 axes x1y1 with yerrorbars ls 4\n";

	outfile3 << "unset multiplot\n";
	outfile3 << "set output\n";
	outfile3 << "reset\n";
	////////////////////////////////////////////////////////////////////////////
	//ENDE Gnuplotscript schreiben
	////////////////////////////////////////////////////////////////////////////
	outfile3.close();
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript ausführen
	////////////////////////////////////////////////////////////////////////////
	string befehl;
	befehl = "gnuplot " + Temp_Skript_Name;
	system(befehl.c_str());
	////////////////////////////////////////////////////////////////////////////
	// Gnuplotscript löschen
	////////////////////////////////////////////////////////////////////////////
	remove(Temp_Skript_Name.c_str());
	remove(Rohdaten_Name.c_str());
	return 0;

}
////////////////////////////////////////////////////////////////////////////////
// ENDE Plot_Quotient_mit_Fehler
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Plots_Zusammenfassen
////////////////////////////////////////////////////////////////////////////////
int Plots_Zusammenfassen(string Pfad_multips2pdf, string Pfad_multips2ps,
		string Name_pdf_Datei, vector<string> Liste_der_ps_Dateinamen)
{
	// Die Plots in der Liste mit den ps Dateien werden in eine PDF Datei
	// zusammengefasst und gelöscht Das skriptfile multips2pdf a.pdf b.ps...z.ps
	// usw Sorgt für die Umwandlung...es muss also lediglich die Kommandozeile
	// geschrieben werden Es können wohl nicht beliebig viele
	// Kommandozeilenparameter übergeben werden...also probieren wir mal 512=2^8
	// Da noch noch die exe und der aufruf als parameter da sind muss weniger
	// als 2^8 genommen werden
	unsigned int M = 200;
	vector<string> Liste_der_grossen_pdf;
	vector<string>::iterator sit;
	string Befehlszeile;
	unsigned int Max_Zahl_grosse_pdf = Liste_der_ps_Dateinamen.size() / M;
	//Achtung integerdivision ist absicht
	if ((Liste_der_ps_Dateinamen.size() % M) != 0) {
		Max_Zahl_grosse_pdf++;
	}
	// Zuerst viele ps in wenigen grossen pdf zusammenfassen
	cout << "Liste_der_ps_Dateinamen.size(): "
		 << Liste_der_ps_Dateinamen.size() << "\n";
	cout << "Max_Zahl_grosse_pdf: " << Max_Zahl_grosse_pdf << "\n";
	cout << "Erzeuge_große_PDF\n";
	for (unsigned int k = 0; k < Max_Zahl_grosse_pdf; k++) {
		cout << k << "te große ps Datei\n";
		stringstream buf;
		buf << Name_pdf_Datei << k << ".pdf";
		string Name_grosse_pdf(buf.str());
		Liste_der_grossen_pdf.push_back(Name_grosse_pdf);
		Befehlszeile = Pfad_multips2pdf + " " + Name_grosse_pdf;
		for (unsigned int i = 0;
			  (i < M) && ((i + k * M) < Liste_der_ps_Dateinamen.size()); i++) {
			//statt size mal potenzen von 2 Probieren
			Befehlszeile += " " + Liste_der_ps_Dateinamen[i + k * M];
		}
		system(Befehlszeile.c_str());
	}
	cout << "Erzeuge große pdf\n";
	// die großen pdf zu einer pdf zusammenfassen
	Befehlszeile = Pfad_multips2pdf + " " + Name_pdf_Datei;
	for (sit = Liste_der_grossen_pdf.begin();
			sit != Liste_der_grossen_pdf.end(); ++sit ) {
		Befehlszeile += " " + *sit;
	}
	//cout<<Liste_der_ps_Dateinamen.size()<<"\n";
	system(Befehlszeile.c_str());
	//Die ps sind jetzt im pdf drin, können also gelöscht werden
	cout << "lösche die ps\n";
	for (sit = Liste_der_ps_Dateinamen.begin();
			sit != Liste_der_ps_Dateinamen.end(); ++sit ) {
		remove(sit->c_str());
	}
	for (sit = Liste_der_grossen_pdf.begin();
			sit != Liste_der_grossen_pdf.end(); ++sit ) {
		remove(sit->c_str());
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// ENDE Plots_Zusammenfassen
////////////////////////////////////////////////////////////////////////////////
