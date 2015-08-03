/*
 * Messung_Nadir.cpp
 *
 *  Created on: 30.04.2010
 *      Author: martin
 */

#include <cstdio>
#include "Messung_Nadir.h"
#include <cmath>
#include <cstdlib>
#include "Ausgewertete_Messung_Nadir.h"
#include <cstdio>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include "Ausdrucke.h"
#include "Speziesfenster.h"
#include "Dateinamensteile_Bestimmen.h"

using std::cout;
using std::string;
using std::vector;
//////////////////////////////////////////////////
//constructor
//////////////////////////////////////////////////
Messung_Nadir::Messung_Nadir(std::string filename) :
	Messung(filename) {}
//////////////////////////////////////////////////
//constructor ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Methoden
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////
//Zeilendichte_Bestimmen
//////////////////////////////////////////////////
void Messung_Nadir::Zeilendichte_Bestimmen(Speziesfenster &Spezfenst, int Index,
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
	// beinhalten darf.
	// Bis Hierhin haben wir I/(piFGamma) berechnet. Nun wollen wir die
	// "normierte" Intensität für einen Peak berechnen Da der Peak gerade im
	// Kanal zumeist auf einem relativ großen Untergrundsignal sitzt, wird
	// dieses zunächst abgezogen.
	// Dafür wird im Bereich um die Linie herum, wo keine signifikanten anderen
	// Linien liegen(Basisfenster) ein linearer Fit der Grundlinie
	// durchgeführt, welcher vom der "normierten" Intensität abgezogen wird.
	// Danach wird im Bereich des Peakfensters die Fläche des Peaks bestimmt.
	// Da die Linie zumeist nur aus 3 Punkten besteht,  hat das Profil im
	// wesentlichen die Form der Auflösungsfunktion, welche eine hyperboloide
	// 4ten Grades ist. Diese hat im wesentlichen 3 Parameter: Wellenlänge des
	// Peaks, Breite und Fläche.
	// Die robusteste Variante ist es, die Wellenlänge und die Halbwertsbreite
	// der Linie vorzugeben und nur die Fläche anzufitten.  Diese Variante kann
	// durch einen linearen Parameterfit erreicht werden und ist zum einen
	// schnell und was hier viel wichtiger ist
	// robust. Andere nichtlineare Verfahren sind entweder nicht robust
	// (Gauss Newton, insbesondere bei dem schlechten Signal/Noise Verhältnis)
	// oder deutlich langsamer (simulated annealing bzw. ist das aufwändiger
	// dies noch schnell zu implementieren... das wäre Schritt 2)
	// Die ermittelte Fläche ist dann unsere gesuchte Säulenzeilendichte.

	// I/(piFGamma)=integral(AMF n ds) mit AMF = s exp(-tau)
	// ...aber zu der Formel später nochmal zurück
	// Das spätere Retrieval ermittelt dann die Dichte n aus der rechten Seite

	//Zunächst Indizes der Wellenlaengen der Basisfensterbestimmen
	int Index_Basisfenster_links_min = Get_Index(Spezfenst.m_Basisfenster_links_WLmin[Index]);
	int Index_Basisfenster_links_max = Get_Index(Spezfenst.m_Basisfenster_links_WLmax[Index]);
	int Index_Basisfenster_rechts_min = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmin[Index]);
	int Index_Basisfenster_rechts_max = Get_Index(Spezfenst.m_Basisfenster_rechts_WLmax[Index]);
	int Index_Peakfenster_min = Get_Index(Spezfenst.m_Peakfenster_WLmin[Index]);
	int Index_Peakfenster_max = Get_Index(Spezfenst.m_Peakfenster_WLmax[Index]);
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
			this->m_Intensitaeten_durch_piF_Gamma[Index_Basisfenster_links_min + i];
	}
	for (int i = 0; i < Bas_r; i++) {
		Basisfenster_WL[Bas_l + i] = this->m_Wellenlaengen[Index_Basisfenster_rechts_min + i];
		Basisfenster_Intensitaet[Bas_l + i] =
			this->m_Intensitaeten_durch_piF_Gamma[Index_Basisfenster_rechts_min + i];
	}
	//Peakfenster WL und I auffüllen
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_WL[i] = m_Wellenlaengen[Index_Peakfenster_min + i];
		Peakfenster_Intensitaet[i] = m_Intensitaeten_durch_piF_Gamma[Index_Peakfenster_min + i];
	}
	// linearen Fit des Basisfensters durchführen
	// Proto: Fit_Linear(double* x,double* y, double& a0, double& a1,int Anfangsindex, int Endindex)
	double a0, a1, rms_err_base;
	Fit_Linear(Basisfenster_WL, Basisfenster_Intensitaet, a0, a1, rms_err_base,
			0, N_Basis - 1);
	// lineare Funktion von Intensitäten des Peakfenster abziehen
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_Intensitaet[i] -= a0 + a1 * Peakfenster_WL[i];
	}
	// Hyperboloiden an Peakfenster anfitten
	// Proto:
	// Fit_Peak_hyperbolic(double* x,double* y,double x0, double FWHM,
	//   double& A, int Anfangsindex, int Endindex)
	double Flaeche;
	Fit_Peak_hyperbolic(Peakfenster_WL, Peakfenster_Intensitaet,
						Spezfenst.m_Wellenlaengen[Index],
						Spezfenst.m_FWHM, Flaeche, 0, N_Peak - 1);
	//Fehler des Fits bestimmen...
	//da Peakfenster_Intensitaet nicht mehr gebraucht wird
	//die Basislinie für die Fehlerberechnung wieder aufaddiert
	m_Zeilendichte = Flaeche;
	for (int i = 0; i < N_Peak; i++) {
		Peakfenster_Intensitaet[i] += a0 + a1 * Peakfenster_WL[i];
	}
	// Funktion double Messung_Limb::Evaluate_Error_primitive(double* x,
	//   double* y, double a0,double a1, double A, double FWHM, double x0,
	//   int Anfangsindex, int Endindex)
	m_Fehler_Zeilendichten =
		Evaluate_Error_primitive(Peakfenster_WL, Peakfenster_Intensitaet,
				a0, a1, m_Zeilendichte, Spezfenst.m_FWHM,
				Spezfenst.m_Wellenlaengen[Index], 0, N_Peak - 1);
	////////////////////////////////////////////////////////////////////////////
	// Hier kann man zur Testzwecken noch einen Plot machen  ///////////////////
	if (mache_Fit_Plots == "ja" && Spezfenst.plot_fit) {
		//TODO das als Funktion implementieren
		vector<double> Funktion;

		for (int i = 0; i < N_Peak; i++) {
			double Basis = a0 + a1 * Peakfenster_WL[i];
			double Peak = m_Zeilendichte *
				slit_func(Spezfenst.m_FWHM,
						Spezfenst.m_Wellenlaengen[Index], Peakfenster_WL[i]);
			//cout<<Peak<<"\n";
			//cout<<m_Zeilendichte<<"\n";
			Funktion.push_back(Peak + Basis);
		}

		string s_OrbNum;
		std::stringstream buf;
		//TODO immer prüfen, ob Dateienamenlänge noch stimmt...
		// falls / im Namen ist das schlecht
		string Datnam = sb_basename(m_Dateiname_L1C);

		//TODO Pfad anpassen
		string plot_dir = Arbeitsverzeichnis + "/Plots";
		mkdir(plot_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		buf << Datnam.c_str() << "_" << Spezfenst.m_Spezies_Name.c_str()
			<< "_" << m_Messung_ID << "_" << Index << ".ps";
		string new_datnam(buf.str());
		string s1(plot_dir + "/" + new_datnam);
		//s1 ist der Volle Pfad der Datei...
		//diesen kann man wegspeichern,
		//um später die .ps files in ein großes pdf zu packen
		Spezfenst.m_Liste_der_Plot_Dateinamen.push_back(s1);
		//Orbitnummer ermitteln/////
		// egal, wie die Datei heißt...
		// die Orbitnummer sind die 5 Zeichen vor .dat
		size_t pos_suffix = 0;
		pos_suffix = Datnam.find(".dat");
		if (pos_suffix == string::npos) {
			cout << " kein .dat in Nadirdateiname...Orbitnummer nicht findbar\n";
			s_OrbNum = "xxxxx";
		} else {
			s_OrbNum = Datnam.substr(pos_suffix - 5, 5);
		}
		//Orbitnummer ermittelt///////
		buf.str(string());
		buf << "Orbit " << s_OrbNum.c_str() << ", Nadir GP:"
			<< " Lat: " << m_Latitude_Ground << " deg,"
			<< " Lon: " << m_Longitude_Ground << " deg; Sat:"
			<< " Lat: " << m_Latitude_Sat << " deg,"
			<< " Lon: " << m_Longitude_Sat << " deg.";
		string s2(buf.str());
		//cout<<s1<<"\n";
		//int Plot_2xy(string Dateiname,string title, string xlabel,
		//  string ylabel,double* x1,double*y1, double* x2,double* y2,
		//  int Startindex,int Endindex);
		//Plot_2xy(s1.c_str(),s1.substr(s1.size()-50,50).c_str(),"$\\lambda$ in nm",
		// "$\\frac{I}{\\piF\\gamma}$",Peakfenster_WL,Peakfenster_Intensitaet,
		// Peakfenster_WL,Funktion,0,N_Peak-1); //-> Fit geht
		Plot_2xy(Arbeitsverzeichnis.c_str(), s1.c_str(), s2.c_str(),
				"Wellenlaenge in nm",
				"Schraege Saeule bei Peakposition in cm^{-2}/nm",
				 Peakfenster_WL, Peakfenster_Intensitaet, Peakfenster_WL,
				 Funktion, 0, N_Peak - 1,
				 m_Zeilendichte, m_Fehler_Zeilendichten);
	}
	// Ende Plot ///////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
}//int Zeilendichte_Bestimmen() ende
//////////////////////////////////////////////////
//Zeilendichte_Bestimmen ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//Ergebnis_Zusammenfassen
//////////////////////////////////////////////////
Ausgewertete_Messung_Nadir  Messung_Nadir::Ergebnis_Zusammenfassen()
{
	Ausgewertete_Messung_Nadir aus;
	// Nullinitialisierung...
	// die Wellenlänge des Übergangs steckt im Speziesfenster
	aus.m_Wellenlaenge = 0;
	//Ergebnisse
	aus.m_Zeilendichte = this->m_Zeilendichte;
	aus.m_Fehler_Zeilendichten = this->m_Fehler_Zeilendichten;
	// Zwischenergebnisse
	aus.m_Deklination = this->m_Deklinationswinkel;
	aus.m_Sonnen_Longitude = this->m_Sonnen_Longitude;
	// Datum
	aus.m_Jahr = this->m_Jahr;
	aus.m_Monat = this->m_Monat;
	aus.m_Tag = this->m_Tag;
	aus.m_Stunde = this->m_Stunde;
	aus.m_Minute = this->m_Minute;
	// Geolokationen
	aus.m_Latitude_Sat = this->m_Latitude_Sat;
	aus.m_Longitude_Sat = this->m_Longitude_Sat;
	aus.m_Hoehe_Sat = this->m_Hoehe_Sat;
	aus.m_Latitude_Ground = this->m_Latitude_Ground;
	aus.m_Longitude_Ground = this->m_Longitude_Ground;
	aus.m_Erdradius = this->m_Erdradius;
	return aus;
}//Ausgewertete_Messung_Nadir  Ergebnis_Zusammenfassen() ENDE
//////////////////////////////////////////////////
//Ergebnis_Zusammenfassen ENDE
//////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Methoden ENDE
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//
// Wartungsfunktionen
//
////////////////////////////////////////////////////////////////////////////////
void Messung_Nadir::Ausgabe_in_Datei(string Dateiname)
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
	//checken, ob Datei auch offen fehlt...
	//aber ok Funktion wird eh beim debuggen eingesetzt...
	//da kriegt man das schon raus..hoffentlich
	//Datei schreiben
	///////////////////////////////////////////////////////////
	//zunächst die randdaten in den header
	// Ergebnisse
	fprintf(outfile, "%12s %1.5E\n", "m_Zeilendichte: ", m_Zeilendichte);
	fprintf(outfile, "%12s %1.5E\n", "m_Fehler_Zeilendichten: ", m_Fehler_Zeilendichten);
	// Zwischenergebnisse
	fprintf(outfile, "%12s %1.5E\n", "m_Deklinationswinkel: ", m_Deklinationswinkel);
	fprintf(outfile, "%20s %1.5E\n", "m_Sonnen_Longitude: ", m_Sonnen_Longitude);
	// Herkunftsmerkmale
	fprintf(outfile, "%12s %20s\n", "m_Dateiname_L1C: ", m_Dateiname_L1C.c_str());
	fprintf(outfile, "%12s %i\n", "m_Messung_ID: ", m_Messung_ID);
	//Datum
	fprintf(outfile, "%12s %i\n", "m_Jahr: ", m_Jahr);
	fprintf(outfile, "%12s %i\n", "m_Monat: ", m_Monat);
	fprintf(outfile, "%12s %i\n", "m_Tag: ", m_Tag);
	fprintf(outfile, "%12s %i\n", "m_Stunde: ", m_Stunde);
	fprintf(outfile, "%12s %i\n", "m_Minute: ", m_Minute);
	// Geolokationen
	fprintf(outfile, "%12s %1.5E\n", "m_Latitude_Sat: ", m_Latitude_Sat);
	fprintf(outfile, "%12s %1.5E\n", "m_Longitude_Sat: ", m_Longitude_Sat);
	fprintf(outfile, "%12s %1.5E\n", "m_Hoehe_Sat: ", m_Hoehe_Sat);
	fprintf(outfile, "%12s %1.5E\n", "m_Latitude_Ground: ", m_Latitude_Ground);
	fprintf(outfile, "%12s %1.5E\n", "m_Longitude_Ground: ", m_Longitude_Ground);
	fprintf(outfile, "%12s %1.5E\n", "m_Erdradius: ", m_Erdradius);
	//Füllbare Felder
	fprintf(outfile, "%12s %i\n", "m_Number_of_Wavelength: ", m_Number_of_Wavelength);
	//Überschrift
	fprintf(outfile, "%12s %12s %12s %12s\n",
			"m_Wellenlaengen", "m_Intensitaeten",
			"m_Intensitaeten_durch_piF", "m_Intensitaeten_durch_piF_Gamma");
	//großes Feld
	for (int i = 0; i < m_Number_of_Wavelength; i++) {
		//die letzte Zeile der Datei ist so leer, das \n in der Vorletzten steht
		fprintf(outfile, "%1.5E %1.5E %1.5E %1.5E\n",
				m_Wellenlaengen[i], m_Intensitaeten[i],
				m_Intensitaeten_durch_piF[i], m_Intensitaeten_durch_piF_Gamma[i]);
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

