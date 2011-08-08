/*******************************************************************************
main.cpp

Dieses Programm ist eine C++ Version des Retrievalprogramms SCIA2D.f90,
welches von Marco Scharringhausen entwickelt wurde.
Zur Vergleichbarkeit sind die durchgeführten Rechnungen gleich, lediglich
das Layout des Programms soll etwas anders sein.
( Das stimmt schon garnicht mehr....
  Streuwinkel werden anders berechnet, SZA auch, Das Raytracing ist 3D, d.h. die
  längenänderung wird für die Schritte berücksichtigt
  usw.....am Ende kommt aber dasselbe raus, d.h. die kleinen Bugs, die vorher
  drin waren waren nicht sehr gravierend
  Bei stärkerer Absorption wäre das aber schon der Fall)


Ich benutzte hier als integrierte Entwicklungsumgebung eclipse (kostenlos
erhältlich, für windows und linux) Auf Windows gibt es auch kostenlose
Versionen des Visual Studios (Express Versionen....da ham einen
verständlicheren Debugger) Natürlich kann man auch andere Editoren wie Kate,
oder vi verwenden, je nach Geschmack

Ein Ziel dieses Programms ist es, schnell zu sein. Ein anderes ist
Flexibilität, was in diesem Fall die möglichst schnelle und einfache
Anpassungfähigkeit des Programms an neue Spezies bedeutet.  Letzteres soll
dadurch erreicht werden, dass diese wichtigen Veränderung nur im Hauptprogramm
durchgeführt werden müssen.  ( Das hat auch geklappt, nur bei den Spezies in
Limbauswertung und Nadirauswertung, müssen einige Teile ergänzt werden, die
aber gleich zu den anderen aussehen) Weitere Teile, die ich oft ändere sind das
Retrievalgitter In Retrievaliteration, kann man die maximal Iterationszahl, bzw
den treshold heruntersetzen, wenns zu lange dauert

 Um die kritischen Stellen leicht zu finden, soll der Code im Hauptprogramm
 entsprechend kurz sein, das heißt, dass längere Passagen in möglichst klar
 benannten Funktionen untergebracht sind, sodass im Hauptprogramm auch die
 grobe Struktur der Lösung des Problems schnell nachvollzogen werden kann.

Schnelligkeit ist z.B. dadurch zu erreichen, dass Dateien, die mehrere
Datensätze enthalten, auch wenn dies nicht der  intuitiven Abarbeitungslogik
entspricht, trotzdem nur einmal geladen, da dies sehr lange dauert.  Das dauert
vor allem lange, wenn die Dateien nicht auf der richtigen Platte liegen, also
lokal speichern

Es sollte auch unbedingt vermieden werden, zu viel mit den Großen Feldern, die
aus den Dateien geladen werden "herumzujonglieren"(d.h. möglichst die großen
Felder nicht kopieren).....Das ist doch nicht so schlimm, das Laden so wie es
jetzt ist, dauert höchstens eine Sekunde, da kann man flexibler sein

Das Programm sollte zudem möglichst parallelisierbar sein, d.h. Gesamtprobleme
sollten so in Teilprobleme aufgeteilt werden, dass zur Lösung der Teilprobleme
die anderen gleichartigen Teilprobleme nicht bekannt sein müssen.  Das Programm
wertet nur für einen Orbit aus.......  Die einfachste Art der Parallelisierung,
ist mehrere Orbits parallel laufen zu lassen.  Dadurch muss am Programm selbst
nicht so viel(eigentlich garnix) geändert werden. Ein umhüllendes Programm kann
dann geschrieben werden, welches die Orbits in gleich große Listen teilt, und
jeweils für teillisten seriell das Programm startet (siehe shell befehle fork,
exec, waitpid, system  , für kompliziertere parallelisierungen mpi u.ä.
verwenden (ist aber glaube ich gar nicht nötig)

( Die veränderlichen Parameter in der Konf Datei müssen dem Programm Übergebbar
sein) Die Konf Datei ist ein Relikt aus der Umsetzung von Marcos Programm und
muss immer im Ordner sein, wo auch main.cpp liegt) Einträge dort werden aber
zum Teil überschrieben

Die Geschwindigkeitsoptimierung sollte nur die Flaschenhälse betreffen
Die Dateien auf der eigenen Platte zu haben ist für die Geschwindigkeit
bis jetzt (september 2010) der wichtigste Geschwindigkeitsfaktor damit wurde
auf ein BINÄRDATEIENSYSTEM UMGESTIEGEN (einfache Umwandlungsfunktionen hab ich
auch-> mail)
 - unglücklicherweise habe ich Schrägen Säulendichten anfangs Zeilendichten
   genannt...ist also hier dasselbe
 - Das Raytracing in Teil 3 dauert recht lange. Da die Absorption Quasi keine
   Rolle spielt, könnte man auch nur die Lichtwege nehmen.  Die kann man
   ausrechnen und braucht kein Raytracing
 - Die Matrix selbst sollte auch nicht Luftmassenfaktoren Matrix oder Air
   mass faktor Matrix heißen, weil das eigentlich ein andere Sachverhalt ist
   (Quotient aus Schräger SÄule und Vertikaler Säule, oder so)

Dort, wo sich die Flexibilität mit der Geschwindigkeit beißt, muss abgewogen
werden, was wichtiger ist
 //TODO Matrixmultiplikationen überprüfen, ob die Matrizen auch in einer
 //günstigen Reihenfolge Multipliziert werden
 //TODO Matrixoperationen sind z.t. für große Matrizen langsam, für 200*200
 //Matrizen ist das Retrieval aber noch schnell Abhilfe schaffen hier
 //Matrizen der LAPACK Bibliothek....da hab ich auch schon ein Programm was
 //LAPACK und ATLAS(BLAS)
 // benutzt und für eine 4000*4000 Matrix die LU-Zerlegung in 5 sekunden
 // schafft (Matlab braucht da 8 für A\x)

 //Sonnenzenitwinkel und Streuwinkel werden für jeden Gitterpunkt berechnet,
 //und nicht für jeden Messschritt Da die Winkel innerhalb eines Gitterpunktes
 //nahezu konstant sind und unterschiede bei der späteren verwendung in sin und
 //cos eher noch kontrahiert werden...also vernachlässigbar sind....das
 //betrifft Teil 3...wo die Absorption wie gesagt, eh fast zu vernachlässigen
 //sit


Für überzeugende Plots, ob das Raytracing auch funktioniert ist es besser TAU
zu plotten als AMF Tau steigt nämlich mit zunehmendem Weg Monoton an...(bereits
implementiert)

Damit das Programm läuft, müssen einige Bibliotheken mit eingebunden werden
(am Besten in der Reihenfolge, also gfortran zuerst):
LimbNadir_IO
lapack
cblas
f77blas
atlas
gfortran

Der Haken ist, dass man ATLAS auf seinem eigenen Rechner installieren muss (was
bei Linux Programmen eine ziemliche Qual ist. Es gibt aber ein Manual
auf deren Seite dem man folgen kann und wenn man Glück hat, dann
funktioniert das sogar), da das auf die eigene Hardware optimiert, das
funktioniert auch nur unter Linux, noch schneller wäre intel mkl, kostet aber
geld

Ein letzter Hinweis,: nicht alle Kommentare in dem Quelltext müssen noch
aktuell sein, da so ein Programm ständigen Änderungen unterworfen ist Falls
irgendetwas nicht funktioniert, mich kontaktieren

*******************************************************************************/
//eingebundene Headerfiles
//standard
#include <string>           // string
#include <cstdio>         // cout
#include <iostream>      // cout
#include <cstdlib>        // für atoi
#include <sstream>
//eigene
#include "Konfiguration.h"
#include "Liniendaten.h"
#include "Messung_Limb.h"
#include "Datei_IO.h"
#include "Speziesfenster.h"
#include "Orbitliste.h"
#include "Sonnenspektrum.h"

#include "Limbauswertung.h"
#include "Nadirauswertung.h"
#include "Nachricht_Schreiben.h"
#include <ctime>  //langsames timing...so auf Sekunden genau
#include "Retrievalgitter.h"
#include "MPL_Matrix.h"
#include "Matrizen_Aufbauen.h"
#include "Retrievaliteration.h"
#include "Retrievalfehler_Abschaetzung.h"
#include "Ausdrucke.h"                              // Plots_Zusammenfassen
#include "Glaetten_2D.h"
#include "NO_emiss.h"

#include "Dateinamensteile_Bestimmen.h"
#include "Ausgewertete_Messung_Limb.h"
#include "Ausgewertete_Messung_Nadir.h"
//===========================================================
// eigene externe Bibliotheken-> Müssen mitgeliefert werden,
// falls nicht dabei, mail an martin
// LimbNadir_IO // Bibliothek zum laden und speichern von Scia L1C Daten in
//              // Ascii oder martins binärformat
//              // wird beim einlesen der Rohdaten in Limbauswertung und
//              // Nadirauswertung verwendet
//              // und macht so die Funktionenen in Datei_IO.h z.T. überflüssig
//===========================================================
// den std namespace verwenden damit nicht immer std::
// geschrieben werden muss
using namespace std;

///////////////////////////////////////
// START globale Variablen
//Nachrichtenprioritäten
int Prioritylevel = 0;
// ENDE globale Variablen
///////////////////////////////////////
//
// the SCIAMACHY slit function
double slit_func(double fwhm, double x0, double x)
{
	const double fwhm2 = fwhm * fwhm;
	const double cnorm = 4. * M_PI * M_SQRT2 / (fwhm2 * fwhm);
	// (0.5 * FWHM)^4
	const double fwhm2to4 = 0.0625 * fwhm2 * fwhm2;

	return 1. / (cnorm * (fwhm2to4 + pow(x0 - x, 4)));
}
double slit_func_gauss(double fwhm, double x0, double x)
{
	double s = fwhm / (2. * std::sqrt(2. * std::log(2.)));
	double n = 1. / (std::sqrt(2. * M_PI) * s);
	double dxs = (x - x0) / s;
	return n * std::exp(-0.5 * dxs * dxs);
}

// argc ist Die Anzahl der Kommandozeilenparameter
// argv[i] enthält den i-ten Kommandozeilenparameter argv[0] ist Programmname
int main(int argc, char *argv[])
{
	//////////////////////////////////////////////////////////////////////////
	// Übernahme der Kommandozeilenargumente

	if (argc != 7) {
		cout << "Falscher Programmaufruf von SCIA_RETRIEVAL_2D\n";
		cout << "Aufruf: SCIA_RETRIEVAL_2D Orbitlistenpfad "
			 << "Pfad_temporäres_Arbeitsverzeichnis "
			 << "Pfad_SCIA_Sonnenspektrum "
			 << "Pfad_Sonnenrefenzspektrum "
			 << "Pfad_multips2pdf Pfad_multips2ps\n";
		cout << "Bsp: SCIA_RETRIEVAL_2D /tmp/orbit.list /tmp "
			 << "/home/meso/SCIA-DATA/SOLAR "
			 << "/home/meso/SCIA-DATA/sao_solar_ref.dat "
			 << "/home/martin/Skripts/multips2pdf "
			 << "/home/martin/Skripts/multips2ps\n";
		cout << "Programm wird abgebrochen\n";
		return 1;
	}


	string Orbitlistenpfad = argv[1];     // um Konf zu überschreiben
	string Arbeitsverzeichnis = argv[2]; // für Ausgaben in Datei(5 mal pro Spezies)
	string Solarpfad = argv[3];            // um Konf zu überschreiben
	string sol_refname = argv[4];    // solar ref for NO emission calculation
	string Pfad_multips2pdf = argv[5];
	string Pfad_multips2ps = argv[6];


	cerr << "Orbitlistenpfad: " << Orbitlistenpfad << "\n";
	cerr << "Arbeitsverzeichnis: " << Arbeitsverzeichnis << "\n";
	cerr << "Solarpfad: " << Solarpfad << "\n";
	cerr << "Solarreferenzpfad: " << sol_refname << "\n";
	cerr << "Pfad_multips2pdf: " << Pfad_multips2pdf << "\n";
	cerr << "Pfad_multips2ps: " << Pfad_multips2ps << "\n";

	//vector<string> Zusaetzliche_Orbits;
	//for (int i=0;i<Anzahl_zusaetzlicher_Orbits;i++) {
	//	Zusaetzliche_Orbits.push_back(argv[i+7]);
	//}

	//TODO ACHTUNG ACHUTNG
	// Falls ja läuft das Laden der Limbdaten anders ab und
	// die Ausschlusskriterien müssen anders angewendet werden;
	string untersuche_limb_mesothermo_states = "ja";

	string mache_Fit_Plots_limb = "nein";
	string mache_Fit_Plots_nadir = "nein";
//   string mache_Fit_Plots_MgI="ja";          // switches für einzelne Spezies
//    Plotterei als Spezieseigenschaft, implementieren volles Retrieval auch
//    string mache_Fit_Plots_MgII="nein";
//    string mache_Fit_Plots_unknown="nein";
	string mache_volles_Retrieval = "ja"; // Falls nicht, nach Teil 2 abbrechen
	string mache_volles_Retrieval_MgI = "nein"; // switches für einzelne Spezies
	string mache_volles_Retrieval_MgII = "nein";
	string mache_volles_Retrieval_unknown = "nein";
	string mache_volles_Retrieval_FeI = "nein";
	string mache_volles_Retrieval_NO = "ja";


	// Zu Arbeitsverzeichnis /////// für jede Spezies sollen 5 Dateien entstehen
	//sssss_orbit_xxxxx_yyyymmdd_hhmm_Dichten.txt
	//sssss_orbit_xxxxx_yyyymmdd_hhmm_Sx.txt
	//sssss_orbit_xxxxx_yyyymmdd_hhmm_AKM.txt
	//sssss_orbit_xxxxx_yyyy_mm_dd_hh_mm_nadir_Saeulen.txt
	//sssss_orbit_xxxxx_yyyy_mm_dd_hh_mm_0limb_Saeulen.txt

	// sssss -Spezies mit führenden 0en z.b. 00MgI
	// xxxxx -Orbit_ID  kriegt man aus Orbitlistenpfad ......../xxxxx.temp.list
	// ...also von hinten zählen und dann substr
	//yyyyymmdd_hhmm kriegt man aus dem Namen der ersten Datei in der Orbitliste
	////////////////////////////////////////////////////////////////////////////

	time_t start_zeit, timer1, deltaT;
	time_t Teil1_Start, Teil1_End, T1_Dauer;
	time_t Teil2_Start, Teil2_End, T2_Dauer;
	time_t Teil3_Start, Teil3_End, T3_Dauer;
	time_t Teil4_Start, Teil4_End, T4_Dauer;
	time_t Teil5_Start, Teil5_End, T5_Dauer;
	time_t Teil6_Start, Teil6_End, T6_Dauer;

	time(&start_zeit);
	Nachricht_Schreiben("Starte Hauptprogramm...", 3, Prioritylevel);
	/***************************************************************************
	 TEIL 1 VORBEREITUNG

	 IN DIESEM PROGRAMMTEIL WERDEN DIE DATEN DER ZU UNTERSUCHENEN SPEZIES ANGEGEBEN.
	 DIE ORBITLISTE WIRD EINGELESEN UND DAS SONNENSPEKTRUM DES MESSTAGES WIRD EINGELESEN.

	 **************************************************************************/
	cerr << "Teil1\n";
	time(&Teil1_Start);
	double wl; // central peak position wavelength
	unsigned int i, j, k, l; // Meine Zählvariablen
	i = 0;
	j = 0;
	k = 0;
	l = 0;
	Konfiguration Konf;
	Konf.Konfiguration_einlesen();
	Konf.Konfiguration_anzeigen();  //-> Test erfolgreich
	//Konf mit argv Argumenten
	Konf.m_Pfad_Datei_mit_Dateinamen_fuer_Messungen_eines_Orbits = Orbitlistenpfad;
	Konf.m_Pfad_Solar_Spektrum = Solarpfad;

	////////////////////////////////////////////////////////////////////////////
	//
	// Sonnenspektrum bestimmen
	//
	////////////////////////////////////////////////////////////////////////////
	Nachricht_Schreiben("Lade Sonnenspektrum...", 3, Prioritylevel);
	Sonnenspektrum Solspec;
	if (Solspec.Laden_SCIA(Konf.m_Pfad_Solar_Spektrum, Konf.m_Pfad_Solar_Fallback_Spektrum) != 0) {
		cout << "Programmabbruch\n";
		return -1;
	}
	//Überprüfen, ob einlesen erfolgreich war
	//Solspec.Speichern_was_geladen_wurde("CHECKDATA/Sonne_so_wie_geladen.txt");
	//ok -> funktioniert
	Nachricht_Schreiben("Lade Referenzsonnenspektrum...", 3, Prioritylevel);
	Sonnenspektrum sol_ref;
	if (sol_ref.Laden_SCIA(sol_refname, Konf.m_Pfad_Solar_Fallback_Spektrum) != 0) {
		cout << "Programmabbruch\n";
		return -1;
	}
	////////////////////////////////////////////////////////////////////////////
	//
	// Sonnenspektrum ist bestimmt
	//
	////////////////////////////////////////////////////////////////////////////

	// Fitparameter für alle Spezies vorbeiten
	Speziesfenster Spez;
	Liniendaten Lindat;
	vector<Speziesfenster> Spezies_Fenster;
	// Linienparameter für alle Spezies einlesen *******************************
	// Dies könnte man später aus einer Datei auslesen /////////////////////////
	// Die Wahl der Fensterdaten sollte Anhand von ausgewählten Plots von I/F
	// um die Linie herum getroffen werden. Das in diesem Fall sehr viel
	// einfacher, als eine Automatisierung
	//
	// TODO Falls das nochmehr spezies werden, dann Vektor von Vektor von
	// Speziesfenster Hier werden die Wellenlängen in nm angegeben.(das nm hebt
	// sich nacher mit der Multiplikation mit der Intervallbreite weg)
	//
	Nachricht_Schreiben("Bestimme Speziesdaten...", 3, Prioritylevel);
	//Füllen z.B. für Magnesium+ *******************************************
	Spez.m_Spezies_Name = "MgII"; //bitte zusammenhängend
	Spez.m_FWHM = Konf.m_FWHM;
	//Wellenlaenge 1
	wl = 279.553;
	Spez.m_Wellenlaengen.push_back(wl);
	Spez.m_Basisfenster_links_WLmin.push_back(wl - 3);
	Spez.m_Basisfenster_links_WLmax.push_back(wl - 1);
	Spez.m_Basisfenster_rechts_WLmin.push_back(wl + 1);
	Spez.m_Basisfenster_rechts_WLmax.push_back(wl + 3);
	Spez.m_Peakfenster_WLmin.push_back(wl - 3);
	Spez.m_Peakfenster_WLmax.push_back(wl + 3);
	Lindat.Einlesen(Konf.m_Pfad_Linienparameter_Metalle, wl);
	Spez.m_Liniendaten.push_back(Lindat);
	/*
	//Wellenlänge 2
	wl = 280.270;
	Spez.m_Wellenlaengen.push_back(wl);
	Spez.m_Basisfenster_links_WLmin.push_back(wl - 3);
	Spez.m_Basisfenster_links_WLmax.push_back(wl - 1);
	Spez.m_Basisfenster_rechts_WLmin.push_back(wl + 1);
	Spez.m_Basisfenster_rechts_WLmax.push_back(wl + 3);
	Spez.m_Peakfenster_WLmin.push_back(wl - 3);
	Spez.m_Peakfenster_WLmax.push_back(wl + 3);
	Lindat.Einlesen(Konf.m_Pfad_Linienparameter_Metalle, wl);
	Spez.m_Liniendaten.push_back(Lindat);
	// */
	// weitere Wellenlängen
	Spezies_Fenster.push_back(Spez);
	// SpezVektoren wieder leeren
	// hier könnte man in Spez auch eine Methode für schreiben
	Spez.m_Wellenlaengen.resize(0);
	Spez.m_Basisfenster_links_WLmin.resize(0);
	Spez.m_Basisfenster_links_WLmax.resize(0);
	Spez.m_Basisfenster_rechts_WLmin.resize(0);
	Spez.m_Basisfenster_rechts_WLmax.resize(0);
	Spez.m_Peakfenster_WLmin.resize(0);
	Spez.m_Peakfenster_WLmax.resize(0);
	Spez.m_Liniendaten.resize(0);
	//weitere Spezies
	// Magnesium ************************************************************
	Spez.m_Spezies_Name = "MgI"; //bitte zusammenhängend
	Spez.m_FWHM = Konf.m_FWHM;
	//Wellenlaenge 1
	wl = 285.213;
	Spez.m_Wellenlaengen.push_back(wl);
	Spez.m_Basisfenster_links_WLmin.push_back(wl - 3);
	//Spez.m_Basisfenster_links_WLmin.push_back(Spez.m_Wellenlaengen[0]-10);
	Spez.m_Basisfenster_links_WLmax.push_back(wl - 1);
	Spez.m_Basisfenster_rechts_WLmin.push_back(wl + 1);
	Spez.m_Basisfenster_rechts_WLmax.push_back(wl + 3);
	//Spez.m_Basisfenster_rechts_WLmax.push_back(Spez.m_Wellenlaengen[0]+4);
	Spez.m_Peakfenster_WLmin.push_back(wl - 3);
	Spez.m_Peakfenster_WLmax.push_back(wl + 3);
	Lindat.Einlesen(Konf.m_Pfad_Linienparameter_Metalle, wl);
	Spez.m_Liniendaten.push_back(Lindat);
	// weitere Wellenlängen
	//**************************************************************************
	Spezies_Fenster.push_back(Spez);
	//weitere Spezies
	// SpezVektoren wieder leeren
	// hier könnte man in Spez auch eine Methode für schreiben
	Spez.m_Wellenlaengen.resize(0);
	Spez.m_Basisfenster_links_WLmin.resize(0);
	Spez.m_Basisfenster_links_WLmax.resize(0);
	Spez.m_Basisfenster_rechts_WLmin.resize(0);
	Spez.m_Basisfenster_rechts_WLmax.resize(0);
	Spez.m_Peakfenster_WLmin.resize(0);
	Spez.m_Peakfenster_WLmax.resize(0);
	Spez.m_Liniendaten.resize(0);
	Spez.clear();  // Instanz leeren
	// Unbekannte Spezies ******************************************************
	Spez.m_Spezies_Name = "unknown"; //bitte zusammenhängend
	Spez.m_FWHM = Konf.m_FWHM;
	//Wellenlaenge 1
	wl = 288.2;
	Spez.m_Wellenlaengen.push_back(wl);
	Spez.m_Basisfenster_links_WLmin.push_back(wl - 2);
	//Spez.m_Basisfenster_links_WLmin.push_back(Spez.m_Wellenlaengen[0]-10);
	Spez.m_Basisfenster_links_WLmax.push_back(wl - 1);
	Spez.m_Basisfenster_rechts_WLmin.push_back(wl + 1);
	Spez.m_Basisfenster_rechts_WLmax.push_back(wl + 3);
	//Spez.m_Basisfenster_rechts_WLmax.push_back(Spez.m_Wellenlaengen[0]+4);
	Spez.m_Peakfenster_WLmin.push_back(wl - 2);
	Spez.m_Peakfenster_WLmax.push_back(wl + 3);
	// Lindat sind unbekannt
	Lindat.m_Wellenlaenge = wl;
	Lindat.m_rel_Einstein = 1; // stimmt nicht
	Lindat.m_f_Wert = 0.162; //0.162;  // für Si
	Lindat.m_E1 = 0; //0.01;         //für Si
	Lindat.m_E2 = 1; //0.99;         //für Si
	// Dies sind beliebige Werte, die auch nicht unbedingt auf I/F führen(prop
	// faktoren wie pi r_elektron)...sondern auf beliebige werte
	Spez.m_Liniendaten.push_back(Lindat);
	// weitere Wellenlängen
	//**************************************************************************
	Spezies_Fenster.push_back(Spez);
	//cout<<Spezies_Fenster.size();
	// cout<<Spezies_Fenster[1].m_Liniendaten[0].m_E1<<"\t";
	//weitere Spezies
	// SpezVektoren wieder leeren
	// hier könnte man in Spez auch eine Methode für schreiben
	Spez.m_Wellenlaengen.resize(0);
	Spez.m_Basisfenster_links_WLmin.resize(0);
	Spez.m_Basisfenster_links_WLmax.resize(0);
	Spez.m_Basisfenster_rechts_WLmin.resize(0);
	Spez.m_Basisfenster_rechts_WLmax.resize(0);
	Spez.m_Peakfenster_WLmin.resize(0);
	Spez.m_Peakfenster_WLmax.resize(0);
	Spez.m_Liniendaten.resize(0);
	Spez.clear();  // Instanz leeren
	// FeI ************************************************************
	Spez.m_Spezies_Name = "FeI"; //bitte zusammenhängend
	Spez.m_FWHM = Konf.m_FWHM;
	//Wellenlaenge 1
	wl = 248.32707;
	Spez.m_Wellenlaengen.push_back(wl);
	Spez.m_Basisfenster_links_WLmin.push_back(wl - 1);
	Spez.m_Basisfenster_links_WLmax.push_back(wl - 0.5);
	Spez.m_Basisfenster_rechts_WLmin.push_back(wl + 1);
	Spez.m_Basisfenster_rechts_WLmax.push_back(wl + 2);
	Spez.m_Peakfenster_WLmin.push_back(wl - 1);
	Spez.m_Peakfenster_WLmax.push_back(wl + 2);
	// Liniendaten
	Lindat.m_Wellenlaenge = wl;
	Lindat.m_rel_Einstein = 0.19; // stimmt nicht
	Lindat.m_f_Wert = 0.543;
	// ^5D_4 nach ^5F*_5 übergang j=4, dj=1 (Funktion für Matlab geschrieben)
	Lindat.m_E1 = 0.1733;
	Lindat.m_E2 = 0.8267;
	// Dies sind beliebige Werte, die auch nicht unbedingt auf I/F führen(prop
	// faktoren wie pi r_elektron)...sondern auf beliebige werte
	Spez.m_Liniendaten.push_back(Lindat);
	// weitere Wellenlängen
	//**************************************************************************
	Spezies_Fenster.push_back(Spez);
	// SpezVektoren wieder leeren
	// hier könnte man in Spez auch eine Methode für schreiben
	Spez.m_Wellenlaengen.resize(0);
	Spez.m_Basisfenster_links_WLmin.resize(0);
	Spez.m_Basisfenster_links_WLmax.resize(0);
	Spez.m_Basisfenster_rechts_WLmin.resize(0);
	Spez.m_Basisfenster_rechts_WLmax.resize(0);
	Spez.m_Peakfenster_WLmin.resize(0);
	Spez.m_Peakfenster_WLmax.resize(0);
	Spez.m_Liniendaten.resize(0);
	Spez.clear();  // Instanz leeren
	//weitere Spezies

	// NO stuff
	// from config file
	for (i = 0; i < Konf.no_NO_transitions; i++) {
		Spez.m_Spezies_Name = "NO";
		NO_emiss NO(Konf.NO_v_u.at(i), Konf.NO_v_l.at(i),
				Konf.NO_v_l_abs.at(i), Konf.atmo_Temp);
		NO.get_solar_data(sol_ref);
		NO.read_luque_data_from_file("DATA/Luqueetal.dat");
		NO.calc_excitation();
		NO.calc_line_emissivities();
		std::cout << "NO transition " << i
			<< ": v_u = " << NO.get_vu()
			<< ", v_l = " << NO.get_vl()
			<< ", v_l_abs = " << NO.get_vl_abs()
			<< " at (initial) " << Konf.atmo_Temp << " K" << std::endl;
		//
		Spez.NO_vec.push_back(NO);
		wl = 246.9; // dummy, will be set later more accurately
		Spez.m_Wellenlaengen.push_back(wl);
		Spez.m_Basisfenster_links_WLmin.push_back(wl - 1);
		Spez.m_Basisfenster_links_WLmax.push_back(wl - 0.5);
		Spez.m_Basisfenster_rechts_WLmin.push_back(wl + 0.5);
		Spez.m_Basisfenster_rechts_WLmax.push_back(wl + 1);
		Spez.m_Peakfenster_WLmin.push_back(wl - 1);
		Spez.m_Peakfenster_WLmax.push_back(wl + 1);
		Spez.m_FWHM = Konf.m_FWHM;
		// Liniendaten
		Lindat.m_Wellenlaenge = wl;
		Lindat.m_rel_Einstein = 1;
		Lindat.m_f_Wert = 0.162;
		Lindat.m_E1 = 0;
		Lindat.m_E2 = 1;
		Spez.m_Liniendaten.push_back(Lindat);
	}
	Spezies_Fenster.push_back(Spez);

	// SpezVektoren wieder leeren
	// hier könnte man in Spez auch eine Methode für schreiben
	Spez.m_Wellenlaengen.resize(0);
	Spez.m_Basisfenster_links_WLmin.resize(0);
	Spez.m_Basisfenster_links_WLmax.resize(0);
	Spez.m_Basisfenster_rechts_WLmin.resize(0);
	Spez.m_Basisfenster_rechts_WLmax.resize(0);
	Spez.m_Peakfenster_WLmin.resize(0);
	Spez.m_Peakfenster_WLmax.resize(0);
	Spez.m_Liniendaten.resize(0);
	Spez.clear();  // Instanz leeren
	//**************************************************************************

	////////////////////////////////////////////////////////////////////////////
	// Linienparameter für alle Spezies einlesen  ende *************************
	Nachricht_Schreiben("Speziesdaten bestimmt...", 3, Prioritylevel);
	// Für jede Spezies einen Vector mit Messergebnissen erstellen
	// Mg und Mg+
	vector<Ausgewertete_Messung_Limb> Ausgewertete_Limbmessung_MgI;
	vector<Ausgewertete_Messung_Limb> Ausgewertete_Limbmessung_MgII;
	vector<Ausgewertete_Messung_Nadir> Ausgewertete_Nadirmessung_MgI;
	vector<Ausgewertete_Messung_Nadir> Ausgewertete_Nadirmessung_MgII;

	vector<Ausgewertete_Messung_Limb> Ausgewertete_Limbmessung_unknown;
	vector<Ausgewertete_Messung_Nadir> Ausgewertete_Nadirmessung_unknown;
	vector<Ausgewertete_Messung_Limb> Ausgewertete_Limbmessung_FeI;
	vector<Ausgewertete_Messung_Nadir> Ausgewertete_Nadirmessung_FeI;
	vector<Ausgewertete_Messung_Limb> Ausgewertete_Limbmessung_NO;
	vector<Ausgewertete_Messung_Nadir> Ausgewertete_Nadirmessung_NO;
	////////////////////////////////////////////////////////////////////////////
	// Orbitliste Laden
	//
	// Die Orbitliste ist so zu erstellen, dass die erste Datei eine Limb Datei
	// ist Da dieses Verfahren eh nur sinnvoll ist, wenn Limbmessungen
	// vorhanden sind, sollte dies kein Problem sein
	// TODO Die Orbitliste sollte nach Limb und Nadir sortiert werden. Danach
	// muss sie nach Zeit sortiert sein( Wichtig für Einteilung der Boxen
	// etc.... falls Ordner nach Dateinamen sortiert, so sollte das immer
	// stimmen) Die Frage ist, ob man die Orbitliste so übernimmt, oder ob man
	// hier nochmal sortiert. Bis jetzt ist die Orbitliste schon automatisch
	// so. Allerdings wär eine Sortierung der Orbitliste ein geringer
	// Zeitaufwand
	//
	////////////////////////////////////////////////////////////////////////////
	Nachricht_Schreiben("Lade Orbitliste...", 3, Prioritylevel);
	Orbitliste Orbitlist;
	if (Orbitlist.Liste_Laden(Konf.m_Pfad_Datei_mit_Dateinamen_fuer_Messungen_eines_Orbits) != 0) {
		cout << "Programmabbruch\n";
		return -1;
	}
	// Zur Überprüfung Orbitliste in Datei schreiben
	//Orbitlist.In_Datei_Speichern("CHECKDATA/Orbitliste_Ueberpruefung.txt");
	// der ordner CHECKDATA muss vorher existieren
	//cout<<Orbitlist.m_Dateinamen[Orbitlist.m_Dateinamen.size()-1];//
	//TODO Das ist leer...Es muss beim erstellen der Orbitliste
	// sichergestellt werden, dass die letzte Zeile kein\n enthält!!!!
	// -> Funktion funktioniert... aber Orbitliste Datei muss  ordentlich
	// erzeugt sein
	////////////////////////////////////////////////////////////////////////////
	//
	// Orbitliste Ist geladen
	//
	////////////////////////////////////////////////////////////////////////////

	// convolve high resolution solar spectra with the sciamachy resolution function
	// TODO: replace the detection by a more sophisticated method.
	if (Solspec.m_Anzahl_WL > 4000)
		Solspec.saoref_to_sciamachy();

	time(&Teil1_End);
	T1_Dauer = Teil1_End - Teil1_Start;
	/***************************************************************************
	 ENDE TEIL 1
	 **************************************************************************/
	/***************************************************************************
	TEIL 2 BESTIMMUNG DER SÄULENDICHTEN
	FÜR ALLE LIMB- UND NADIRMESSUNGEN EINES ORBITS WERDEN DIE SÄULENDICHTEN FÜR
	JEDE LINIE BESTIMMT UND IN EINE EINFACHERE STRUKTUR(Ausgewertete_Messungen)
	GESTECKT.
	***************************************************************************/
	cerr << "Teil2\n";
	time(&Teil2_Start);
	////////////////////////////////////////////////////////////////////////////
	//
	// Säulendichtenbestimmung
	//
	////////////////////////////////////////////////////////////////////////////
	Nachricht_Schreiben("Bestimme Säulendichten...", 3, Prioritylevel);
	// Aus jeder Messung egal ob Limb oder Nadir die Säulendichte für jede
	// Linie jeder Spezies bestimmen
	//counter für Qualitätscheck der Messung initialisieren
	int counter_Nachtmessungen = 0;
	//TODO ausgabe der counter später implementieren
	int counter_NLC_detektiert = 0;
	int counter_Richtungsvektor_nicht_ok = 0;
	int counter_Nachtmessungen_Nadir = 0;
	int counter_Nadir_Nacht_Dateien = 0;
	// Aus jeder Messung egal ob Limb oder Nadir die Säulendichte für jede
	// Linie jeder Spezies bestimmen
	for (l = 0; l < Orbitlist.m_Dateinamen.size(); l++) {
		string L1CDatei = Orbitlist.m_Dateinamen[l];
		//cout<<L1CDatei<<"\n";

		//Falls Limb-> Limbauslesung
		if (Orbitlist.Ist_Messung_Limbmessung(l) == true) {
			/////////////////////////////////////
			//Limbauswertung
			////////////////////////////////////
			// Die erste Limbdatei ist wichtig für die Interpolation des
			// Sonnenspektrums... es ist also wichtig, dass für l=0 eine
			// Limbmessung vorliegt
			//cerr<<"limbauswertung start\n";
			//NotTODO Die Säulendichtebestimmung kann deutlich schneller
			//geschehen, Verbesserungen hier sind aber irrelevant
			Limb_Auswertung(Orbitlist, l, Solspec, Spezies_Fenster,
							counter_Nachtmessungen, counter_NLC_detektiert,
							counter_Richtungsvektor_nicht_ok,
							Arbeitsverzeichnis, mache_Fit_Plots_limb,
							untersuche_limb_mesothermo_states,
							Ausgewertete_Limbmessung_MgI,
							Ausgewertete_Limbmessung_MgII,
							Ausgewertete_Limbmessung_unknown,
							Ausgewertete_Limbmessung_FeI,
							Ausgewertete_Limbmessung_NO,
							Konf);
			//cerr<<"limbauswertung Ende\n";
			// Die Zwischenergebnisse stehen nun in
			// Ausgewertete_Limbmessung_MgI und  Ausgewertete_Limbmessung_MgII
			/////////////////////////////////////
			//Limbauswertung Ende
			////////////////////////////////////
		}//ende if (Orbitlist.Ist_Messung_Limbmessung(i)==0)
		//Falls Nadir-> Nadirauslesung
		//TODO Ladeunterroutine überprüfen
		if (Orbitlist.Ist_Messung_Nadirmessung(l) == true) {
			/////////////////////////////////////
			//Nadirauswertung
			/////////////////////////////////////
			// TODO Die Zeilendichtebestimmung kann deutlich schneller
			// geschehen...da werden viel zu viele unnötig große Schleifen
			// verwendet
			//cerr<<"Nadirauswertung start\n";
			Nadir_Auswertung(Orbitlist, l, Solspec, Spezies_Fenster,
							 counter_Nachtmessungen_Nadir,
							 counter_Nadir_Nacht_Dateien,
							 Arbeitsverzeichnis, mache_Fit_Plots_nadir,
							 Ausgewertete_Nadirmessung_MgI,
							 Ausgewertete_Nadirmessung_MgII,
							 Ausgewertete_Nadirmessung_unknown,
							 Ausgewertete_Nadirmessung_FeI,
							 Ausgewertete_Nadirmessung_NO);
			/////////////////////////////////////
			//Nadirauswertung Ende
			////////////////////////////////////
		}//ENDE if(Orbitlist.Ist_Messung_Nadirmessung(l)==0)
	}//ende Schleife über l

	//////// TEIL 2B/////////////////////
	// Glätten der Ausgewerteten Messunge
	/*
	for (int i = 0; i < 3; i++) {
		SCD_Glaettung(Ausgewertete_Limbmessung_MgI, 1, untersuche_limb_mesothermo_states);
		SCD_Glaettung(Ausgewertete_Limbmessung_unknown, 1, untersuche_limb_mesothermo_states);
		SCD_Glaettung(Ausgewertete_Limbmessung_MgII, 2, untersuche_limb_mesothermo_states);
		SCD_Glaettung(Ausgewertete_Limbmessung_FeI, 1, untersuche_limb_mesothermo_states);
	}
	// */

	//////// TEIL 2B/////////////////////
	////////////////////////////////////////////////////////////////////////////
	//
	// Zeilendichtenbestimmung ist abgeschlossen
	// Bis hierhin braucht das Programm keine Sekunde, wenn die ROHDATEN auf
	// der LOKALEN Platte liegen.
	//
	////////////////////////////////////////////////////////////////////////////
	Nachricht_Schreiben("Zeilendichten sind bestimmt....", 3, Prioritylevel);
	////////////////////////////////////////////////////////////////////////////
	//Ausgabe der Vektoren mit den Zwischenergebnissen in Dateien
	// TODO sieht Ok aus...evtl später mal ein par fits angucken
	string xxxxx = xxxxx_Bestimmen(Orbitlistenpfad);
	string yyyymmdd_hhmm = yyyymmdd_hhmm_Bestimmen(Orbitlist.m_Dateinamen[0]);
	string sssss_MgI = "00MgI";
	string sssss_MgII = "0MgII";
	string sssss_unknown = "unkno";
	string sssss_FeI = "00FeI";
	string sssss_NO = "000NO";
	string Endung_Limb = "_0limb_Saeulen.txt";
	string Endung_Limb_back = "_1limb_Saeulen.txt";
	string Endung_Nadir = "_nadir_Saeulen.txt";
	string Dateiout_Mittelteil = "_orbit_" + xxxxx + "_" + yyyymmdd_hhmm;
	//sssss_orbit_xxxxx_yyyy_mm_dd_hh_mm_0limb_Saeulen.txt
	string Pfad_Saeulen_Limb_MgI = Arbeitsverzeichnis + "/" + sssss_MgI + Dateiout_Mittelteil + Endung_Limb;
	string Pfad_Saeulen_Limb_MgII = Arbeitsverzeichnis + "/" + sssss_MgII + Dateiout_Mittelteil + Endung_Limb;
	string Pfad_Saeulen_Limb_unknown = Arbeitsverzeichnis + "/" + sssss_unknown + Dateiout_Mittelteil + Endung_Limb;
	string Pfad_Saeulen_Limb_FeI = Arbeitsverzeichnis + "/" + sssss_FeI + Dateiout_Mittelteil + Endung_Limb;
	string Pfad_Saeulen_Limb_NO = Arbeitsverzeichnis + "/" + sssss_NO + Dateiout_Mittelteil + Endung_Limb;
	string Pfad_Saeulen_Limb_NO_back = Arbeitsverzeichnis + "/" + sssss_NO
		+ Dateiout_Mittelteil + Endung_Limb_back;
	string Pfad_Saeulen_Nadir_MgI = Arbeitsverzeichnis + "/" + sssss_MgI + Dateiout_Mittelteil + Endung_Nadir;
	string Pfad_Saeulen_Nadir_MgII = Arbeitsverzeichnis + "/" + sssss_MgII + Dateiout_Mittelteil + Endung_Nadir;
	string Pfad_Saeulen_Nadir_unknown = Arbeitsverzeichnis + "/" + sssss_unknown + Dateiout_Mittelteil + Endung_Nadir;
	string Pfad_Saeulen_Nadir_FeI = Arbeitsverzeichnis + "/" + sssss_FeI + Dateiout_Mittelteil + Endung_Nadir;
	string Pfad_Saeulen_Nadir_NO = Arbeitsverzeichnis + "/" + sssss_NO + Dateiout_Mittelteil + Endung_Nadir;
	//cerr<<"Ausgabe_Saeulendichten\n";
	Ausgabe_Saeulendichten(Pfad_Saeulen_Limb_MgI.c_str(), Ausgewertete_Limbmessung_MgI);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Limb_MgII.c_str(), Ausgewertete_Limbmessung_MgII);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Limb_unknown.c_str(), Ausgewertete_Limbmessung_unknown);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Limb_FeI.c_str(), Ausgewertete_Limbmessung_FeI);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Limb_NO.c_str(), Ausgewertete_Limbmessung_NO);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Nadir_MgI.c_str(), Ausgewertete_Nadirmessung_MgI);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Nadir_MgII.c_str(), Ausgewertete_Nadirmessung_MgII);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Nadir_unknown.c_str(), Ausgewertete_Nadirmessung_unknown);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Nadir_FeI.c_str(), Ausgewertete_Nadirmessung_FeI);
	Ausgabe_Saeulendichten(Pfad_Saeulen_Nadir_NO.c_str(), Ausgewertete_Nadirmessung_NO);
	//cerr<<"Plotdateien zusammenfassen\n";
	//Plotdateien zusammenfassen
	if ((mache_Fit_Plots_limb == "ja") || (mache_Fit_Plots_nadir == "ja")) {
		vector<Speziesfenster>::iterator sfit;
		for (sfit = Spezies_Fenster.begin(); sfit != Spezies_Fenster.end(); ++sfit) {
			// Namen der Ausgabedatei zusammenschustern
			string pdf_Datei = Arbeitsverzeichnis + "/Plots/Orbit_" + xxxxx
				+ "_" + sfit->m_Spezies_Name + "Fits.pdf";
			Plots_Zusammenfassen(Pfad_multips2pdf, Pfad_multips2ps, pdf_Datei,
					sfit->m_Liste_der_Plot_Dateinamen);
		}
	}
	////////////////////////////////////////////////////////////////////////////
	time(&Teil2_End);
	T2_Dauer = Teil2_End - Teil2_Start;
	if (mache_volles_Retrieval != "ja") { // Falls nicht, nach Teil 2 abbrechen
		// Hier könnt man echtmal ne sprunganweisung machen...ich verlass mich
		// erstmal auf den garbage collector von c++
		cerr << "Programm wird vorzeitig nach Säulendichtenbestimmung beendet\n";
		return 0;
	}

	/***************************************************************************
	ENDE TEIL 2

	***************************************************************************/

	/***************************************************************************
	TEIL 3 AUFBAU DER FÜR DAS RETRIEVAL BENÖTIGTEN MATRIZZEN
	MEHRERE EINFACHE MATRIZZEN WERDEN BENÖTIGT.
	DIE WICHTIGSTE UND KOMPLIZIERTESTE MATRIX IST DIE MATRIX FÜR DIE
	LUFTMASSENFAKTOREN AMF.
	***************************************************************************/
	cerr << "Teil3\n";
	time(&Teil3_Start);
	////////////////////////////////////////////////////////////////////////////
	// Matrizen und Vektoren aufbauen
	////////////////////////////////////////////////////////////////////////////
	/***************************************************************************
	// Wir benötigen das Gitter, auf dem wir die Dichte ausrechnen
	// Ab nun Vektoren und Matrizen für jede Spezies bestimmen
	// einen Vektor der Dichte_n
	// einen Vektor x_a für die a-priori-Lösung der Dichte (oder Startwert usw)
	// einen Vektor der Zeilendichten
	// Die Gesamtwichtungsfaktoren Lambda_H und Lambda_PHI
	// zwei Matrizen S_H und S_LAT mit den Unrterschieden der
	// Nachbarpunkte(zur Glättung) die Matrix K, in der das Absorptionsgesetz
	// für die Lichtwege drinsteckt
	// Kovarianzmatrizen S_y und S_a...die man aber als Diagonalvektor benutzt
	***************************************************************************/
	//das hier weg( das erlaubt nur Limbmessungen)
	Ausgewertete_Nadirmessung_MgI.resize(0);
	//das hier weg( das erlaubt nur Limbmessungen)
	Ausgewertete_Nadirmessung_MgII.resize(0);
	Ausgewertete_Nadirmessung_unknown.resize(0);
	Ausgewertete_Nadirmessung_FeI.resize(0);
	Retrievalgitter Grid;
	//5.0 ist vernünftig; 2.0 ist gut für TAU_LOS_plot
	double Mindestabstand_Lat_in_Grad = 4.0;
	Grid.Retrievalgitter_erzeugen(Ausgewertete_Limbmessung_MgI,
			Mindestabstand_Lat_in_Grad, Konf);
	//Folgende Ausgabe sieht ok aus 28.9.2010 (Durchstoßpunkte werden erst
	//später ermittelt)
	string Pfad_Grid = Arbeitsverzeichnis + "/" + "Gitter.txt";
	Grid.In_Datei_Ausgeben(Pfad_Grid);

	// ein Vektor der totalen Anzahldichte
	MPL_Matrix Dichte_n_tot(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_n_tot.Null_Initialisierung();
	prepare_total_density(Grid, Dichte_n_tot, Ausgewertete_Limbmessung_NO);
	////////////////////////////////////////////////////////////////////////////
	// Spezies Mg I //
	////////////////////////////////////////////////////////////////////////////
	// einen Vektor der Dichte_n
	MPL_Matrix Dichte_n_MgI(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_n_MgI.Null_Initialisierung();
	//Dichte_n_MgI.in_Datei_speichern("/tmp/mlangowski/0/Dichte_n_MgI");
	// einen Vektor x_a für die a-priori-Lösung der Dichte (oder Startwert usw)
	MPL_Matrix Dichte_apriori_MgI(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_apriori_MgI.Null_Initialisierung();  // Null als Startwert
	// Vektor mit Zeilendichten  für alle Messungen einer Spezies
	MPL_Matrix Saeulendichten_MgI(Ausgewertete_Limbmessung_MgI.size()
			+ Ausgewertete_Nadirmessung_MgI.size(), 1); //Spaltenvektor
	MPL_Matrix Saeulendichten_Fehler_MgI(Ausgewertete_Limbmessung_MgI.size()
			+ Ausgewertete_Nadirmessung_MgI.size(), 1); //Spaltenvektor
	// verwendete Konfiguration bis 20.1.2011
	//double MgI_Lambda_Hoehe= 5E-7;//5E-6;//5E-6;//5E-6;     gute Werte 5 E-6 
	// TODO das laden der Parameter eindeutiger machen
	// Konf.m_Retrieval_Kovarianzen[0+1*m_Anzahl_der_Emitter];
	//double MgI_Lambda_Breite= 5E-6;//1E-6;//1E-7;//1E-7;      gute Werte 1 E-6
	//double MgI_Lambda_apriori= 5E-5;//5E-6 oder -5;
	//double MgI_Lambda_Hoehe= 5E-7;      // Gut für 3 km Schirtte mit 30 Höhen
	//double MgI_Lambda_Breite= 5E-6;     // Gut für 3 km Schirtte mit 30 Höhen
	//double MgI_Lambda_apriori=5E-5;     // Gut für 3 km Schirtte mit 30 Höhen

	//double MgI_Lambda_Hoehe= 1E-5;      // Gut für 1 km Schirtte mit 82 Höhen
	//double MgI_Lambda_Breite= 1E-5;
	//double MgI_Lambda_apriori=5E-5;
	//double MgI_Lambda_letzte_Hoehe=100;

//   double MgI_Lambda_Hoehe= 1E-5;
//   double MgI_Lambda_Breite= 1E-5;
//   double MgI_Lambda_apriori=5E-5;
//   double MgI_Lambda_letzte_Hoehe=1E-32; //100 für ein

	// Spezies Index for MgI
	int spez_index = 1;

	double MgI_Lambda_letzte_Hoehe = 1E-32; //100 für ein

	// take the lambdas from the config file, no need to re-compile
	double MgI_Lambda_apriori
		= Konf.m_Retrieval_Kovarianzen[spez_index + 0 * Konf.m_Anzahl_der_Emitter];
	double MgI_Lambda_Breite
		= Konf.m_Retrieval_Kovarianzen[spez_index + 1 * Konf.m_Anzahl_der_Emitter];
	double MgI_Lambda_Hoehe
		= Konf.m_Retrieval_Kovarianzen[spez_index + 2 * Konf.m_Anzahl_der_Emitter];

	cout << "Retrieval Kovarianzen MgI:" << endl;
	cout << "L_apriori = " << MgI_Lambda_apriori << ", "
		 << "L_Breite = " << MgI_Lambda_Breite << ", "
		 << "L_Höhe = " << MgI_Lambda_Hoehe << endl;

	// Speziesunabhängige Matrizen, alle sind quadratisch
	MPL_Matrix S_Breite(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_Hoehe(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_letzte_Hoehe_MgI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// Messfehlermatrix S_y y*y
	// quadratische matrix
	MPL_Matrix S_y_MgI(Saeulendichten_MgI.m_Elementanzahl, Saeulendichten_MgI.m_Elementanzahl);
	MPL_Matrix S_apriori_MgI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// AMF_Matrix erstellen // Zeilen wie y spalten wie x
	MPL_Matrix AMF_MgI(Saeulendichten_MgI.m_Elementanzahl, Grid.m_Anzahl_Punkte);
	// S_Breite, S_Hoehe, S_Apriori, S_y, AMF aufbauen
	int IERR = 0;
	//cerr<<"MgI Matrizen aufbaun\n";
	if (mache_volles_Retrieval_MgI == "ja") {
		Matrizen_Aufbauen(S_Breite, S_Hoehe, S_letzte_Hoehe_MgI,
						   MgI_Lambda_letzte_Hoehe,
						   S_apriori_MgI, S_y_MgI,
						   AMF_MgI, MgI_Lambda_apriori,
						   // TODO 1 ist MgI, 0 ist MgII hier wieder korrigieren
						   Spezies_Fenster[spez_index],
						   Grid,
						   Ausgewertete_Limbmessung_MgI,
						   Ausgewertete_Nadirmessung_MgI,
						   Konf, IERR);
		if (IERR != 0) {
			cerr << "Fehler bei Matrizen_Aufbauen\n";
			return -1; //Hauptprogramm beenden
		}
	}
	// Säulendichten und Fehler auffüllen (Fehler für Wichtungsmatrixberechnung)
	// Limb MgI
	//cerr<<"MgI Limb\n";
	for (unsigned int i = 0; i < Ausgewertete_Limbmessung_MgI.size(); i++) {
		Saeulendichten_MgI(i) = Ausgewertete_Limbmessung_MgI[i].m_Zeilendichte;
		Saeulendichten_Fehler_MgI(i)
			= Ausgewertete_Limbmessung_MgI[i].m_Fehler_Zeilendichten;
	}
	// Nadir MgI
	//cerr<<"MgI Nadir\n";
	for (unsigned int i = Ausgewertete_Limbmessung_MgI.size();
			i < Ausgewertete_Limbmessung_MgI.size() + Ausgewertete_Nadirmessung_MgI.size(); i++) {
		int Nadir_i = i - Ausgewertete_Limbmessung_MgI.size();
		Saeulendichten_MgI(i)
			= Ausgewertete_Nadirmessung_MgI[Nadir_i].m_Zeilendichte;
		Saeulendichten_Fehler_MgI(i)
			= Ausgewertete_Nadirmessung_MgI[Nadir_i].m_Fehler_Zeilendichten;
	}
	if (mache_volles_Retrieval_MgI == "ja")
		generate_Sy(S_y_MgI, Saeulendichten_Fehler_MgI);
	//Ende Säulendichten und Fehler auffüllen
	////////////////////////////////////////////////////////////////////////////
	// ENDE Spezies Mg I //
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// Spezies Mg II //
	////////////////////////////////////////////////////////////////////////////
	// einen Vektor der Dichte_n
	MPL_Matrix Dichte_n_MgII(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_n_MgII.Null_Initialisierung();
	//Dichte_n_MgI.in_Datei_speichern("/tmp/mlangowski/0/Dichte_n_MgI");
	// einen Vektor x_a für die a-priori-Lösung der Dichte (oder Startwert usw)
	MPL_Matrix Dichte_apriori_MgII(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_apriori_MgII.Null_Initialisierung();  // Null als Startwert
	// Vektor mit Zeilendichten  für alle Messungen einer Spezies
	//quadratische matrix
	MPL_Matrix S_letzte_Hoehe_MgII(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix Saeulendichten_MgII(Ausgewertete_Limbmessung_MgII.size()
			+ Ausgewertete_Nadirmessung_MgII.size(), 1); //Spaltenvektor
	MPL_Matrix Saeulendichten_Fehler_MgII(Ausgewertete_Limbmessung_MgII.size()
			+ Ausgewertete_Nadirmessung_MgII.size(), 1); //Spaltenvektor
	// Konfiguration bis januar 2011
	//double MgII_Lambda_Hoehe= 5E-6;//5E-6;//5E-6;//5E-6;
	//double MgII_Lambda_Breite= 1E-6;//1E-6;//1E-7;//1E-7;
	//double MgII_Lambda_apriori= 1E-3;//5E-6 oder -5;

	// Spezies Index for MgII
	spez_index = 0;

	double MgII_Lambda_letzte_Hoehe = 1E-32; //100 für ein

	// take the lambdas from the config file, no need to re-compile
	double MgII_Lambda_apriori
		= Konf.m_Retrieval_Kovarianzen[spez_index + 0 * Konf.m_Anzahl_der_Emitter];
	double MgII_Lambda_Breite
		= Konf.m_Retrieval_Kovarianzen[spez_index + 1 * Konf.m_Anzahl_der_Emitter];
	double MgII_Lambda_Hoehe
		= Konf.m_Retrieval_Kovarianzen[spez_index + 2 * Konf.m_Anzahl_der_Emitter];

	cout << "Retrieval Kovarianzen MgII:" << endl;
	cout << "L_apriori = " << MgII_Lambda_apriori << ", "
		 << "L_Breite = " << MgII_Lambda_Breite << ", "
		 << "L_Höhe = " << MgII_Lambda_Hoehe << endl;
	// Messfehlermatrix S_y y*y
	// quadratische matrix
	MPL_Matrix S_y_MgII(Saeulendichten_MgII.m_Elementanzahl, Saeulendichten_MgII.m_Elementanzahl);
	MPL_Matrix S_apriori_MgII(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// AMF_Matrix erstellen // Zeilen wie y spalten wie x
	MPL_Matrix AMF_MgII(Saeulendichten_MgII.m_Elementanzahl, Grid.m_Anzahl_Punkte);
	// S_Breite, S_Hoehe, S_Apriori, S_y, AMF aufbauen
	IERR = 0;
	//cerr<<"MgII Matrizen aufbaun\n";
	if (mache_volles_Retrieval_MgII == "ja") {
		Matrizen_Aufbauen(S_Breite, S_Hoehe, S_letzte_Hoehe_MgII,
						   MgII_Lambda_letzte_Hoehe,
						   S_apriori_MgII, S_y_MgII,
						   AMF_MgII, MgII_Lambda_apriori,
						   // TODO 1 ist MgI, 0 ist MgII hier wieder korrigieren
						   Spezies_Fenster[spez_index],
						   Grid,
						   Ausgewertete_Limbmessung_MgII,
						   Ausgewertete_Nadirmessung_MgII,
						   Konf, IERR);
		if (IERR != 0) {
			cerr << "Fehler bei Matrizen_Aufbauen\n";
			return -1; //Hauptprogramm beenden
		}
	}
	// Säulendichten und Fehler auffüllen (Fehler für Wichtungsmatrixberechnung)
	// Limb MgII
	//cerr<<"MgII Limb\n";
	for (unsigned int i = 0; i < Ausgewertete_Limbmessung_MgII.size(); i++) {
		Saeulendichten_MgII(i) = Ausgewertete_Limbmessung_MgII[i].m_Zeilendichte;
		Saeulendichten_Fehler_MgII(i)
			= Ausgewertete_Limbmessung_MgII[i].m_Fehler_Zeilendichten;
	}
	// Nadir MgII
	//cerr<<"MgII Nadir\n";
	for (unsigned int i = Ausgewertete_Limbmessung_MgII.size();
			i < Ausgewertete_Limbmessung_MgII.size() + Ausgewertete_Nadirmessung_MgII.size(); i++) {
		int Nadir_i = i - Ausgewertete_Limbmessung_MgII.size();
		Saeulendichten_MgII(i)
			= Ausgewertete_Nadirmessung_MgII[Nadir_i].m_Zeilendichte;
		Saeulendichten_Fehler_MgII(i)
			= Ausgewertete_Nadirmessung_MgII[Nadir_i].m_Fehler_Zeilendichten;
	}
	if (mache_volles_Retrieval_MgII == "ja")
		generate_Sy(S_y_MgII, Saeulendichten_Fehler_MgII);
	//Ende Säulendichten und Fehler auffüllen
	////////////////////////////////////////////////////////////////////////////
	// ENDE Spezies Mg II //
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// Spezies unbekannte //
	////////////////////////////////////////////////////////////////////////////
	// einen Vektor der Dichte_n
	MPL_Matrix Dichte_n_unknown(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_n_unknown.Null_Initialisierung();
	//Dichte_n_MgI.in_Datei_speichern("/tmp/mlangowski/0/Dichte_n_MgI");
	// einen Vektor x_a für die a-priori-Lösung der Dichte (oder Startwert usw)
	MPL_Matrix Dichte_apriori_unknown(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_apriori_unknown.Null_Initialisierung();  // Null als Startwert
	//quadratische matrix
	MPL_Matrix S_letzte_Hoehe_unknown(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// Vektor mit Zeilendichten  für alle Messungen einer Spezies
	MPL_Matrix Saeulendichten_unknown(Ausgewertete_Limbmessung_unknown.size()
			+ Ausgewertete_Nadirmessung_unknown.size(), 1); //Spaltenvektor
	MPL_Matrix Saeulendichten_Fehler_unknown(Ausgewertete_Limbmessung_unknown.size()
			+ Ausgewertete_Nadirmessung_unknown.size(), 1); //Spaltenvektor

	// Spezies Index for unknown
	spez_index = 2;

	double unknown_Lambda_letzte_Hoehe = 1E-32; //100 für ein

	// take the lambdas from the config file, no need to re-compile
	double unknown_Lambda_apriori
		= Konf.m_Retrieval_Kovarianzen[spez_index + 0 * Konf.m_Anzahl_der_Emitter];
	double unknown_Lambda_Breite
		= Konf.m_Retrieval_Kovarianzen[spez_index + 1 * Konf.m_Anzahl_der_Emitter];
	double unknown_Lambda_Hoehe
		= Konf.m_Retrieval_Kovarianzen[spez_index + 2 * Konf.m_Anzahl_der_Emitter];

	cout << "Retrieval Kovarianzen unknown:" << endl;
	cout << "L_apriori = " << unknown_Lambda_apriori << ", "
		 << "L_Breite = " << unknown_Lambda_Breite << ", "
		 << "L_Höhe = " << unknown_Lambda_Hoehe << endl;

	// Messfehlermatrix S_y y*y
	MPL_Matrix S_y_unknown(Saeulendichten_unknown.m_Elementanzahl,
			Saeulendichten_unknown.m_Elementanzahl); // quadratische matrix
	MPL_Matrix S_apriori_unknown(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// AMF_Matrix erstellen // Zeilen wie y spalten wie x
	MPL_Matrix AMF_unknown(Saeulendichten_unknown.m_Elementanzahl, Grid.m_Anzahl_Punkte);
	// S_Breite, S_Hoehe, S_Apriori, S_y, AMF aufbauen
	IERR = 0;
	//cerr<<"MgI Matrizen aufbaun\n";
	if (mache_volles_Retrieval_unknown == "ja") {
		Matrizen_Aufbauen(S_Breite, S_Hoehe, S_letzte_Hoehe_unknown,
						   unknown_Lambda_letzte_Hoehe,
						   S_apriori_unknown, S_y_unknown,
						   AMF_unknown, unknown_Lambda_apriori,
						   // TODO 1 ist MgI, 0 ist MgII hier wieder korrigieren
						   Spezies_Fenster[spez_index],
						   Grid,
						   Ausgewertete_Limbmessung_unknown,
						   Ausgewertete_Nadirmessung_unknown,
						   Konf, IERR);
		if (IERR != 0) {
			cerr << "Fehler bei Matrizen_Aufbauen\n";
			return -1; //Hauptprogramm beenden
		}
	}
	// Säulendichten und Fehler auffüllen (Fehler für Wichtungsmatrixberechnung)
	// Limb unknown
	//cerr<<"unknown Limb\n";
	for (unsigned int i = 0; i < Ausgewertete_Limbmessung_unknown.size(); i++) {
		Saeulendichten_unknown(i) = Ausgewertete_Limbmessung_unknown[i].m_Zeilendichte;
		Saeulendichten_Fehler_unknown(i)
			= Ausgewertete_Limbmessung_unknown[i].m_Fehler_Zeilendichten;
	}
	// Nadir unknown
	//cerr<<"unknown Nadir\n";
	for (unsigned int i = Ausgewertete_Limbmessung_unknown.size();
			i < Ausgewertete_Limbmessung_unknown.size() + Ausgewertete_Nadirmessung_unknown.size(); i++) {
		int Nadir_i = i - Ausgewertete_Limbmessung_unknown.size();
		Saeulendichten_unknown(i)
			= Ausgewertete_Nadirmessung_unknown[Nadir_i].m_Zeilendichte;
		Saeulendichten_Fehler_unknown(i)
			= Ausgewertete_Nadirmessung_unknown[Nadir_i].m_Fehler_Zeilendichten;
	}
	if (mache_volles_Retrieval_unknown == "ja")
		generate_Sy(S_y_unknown, Saeulendichten_Fehler_unknown);
	////////////////////////////////////////////////////////////////////////////
	// ENDE Spezies unknown
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// Spezies FeI //
	////////////////////////////////////////////////////////////////////////////
	// einen Vektor der Dichte_n
	MPL_Matrix Dichte_n_FeI(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_n_FeI.Null_Initialisierung();
	//Dichte_n_MgI.in_Datei_speichern("/tmp/mlangowski/0/Dichte_n_MgI");
	// einen Vektor x_a für die a-priori-Lösung der Dichte (oder Startwert usw)
	MPL_Matrix Dichte_apriori_FeI(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_apriori_FeI.Null_Initialisierung();  // Null als Startwert
	//quadratische matrix
	MPL_Matrix S_letzte_Hoehe_FeI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// Vektor mit Zeilendichten  für alle Messungen einer Spezies
	MPL_Matrix Saeulendichten_FeI(Ausgewertete_Limbmessung_FeI.size()
			+ Ausgewertete_Nadirmessung_FeI.size(), 1); //Spaltenvektor
	MPL_Matrix Saeulendichten_Fehler_FeI(Ausgewertete_Limbmessung_FeI.size()
			+ Ausgewertete_Nadirmessung_FeI.size(), 1); //Spaltenvektor

	// Spezies Index for FeI
	spez_index = 3;

	double FeI_Lambda_letzte_Hoehe = 1E-32; //100 für ein

	// take the lambdas from the config file, no need to re-compile
	double FeI_Lambda_apriori
		= Konf.m_Retrieval_Kovarianzen[spez_index + 0 * Konf.m_Anzahl_der_Emitter];
	double FeI_Lambda_Breite
		= Konf.m_Retrieval_Kovarianzen[spez_index + 1 * Konf.m_Anzahl_der_Emitter];
	double FeI_Lambda_Hoehe
		= Konf.m_Retrieval_Kovarianzen[spez_index + 2 * Konf.m_Anzahl_der_Emitter];

	cout << "Retrieval Kovarianzen FeI:" << endl;
	cout << "L_apriori = " << FeI_Lambda_apriori << ", "
		 << "L_Breite = " << FeI_Lambda_Breite << ", "
		 << "L_Höhe = " << FeI_Lambda_Hoehe << endl;

	// Messfehlermatrix S_y y*y
	MPL_Matrix S_y_FeI(Saeulendichten_FeI.m_Elementanzahl,
			Saeulendichten_FeI.m_Elementanzahl); // quadratische matrix
	MPL_Matrix S_apriori_FeI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	// AMF_Matrix erstellen // Zeilen wie y spalten wie x
	MPL_Matrix AMF_FeI(Saeulendichten_FeI.m_Elementanzahl, Grid.m_Anzahl_Punkte);
	// S_Breite, S_Hoehe, S_Apriori, S_y, AMF aufbauen
	IERR = 0;
	//cerr<<"MgI Matrizen aufbaun\n";
	if (mache_volles_Retrieval_FeI == "ja") {
		Matrizen_Aufbauen(S_Breite, S_Hoehe, S_letzte_Hoehe_FeI,
						   FeI_Lambda_letzte_Hoehe,
						   S_apriori_FeI, S_y_FeI,
						   AMF_FeI, FeI_Lambda_apriori,
						   Spezies_Fenster[spez_index], // TODO Index setzten
						   Grid,
						   Ausgewertete_Limbmessung_FeI,
						   Ausgewertete_Nadirmessung_FeI,
						   Konf, IERR);
		if (IERR != 0) {
			cerr << "Fehler bei Matrizen_Aufbauen\n";
			return -1; //Hauptprogramm beenden
		}
	}
	// Säulendichten und Fehler auffüllen (Fehler für Wichtungsmatrixberechnung)
	for (unsigned int i = 0; i < Ausgewertete_Limbmessung_FeI.size(); i++) {
		Saeulendichten_FeI(i) = Ausgewertete_Limbmessung_FeI[i].m_Zeilendichte;
		Saeulendichten_Fehler_FeI(i)
			= Ausgewertete_Limbmessung_FeI[i].m_Fehler_Zeilendichten;
	}
	for (unsigned int i = Ausgewertete_Limbmessung_FeI.size();
			i < Ausgewertete_Limbmessung_FeI.size() + Ausgewertete_Nadirmessung_FeI.size(); i++) {
		int Nadir_i = i - Ausgewertete_Limbmessung_FeI.size();
		Saeulendichten_FeI(i)
			= Ausgewertete_Nadirmessung_FeI[Nadir_i].m_Zeilendichte;
		Saeulendichten_Fehler_FeI(i)
			= Ausgewertete_Nadirmessung_FeI[Nadir_i].m_Fehler_Zeilendichten;
	}
	if (mache_volles_Retrieval_FeI == "ja")
		generate_Sy(S_y_FeI, Saeulendichten_Fehler_FeI);
	////////////////////////////////////////////////////////////////////////////
	// ENDE Spezies FeI
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// Spezies NO //
	////////////////////////////////////////////////////////////////////////////
	// einen Vektor der Dichte_n
	MPL_Matrix Dichte_n_NO(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_n_NO.Null_Initialisierung();

	// einen Vektor x_a für die a-priori-Lösung der Dichte (oder Startwert usw)
	MPL_Matrix Dichte_apriori_NO(Grid.m_Anzahl_Punkte, 1); //Spaltenvektor
	Dichte_apriori_NO.Null_Initialisierung();  // Null als Startwert

	// Vektor mit Zeilendichten für alle Messungen einer Spezies
	//quadratische matrix
	MPL_Matrix S_letzte_Hoehe_NO(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix Saeulendichten_NO(Ausgewertete_Limbmessung_NO.size()
			+ Ausgewertete_Nadirmessung_NO.size(), 1); //Spaltenvektor
	MPL_Matrix Saeulendichten_Fehler_NO(Ausgewertete_Limbmessung_NO.size()
			+ Ausgewertete_Nadirmessung_NO.size(), 1); //Spaltenvektor

	// Spezies Index for NO
	spez_index = 4;

	double NO_Lambda_letzte_Hoehe = 1E-32; //100 für ein

	// take the lambdas from the config file, no need to re-compile
	double NO_Lambda_apriori
		= Konf.m_Retrieval_Kovarianzen[spez_index + 0 * Konf.m_Anzahl_der_Emitter];
	double NO_Lambda_Breite
		= Konf.m_Retrieval_Kovarianzen[spez_index + 1 * Konf.m_Anzahl_der_Emitter];
	double NO_Lambda_Hoehe
		= Konf.m_Retrieval_Kovarianzen[spez_index + 2 * Konf.m_Anzahl_der_Emitter];

	cout << "Retrieval Kovarianzen NO:" << endl;
	cout << "L_apriori = " << NO_Lambda_apriori << ", "
		 << "L_Breite = " << NO_Lambda_Breite << ", "
		 << "L_Höhe = " << NO_Lambda_Hoehe << endl;

	// Messfehlermatrix S_y y*y
	// quadratische matrix
	MPL_Matrix S_y_NO(Saeulendichten_NO.m_Elementanzahl, Saeulendichten_NO.m_Elementanzahl);
	MPL_Matrix S_apriori_NO(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);

	// AMF_Matrix erstellen
	// Zeilen wie y spalten wie x
	MPL_Matrix AMF_NO(Saeulendichten_NO.m_Elementanzahl, Grid.m_Anzahl_Punkte);
	// S_Breite, S_Hoehe, S_Apriori, S_y, AMF aufbauen
	IERR = 0;
	// cerr << "NO Matrizen aufbauen" << endl;
	if (mache_volles_Retrieval_NO == "ja") {
		Matrizen_Aufbauen(S_Breite, S_Hoehe, S_letzte_Hoehe_NO,
						   NO_Lambda_letzte_Hoehe,
						   S_apriori_NO, S_y_NO,
						   AMF_NO, NO_Lambda_apriori,
						   Spezies_Fenster[spez_index],
						   Grid,
						   Ausgewertete_Limbmessung_NO,
						   Ausgewertete_Nadirmessung_NO,
						   Konf, IERR);
		if (IERR != 0) {
			cerr << "Fehler bei Matrizen_Aufbauen\n";
			return -1; //Hauptprogramm beenden
		}
	}

	// Säulendichten und Fehler auffüllen (Fehler für Wichtungsmatrixberechnung)
	for (unsigned int i = 0; i < Ausgewertete_Limbmessung_NO.size(); i++) {
		Saeulendichten_NO(i) = Ausgewertete_Limbmessung_NO[i].m_Zeilendichte;
		Saeulendichten_Fehler_NO(i)
			= Ausgewertete_Limbmessung_NO[i].m_Fehler_Zeilendichten;
	}
	// Nadir NO
	for (unsigned int i = Ausgewertete_Limbmessung_NO.size();
			i < Ausgewertete_Limbmessung_NO.size()
				+ Ausgewertete_Nadirmessung_NO.size(); i++) {
		int Nadir_i = i - Ausgewertete_Limbmessung_NO.size();
		Saeulendichten_NO(i)
			= Ausgewertete_Nadirmessung_NO[Nadir_i].m_Zeilendichte;
		Saeulendichten_Fehler_NO(i)
			= Ausgewertete_Nadirmessung_NO[Nadir_i].m_Fehler_Zeilendichten;
	}
	if (mache_volles_Retrieval_NO == "ja")
		generate_Sy(S_y_NO, Saeulendichten_Fehler_NO);
	////////////////////////////////////////////////////////////////////////////
	// ENDE Spezies NO
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// Matrizen sind aufgebaut  // Schön wärs 7.7.2010...erstmal weiter, weil
	// es beim Raytracing noch probleme gibt
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// Matrizen sind aufgebaut
	// // wer hat an der Uhr gedreht 14.9.2010; getestet 4.10.2010
	// // jetzt auch für Mg II.... oha...schon 7.12.2010
	////////////////////////////////////////////////////////////////////////////

	time(&Teil3_End);
	T3_Dauer = Teil3_End - Teil3_Start;
	/***************************************************************************
	 ENDE TEIL 3
	***************************************************************************/

	/***************************************************************************
	 TEIL 4 AUSFÜHREN DES RETRIEVALS
	 AUFSTELLEN UND LÖSEN DER NORMALGLEICHUNG
	 **************************************************************************/
	cerr << "Teil4\n";
	time(&Teil4_Start);

	////////////////////////////////////////////////////////////////////////////
	// Iteration des Retrievals
	////////////////////////////////////////////////////////////////////////////
	/////////////////////////
	// Spezies Mg I //
	/////////////////////////
	//Dichte_n_MgI.in_Datei_speichern("/tmp/mlangowski/0/Dichte_n_MgI");
	//MPL_Matrix bla;
	//bla=Dichte_n_MgI;
	// TODO neue Version fit machen
	if (mache_volles_Retrieval_MgI == "ja") {
		switch (Konf.retrieval_algo) {
		case 0:
			Retrievaliteration_old(Dichte_n_MgI, Dichte_apriori_MgI,
							   Saeulendichten_MgI,
							   S_apriori_MgI, S_y_MgI, S_Breite, S_Hoehe,
							   S_letzte_Hoehe_MgI, MgI_Lambda_Breite,
							   MgI_Lambda_Hoehe, AMF_MgI, Konf);
		case 1:
		default:
			Retrievaliteration(Dichte_n_MgI, Dichte_apriori_MgI,
						   Saeulendichten_MgI,
						   S_apriori_MgI, S_y_MgI, S_Breite, S_Hoehe,
						   MgI_Lambda_Breite,
						   MgI_Lambda_Hoehe, AMF_MgI, Konf);
			break;
		}
	}
	/////////////////////////
	// Spezies Mg II//
	/////////////////////////
	if (mache_volles_Retrieval_MgII == "ja") {
		switch (Konf.retrieval_algo) {
		case 0:
			Retrievaliteration_old(Dichte_n_MgII, Dichte_apriori_MgII,
							   Saeulendichten_MgII,
							   S_apriori_MgII, S_y_MgII, S_Breite, S_Hoehe,
							   S_letzte_Hoehe_MgII, MgII_Lambda_Breite,
							   MgII_Lambda_Hoehe, AMF_MgII, Konf);
			break;
		case 1:
		default:
			Retrievaliteration(Dichte_n_MgII, Dichte_apriori_MgII,
						   Saeulendichten_MgII,
						   S_apriori_MgII, S_y_MgII, S_Breite, S_Hoehe,
						   MgII_Lambda_Breite,
						   MgII_Lambda_Hoehe, AMF_MgII, Konf);
			break;
		}
	}
	/////////////////////////
	// Spezies unknown//
	/////////////////////////
	if (mache_volles_Retrieval_unknown == "ja") {
		switch (Konf.retrieval_algo) {
		case 0:
			Retrievaliteration_old(Dichte_n_unknown, Dichte_apriori_unknown,
							   Saeulendichten_unknown,
							   S_apriori_unknown, S_y_unknown, S_Breite, S_Hoehe,
							   S_letzte_Hoehe_unknown, unknown_Lambda_Breite,
							   unknown_Lambda_Hoehe, AMF_unknown, Konf);
			break;
		case 1:
		default:
			Retrievaliteration(Dichte_n_unknown, Dichte_apriori_unknown,
						   Saeulendichten_unknown,
						   S_apriori_unknown, S_y_unknown, S_Breite, S_Hoehe,
						   unknown_Lambda_Breite,
						   unknown_Lambda_Hoehe, AMF_unknown, Konf);
			break;
		}
	}
	/////////////////////////
	// Spezies FeI   //
	/////////////////////////
	if (mache_volles_Retrieval_FeI == "ja") {
		switch (Konf.retrieval_algo) {
		case 0:
			Retrievaliteration_old(Dichte_n_FeI, Dichte_apriori_FeI,
							   Saeulendichten_FeI,
							   S_apriori_FeI, S_y_FeI, S_Breite, S_Hoehe,
							   S_letzte_Hoehe_FeI, FeI_Lambda_Breite,
							   FeI_Lambda_Hoehe, AMF_FeI, Konf);
			break;
		case 1:
		default:
			Retrievaliteration(Dichte_n_FeI, Dichte_apriori_FeI,
						   Saeulendichten_FeI,
						   S_apriori_FeI, S_y_FeI, S_Breite, S_Hoehe,
						   FeI_Lambda_Breite,
						   FeI_Lambda_Hoehe, AMF_FeI, Konf);
			break;
		}
	}
	/////////////////////////
	// Spezies NO
	/////////////////////////
	if (mache_volles_Retrieval_NO == "ja") {
		switch (Konf.retrieval_algo) {
		case 0:
			Retrievaliteration_old(Dichte_n_NO, Dichte_apriori_NO,
						   Saeulendichten_NO,
						   S_apriori_NO, S_y_NO, S_Breite, S_Hoehe,
						   S_letzte_Hoehe_NO,
						   NO_Lambda_Breite,
						   NO_Lambda_Hoehe, AMF_NO, Konf);
			break;
		case 1:
		default:
			Retrievaliteration(Dichte_n_NO, Dichte_apriori_NO,
						   Saeulendichten_NO,
						   S_apriori_NO, S_y_NO, S_Breite, S_Hoehe,
						   NO_Lambda_Breite,
						   NO_Lambda_Hoehe, AMF_NO, Konf);
			break;
		}
	}


	time(&Teil4_End);
	T4_Dauer = Teil4_End - Teil4_Start;
	/***************************************************************************
	TEIL 5 FEHLERABSCHÄTZUNG
	BERECHNUNG DER AVERAGING-KERNEL MATRIX UND DER FEHLERMATRIX
	***************************************************************************/
	cerr << "Teil5\n";
	time(&Teil5_Start);

	// todolist
	// die LHS Seite der Normalgleichung berechnen
	// zum Invertieren 1 Matrix anhängen und Diagonalisieren -> Fehlermatrix
	// die Averaging Kernels erhält man mit Multiplikation der Fehlermatrix mit
	// AMF^T S_y AMF
	//Retrievalfehler_Abschaetzung();

	/////////////////////////
	// Spezies MgI //
	/////////////////////////
	MPL_Matrix S_x_MgI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_x_meas_MgI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	//Averaging Kernel Matrix
	MPL_Matrix AKM_MgI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	if (mache_volles_Retrieval_MgI == "ja") {
		Retrievalfehler_Abschaetzung(S_x_MgI, S_x_meas_MgI, AKM_MgI,
									 S_apriori_MgI, S_y_MgI, S_Breite, S_Hoehe,
									 S_letzte_Hoehe_MgI, MgI_Lambda_Breite,
									 MgI_Lambda_Hoehe, AMF_MgI, Konf);
	}
	/////////////////////////
	// Spezies MgII //
	/////////////////////////
	MPL_Matrix S_x_MgII(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_x_meas_MgII(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	//Averaging Kernel Matrix
	MPL_Matrix AKM_MgII(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	if (mache_volles_Retrieval_MgII == "ja") {
		Retrievalfehler_Abschaetzung(S_x_MgII, S_x_meas_MgII, AKM_MgII,
									 S_apriori_MgII, S_y_MgII, S_Breite, S_Hoehe,
									 S_letzte_Hoehe_MgII, MgII_Lambda_Breite,
									 MgII_Lambda_Hoehe, AMF_MgII, Konf);
	}
	/////////////////////////
	// Spezies unknown //
	/////////////////////////
	MPL_Matrix S_x_unknown(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_x_meas_unknown(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	//Averaging Kernel Matrix
	MPL_Matrix AKM_unknown(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	if (mache_volles_Retrieval_unknown == "ja") {
		Retrievalfehler_Abschaetzung(S_x_unknown, S_x_meas_unknown, AKM_unknown,
									 S_apriori_unknown, S_y_unknown, S_Breite,
									 S_Hoehe, S_letzte_Hoehe_unknown,
									 unknown_Lambda_Breite,
									 unknown_Lambda_Hoehe, AMF_unknown, Konf);
	}
	/////////////////////////
	// Spezies FeI  //
	/////////////////////////
	MPL_Matrix S_x_FeI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_x_meas_FeI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	//Averaging Kernel Matrix
	MPL_Matrix AKM_FeI(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	if (mache_volles_Retrieval_FeI == "ja") {
		Retrievalfehler_Abschaetzung(S_x_FeI, S_x_meas_FeI, AKM_FeI,
									 S_apriori_FeI, S_y_FeI, S_Breite, S_Hoehe,
									 S_letzte_Hoehe_FeI, FeI_Lambda_Breite,
									 FeI_Lambda_Hoehe, AMF_FeI, Konf);
	}
	/////////////////////////
	// Spezies NO //
	/////////////////////////
	MPL_Matrix S_x_NO(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	MPL_Matrix S_x_meas_NO(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	//Averaging Kernel Matrix
	MPL_Matrix AKM_NO(Grid.m_Anzahl_Punkte, Grid.m_Anzahl_Punkte);
	if (mache_volles_Retrieval_NO == "ja") {
		Retrievalfehler_Abschaetzung(S_x_NO, S_x_meas_NO, AKM_NO,
									 S_apriori_NO, S_y_NO, S_Breite, S_Hoehe,
									 S_letzte_Hoehe_NO, NO_Lambda_Breite,
									 NO_Lambda_Hoehe, AMF_NO, Konf);
	}
	time(&Teil5_End);
	T5_Dauer = Teil5_End - Teil5_Start;
	/***************************************************************************
	ENDE TEIL 5
	 **************************************************************************/

	/***************************************************************************
	TEIL 6 ERGEBNISAUSGABE
	DIE ERGEBNISAUSGABE SCHLIEßT DAS PROGRAMM AB
	***************************************************************************/
	cerr << "Teil6\n";
	time(&Teil6_Start);
	//int Ausgabe_Dichten(string Dateiname_out,Retrievalgitter Grid,
	//  MPL_Matrix Dichten, MPL_Matrix S_x, MPL_Matrix AKM)
	string Dateiname_out;
	// MgI /////////////////////
	Dateiname_out = Arbeitsverzeichnis + "/" + sssss_MgI + Dateiout_Mittelteil;
	//cout<<"Dateiname_out: "<<Dateiname_out<<"\n";
	//cout<<"Dateiout_Mittelteil: "<<Dateiout_Mittelteil<<"\n";
	if (mache_volles_Retrieval_MgI == "ja") {
		Ausgabe_Dichten(Dateiname_out, Grid, Dichte_n_MgI, Dichte_n_tot,
				Dichte_apriori_MgI, S_x_MgI, S_x_meas_MgI, AKM_MgI);
	}
	//////////////////////////////
	// MgII ////////////////////
	Dateiname_out = Arbeitsverzeichnis + "/" + sssss_MgII + Dateiout_Mittelteil;
	//cout<<"Dateiname_out: "<<Dateiname_out<<"\n";
	if (mache_volles_Retrieval_MgII == "ja") {
		Ausgabe_Dichten(Dateiname_out, Grid, Dichte_n_MgII, Dichte_n_tot,
				Dichte_apriori_MgII, S_x_MgII, S_x_meas_MgII, AKM_MgII);
	}
	//////////////////////////////
	// unknown//////////////
	Dateiname_out = Arbeitsverzeichnis + "/" + sssss_unknown + Dateiout_Mittelteil;
	//cout<<"Dateiname_out: "<<Dateiname_out<<"\n";
	if (mache_volles_Retrieval_unknown == "ja") {
		Ausgabe_Dichten(Dateiname_out, Grid, Dichte_n_unknown, Dichte_n_tot,
				Dichte_apriori_unknown,
				S_x_unknown, S_x_meas_unknown, AKM_unknown);
	}
	// FeI//////////////
	Dateiname_out = Arbeitsverzeichnis + "/" + sssss_FeI + Dateiout_Mittelteil;
	//cout<<"Dateiname_out: "<<Dateiname_out<<"\n";
	if (mache_volles_Retrieval_FeI == "ja") {
		Ausgabe_Dichten(Dateiname_out, Grid, Dichte_n_FeI, Dichte_n_tot,
				Dichte_apriori_FeI, S_x_FeI, S_x_meas_FeI, AKM_FeI);
	}
	//////////////////////////////
	// NO ////////////////////
	Dateiname_out = Arbeitsverzeichnis + "/" + sssss_NO + Dateiout_Mittelteil;
	if (mache_volles_Retrieval_NO == "ja") {
		MPL_Matrix result = AMF_NO * Dichte_n_NO;
		Ausgabe_Dichten(Dateiname_out, Grid, Dichte_n_NO, Dichte_n_tot,
				Dichte_apriori_NO, S_x_NO, S_x_meas_NO, AKM_NO);
		Ausgabe_Saeulendichten_back(Pfad_Saeulen_Limb_NO_back,
				Ausgewertete_Limbmessung_NO, result);
	}

	time(&Teil6_End);
	T6_Dauer = Teil6_End - Teil6_Start;
	/***************************************************************************
	ENDE TEIL 6
	***************************************************************************/
	int Messungen_pro_Limbdatei = 7;
	if (untersuche_limb_mesothermo_states == "ja") {
		Messungen_pro_Limbdatei = 25;
	}
	cout << "Aussortierte Limbmessungen:\n";
	cout << "Limb Nachtmessungen: " << Messungen_pro_Limbdatei *counter_Nachtmessungen << "\n";
	cout << "Limb Messungen mit NLCs:" << Messungen_pro_Limbdatei *counter_NLC_detektiert << "\n";
	cout << "Limb Messungen mit ungenauen Koordinaten: " << counter_Richtungsvektor_nicht_ok << "\n";
	cout << "Nadir Nachtmessungen: " << counter_Nachtmessungen_Nadir << "\n";
	cout << "Nadir Nachtdateien: " << counter_Nadir_Nacht_Dateien << "\n";
	////////////////////////////////////////////////////////////////////////////
	// Zeitmessung
	time(&timer1);
	deltaT = timer1 - start_zeit;
	stringstream buf;
	buf << "Gesamtlaufzeit des Programms in Sekunden:\t " << deltaT;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Dauer Teilprozess 1 Vorbereitung:\t\t " << T1_Dauer;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Dauer Teilprozess 2 Säulendichtenbestimmung:\t " << T2_Dauer;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Dauer Teilprozess 3 Aufbau des Gitters:\t\t " << T3_Dauer;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Zeit steckt in Teil 3 im wesentlichen in den Raytracing-Schleifen.";
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Dauer Teilprozess 4 Dichte-Retrieval:\t\t " << T4_Dauer;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Dauer Teilprozess 5 Fehlerabschätzung:\t\t " << T5_Dauer;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	buf.str(string());
	buf << "Dauer Teilprozess 6 Ausgabe:\t\t\t " << T6_Dauer;
	Nachricht_Schreiben(buf.str(), 10, Prioritylevel);
	////////////////////////////////////////////////////////////////////////////
	Nachricht_Schreiben("Beende Programm regulär...", 3, Prioritylevel);
	return 0;
} // Ende Hauptprogramm
