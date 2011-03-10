/*
 * Spezies.h
 *
 *  Created on: 13.04.2010
 *      Author: martin
 */

/*******************************************************************
Objekte dieser Klasse enthalten die Wellenlängenfenster, sowie Liniendaten
für jede Linie
Achtung die Liniendaten sind abhängig vom Streuwinkel, sodass die Pasenfunktion,
sowie die Emissivitäten für jede Messung neu bestimmt werden müssen.
 *******************************************************************/

#ifndef SPEZIESFENSTER_HH_
#define SPEZIESFENSTER_HH_

#include"Liniendaten.h"
#include<vector>
#include<string>

using namespace std;

class Speziesfenster
{
public:
	// Funktionen
	void clear();

	//Membervariablen *********************************************************/
	string m_Spezies_Name;                        // z.B. Eisen I oder Eisen II
	// Teile für Zeilendichte Berechnung
	vector<double> m_Wellenlaengen;  // Vektor mit allen zugehörigen Linien
	vector<Liniendaten> m_Liniendaten;
	vector<double> m_Basisfenster_links_WLmin;
	vector<double> m_Basisfenster_links_WLmax;
	vector<double> m_Basisfenster_rechts_WLmin;
	vector<double> m_Basisfenster_rechts_WLmax;
	vector<double> m_Peakfenster_WLmin;
	vector<double> m_Peakfenster_WLmax;
	double m_FWHM;  // FWHM der Peaks, oder zumindest Startwert
	//Für alle Linien erstmal ein Vektor
	vector<string>  m_Liste_der_Plot_Dateinamen;
	//Membervariablen ende ****************************************************/
};

inline void Speziesfenster::clear()
{
	// Die Vektorengrößen auf 0 zurücksetzen
	m_Wellenlaengen.resize(0);
	m_Basisfenster_links_WLmin.resize(0);
	m_Basisfenster_links_WLmax.resize(0);
	m_Basisfenster_rechts_WLmin.resize(0);
	m_Basisfenster_rechts_WLmax.resize(0);
	m_Peakfenster_WLmin.resize(0);
	m_Peakfenster_WLmax.resize(0);
	m_Liniendaten.resize(0);
	m_Liste_der_Plot_Dateinamen.resize(0);
}


#endif /* SPEZIESFENSTER_HH_ */
