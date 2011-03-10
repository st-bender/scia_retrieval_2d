/*
 * Sonnenspektrum.h
 *
 *  Created on: 20.04.2010
 *      Author: martin
 */

#ifndef SONNENSPEKTRUM_HH_
#define SONNENSPEKTRUM_HH_

#include <string>
#include "Messung_Limb.h"

using namespace std;

class Sonnenspektrum
{
    public:
        int Laden_GOME(string Dateiname, string Fallback_Dateiname);  // Diese Funktion ist für die Sonnenspektren von GOME von
                                                                                             // Mark Weber
        int Laden_SCIA(string Dateiname, string Fallback_Dateiname); //Sciamachy Sonnenspektrum des Orbits

        int Interpolieren(Messung_Limb& Messung_Erdschein);
        int nicht_interpolieren();
        double m_Wellenlaengen[850]; // wie lang sind die eigentlich
        double m_Intensitaeten[850];
        // Auf Erdschein-Messungs-Wellenlängen Interpolierte Intensität (Wellenlängen sind dann gleich dem ErdscheinWL)
        double m_Int_interpoliert[826];
        double m_WL_interpoliert[826]; //für Ausgabe
        int m_Anzahl_WL;
        //Wartungsfunktion
        int Speichern(string Dateiname);  //zur Kontrolle
        int Speichern_was_geladen_wurde(string Dateiname);
};


#endif /* SONNENSPEKTRUM_HH_ */
