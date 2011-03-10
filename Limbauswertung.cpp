/*
 * Limbauswertung.cpp
 *
 *  Created on: 28.04.2010
 *      Author: martin
 */

/*************************************************************************
 Was passiert hier:
 Die Limbdateien werden aus der Orbitliste geladen.

 Die Interpolation des Sonnenspektrums auf das Limbmessungsspektrum muss nur einmal erfolgen.
 Da das Limbgitter vom Anfang her dem Nadirgitter entspricht, geschiet dies insgesamt auch nur einmal. Deshalb muss die
 erste Datei aus der Orbitliste eine Limbdatei sein.

 Die gemessene Intensität ist mit der Zeilendichte linear über den sogenannten Emissivitätsfaktor, oder auch g-Factor verbunden.
*************************************************************************/

#include <vector>
#include "Orbitliste.h"
#include "Sonnenspektrum.h"
#include "Speziesfenster.h"
#include "Ausgewertete_Messung_Limb.h"
#include "Messung_Limb.h"
#include "Datei_IO.h"  //ReadL1C_Limb
#include "Messungs_ausschliessende_Tests.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Limb_Auswertung
//////////////////////////////////////////////////////////////////////////////////////
int Limb_Auswertung(Orbitliste Orbitlist,
                                  int l,
                                  Sonnenspektrum& Solspec,
                                  vector<Speziesfenster>& Spezies_Fenster,
                                  int* counter_Nachtmessungen,
                                  int* counter_NLC_detektiert,
                                  int* counter_Richtungsvektor_nicht_ok,
                                  string Arbeitsverzeichnis, string mache_Fit_Plots,
                                  string limb_meso_thermo,       // "ja" oder "nein"
                                  vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_MgI,
                                  vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_MgII,
                                  vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_unknown,
                                  vector<Ausgewertete_Messung_Limb>& Ausgewertete_Limbmessung_FeI)
{
    unsigned int i,j,k;
    //Einmalig die Rohdaten aus der Datei Laden
    vector<Messung_Limb> Rohdaten;  // Achtung das ist noch nicht der entgültige Vektor, weil dieser noch um die Speziesparameter ergänzt werden muss
    Messung_Limb Tropo;  //Hier stecken NUR die Intensitäten der troposphärischen Säule drin....bzw der tiefsten bei Limb_meso_thermo
    Messung_Limb mean_10_20; //Hier stecken NUR die Intensitäten der Säulen 10 bis 20 (ca 30-60km) drin
    //cerr<<"Rohdaten einlesen\n";
    if(limb_meso_thermo!="ja")
    {    Rohdaten=ReadL1C_Limb_mpl_binary(Orbitlist.m_Dateinamen[l], Tropo,mean_10_20);    }
    else
    {    //Rohdaten=ReadL1C_Limb_meso_thermo_mpl_binary(Orbitlist.m_Dateinamen[l], Tropo);
        int Anzahl_Hoehen=25;
        Rohdaten=ReadL1C_Limb_meso_thermo_mpl_binary_reduziert(Orbitlist.m_Dateinamen[l], Tropo, Anzahl_Hoehen);
    }
    //cerr<<Orbitlist.m_Dateinamen[l]<<" wird bearbeitet\n";
    //if(l==0) // einmal das Sonnenspektrum anpassen/interpolieren
    //{
        //Testen, ob ReadL1C ordentlich gearbeitet hat
        //Rohdaten[0].Ausgabe_in_Datei("CHECKDATA/Rohdaten_erste_Limb_Messung.txt"); -> Das geht jetzt
        Solspec.Interpolieren(Rohdaten[0]);
        //Solspec.nicht_interpolieren();

        // Testen, ob die Interpolation erfolgreich war
        //Solspec.Speichern("CHECKDATA/Sonne_interpoliert_auf_826.txt"); ->ok
    //}
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Hier wäre ein guter Ort, um zu prüfen, ob die Rohdaten weiter verwendet werden dürfen
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bool ist_Nachtmessung;
    bool NLC_detektiert;
    if(limb_meso_thermo!="ja")
    {    Test_auf_Nachtmessung_Limb(Tropo, ist_Nachtmessung);                          }
    else
    {    Test_auf_Nachtmessung_Limb_meso_thermo(Tropo, ist_Nachtmessung);    }
    if(ist_Nachtmessung==true)
    {
       // cout<<"counter_Nachtmessungen vorher: "<<*counter_Nachtmessungen <<"\n";
        (*counter_Nachtmessungen)++;
       // cout<<"counter_Nachtmessungen nacher: "<<*counter_Nachtmessungen <<"\n";
        return 1;  //Nachtmessung 1
    }
    /*Test_auf_NLC_Limb(Rohdaten, NLC_detektiert);
    if(NLC_detektiert==true)
    {
        (*counter_NLC_detektiert)++;
        return 2;  //NLC 2
    }*/
    Test_auf_korrekte_geolocations_Limb(Rohdaten,counter_Richtungsvektor_nicht_ok);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // Jedes Element von Rohdaten wird nun in seperaten Schritten mit den Speziesdaten aufgefüllt
    // und die Bestimmung der Zeilendichte wird durchgeführt
    for(i=0;i<Rohdaten.size();i++) //Schleife über alle Rohdaten
    {
        Messung_Limb Messung=Rohdaten[i]; //Übergabe der Rohdaten, sodass diese nicht verändert werden->wegen Parallelisierung
        // sin und cos sind langsame Funktionen...aber werden hierbei auch nicht oft eingesetzt
        Messung.Deklinationswinkel_bestimmen(); //Sonnenlattitude
        Messung.Sonnen_Longitude_bestimmen();
        Messung.Intensitaeten_normieren(Solspec.m_Int_interpoliert);  // m_Intensitaeten enthält nun nichtmehr I sondern I/(piF)  // Das könnte man auch nur für die Par Fenster durchführen
        //cout<<Messung.m_Intensitaeten_durch_piF[0]<<"\n";
        for(j=0;j<Spezies_Fenster.size();j++) //Schleife über alle Spezies wie z.b. Mg oder Mg+
        {
            for(k=0;k<Spezies_Fenster[j].m_Wellenlaengen.size();k++)//Schleife über alle Linien dieser Spezies
            {
                // Aus SZA_TP und SAA_TP lässt sich die Polararisation in den Liniendaten des Speziesfensters
                // ermitteln, und damit die Emissivität berechnen
                //Spezfenst.m_Liniendaten[k].m_theta=Messung.m_Streuwinkel;  der Streuwinkel muss woanders berechnet werden
                // Die Phasenfunktion, steckt so nun in den Slant Coloumns drin
                //cout<<Spezfenst.m_Liniendaten[k].m_theta<<"\n";
                Spezies_Fenster[j].m_Liniendaten[k].Emissivitaet_ermitteln();
                //cout<<Spezfenst.m_Basisfenster_rechts_WLmin[k]<<"\n";
                //cout<<"j:"<<j<<"\t"<<"k:"<<k<<"\t"<<Spezfenst.m_Wellenlaengen[k]<<"\t"<<Spezfenst.m_Liniendaten[k].m_E1<<"\t"<<Spezfenst.m_Liniendaten[k].m_Gamma<<"\n";
                Messung.Intensitaeten_durch_piF_Gamma_berechnen(Spezies_Fenster[j],k);
                Messung.Intensitaeten_durch_piF_Gamma_mal_Gitterabstand_berechnen(Spezies_Fenster[j],k);// In der Formel ist piF in W/(m^2*Wellenlänge) verlangt..also muss noch mit der Kanalbreite multipliziert werden
                //cout<<Messung.m_Intensitaeten_durch_piF_Gamma[1]<<"\n";// ->hmm 4e7...keine Ahnung ob das Sinn macht

                // Jetzt Zeilendichte und Fehler bestimmen
//                if (Spezies_Fenster[j].m_Spezies_Name=="MgI")
//                {   //Messung.Saeulendichte_Bestimmen_MgI285nm(Spezies_Fenster[j],k,Arbeitsverzeichnis,mache_Fit_Plots, mean_10_20.m_Intensitaeten);
//                    Messung.Plots_der_Spektren_erzeugen(Spezies_Fenster[j],k,Arbeitsverzeichnis,mache_Fit_Plots, mean_10_20.m_Intensitaeten);
//                }
//                else
//                {
                    Messung.Zeilendichte_Bestimmen(Spezies_Fenster[j],k, Arbeitsverzeichnis, mache_Fit_Plots); // Hmm hier gibts noch Diskussionsbedarf
//                }

                // Zu Testzwecken fertige Messung in Datei Speichern
                //if((k==0) && (j==0)&&(i==0))
                //{
                //    Messung.Ausgabe_in_Datei("CHECKDATA/Messung_Limb_Fenster0_Hoehe_74km_0teLinie.txt");
                //}
                // Todo...wieder rausnehmen
                //if(Messung.m_Lattidude_TP<5.0) {   continue;    }

                // Ergebnis zusammenfassen
                Ausgewertete_Messung_Limb Ergebnis=Messung.Ergebnis_Zusammenfassen();
                Ergebnis.m_Wellenlaenge=Spezies_Fenster[j].m_Wellenlaengen[k];      // Die braucht man später für die Luftmassenmatrix
                //Ergebnis.Ausgabe_auf_Bildschirm();
                // Zusammenfassung der Zwischenresultate dem Vektor für die jeweilige Spezies zuordnen
                if(Spezies_Fenster[j].m_Spezies_Name=="MgI")
                {   //TODO negative Werte zulassen
                    if(Ergebnis.m_Zeilendichte>0)
                    {                        Ausgewertete_Limbmessung_MgI.push_back(Ergebnis);              }
                    else
                    {    Ergebnis.m_Zeilendichte=0;
                         Ergebnis.m_Fehler_Zeilendichten=Ergebnis.m_Fehler_Zeilendichten;
                         Ausgewertete_Limbmessung_MgI.push_back(Ergebnis);
                    }
                    //cout<<Ausgewertete_Limbmessung_MgI[0].m_Zeilendichte<<"\n";//-das geht
                }
                if(Spezies_Fenster[j].m_Spezies_Name=="MgII")
                {
                    Ausgewertete_Limbmessung_MgII.push_back(Ergebnis);
                }
                if(Spezies_Fenster[j].m_Spezies_Name=="unknown")
                {
                    Ausgewertete_Limbmessung_unknown.push_back(Ergebnis);
                }
                if(Spezies_Fenster[j].m_Spezies_Name=="FeI") { Ausgewertete_Limbmessung_FeI.push_back(Ergebnis); }

            }//ende k Linie
        }//ende j Spezies_Fenster
    }//ende i Rohdaten
    // Die Rohdaten werden nun nichtmehr gebraucht, sodass alle dynamisch allokierten Felder gelöscht werden
    // können

    Rohdaten.resize(0);
    return 0; // keine Probleme
}
//////////////////////////////////////////////////////////////////////////////////////
// ENDE Limb_Auswertung
//////////////////////////////////////////////////////////////////////////////////////
