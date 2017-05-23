/*
 * Glaetten_2D.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2011 Martin Langowski
 *
 * Initial version created on: 24.01.2011
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include <string>
#include <vector>
#include "Ausgewertete_Messung_Limb.h"
#include <iostream>
#include <cstdio>

using std::cerr;
using std::string;
using std::vector;

int SCD_Glaettung(vector<Ausgewertete_Messung_Limb>& AML, int Anzahl_Linien,
		string limbmesothermo)
{
	vector<Ausgewertete_Messung_Limb> AML_OLD = AML;
	int Anzahl_Hoehen;
	if (limbmesothermo == "ja")    {
		Anzahl_Hoehen = 25;/*25*/;
	} else {
		Anzahl_Hoehen = 7;
	}
	int Anzahl_Messungen = AML.size();
	if (Anzahl_Messungen % Anzahl_Hoehen != 0) {
		cerr << "keine SCD Glaettung Anzahl Breiten kein natürliche Zahl....\n";
		return 1;
	}
	for (unsigned int i = 0; i < AML.size(); i++) {
		AML[i].m_Zeilendichte = 0;       // mit 0 resetten
	}
	//cerr<<"Anzahl_Hoehen: "<<Anzahl_Hoehen<<"\n";
	//cerr<<"Anzahl_Linien: "<<Anzahl_Linien<<"\n";
	//cerr<<"Anzahl_Messungen: "<<Anzahl_Messungen<<"\n";

	//Es wird über die direkten Nachbarn gemittelt
	// Besondere Punkte sind:
	// Die ersten Anzahl_Hoehen und die letzten Anzahl_Hoehen Messungen
	// die 0te und die Anzahl_Hoehen_te Messung
	// Zuerst 4 Ecken Sonderbehandlung, dann 4 Kanten und schließlich den Rest
	// der Matrix
	for (int k = 0; k < Anzahl_Linien; k++) {
		// 4 ECKEN ///////////////////////
		int Messung_Nr;
		// ECKE 1
		Messung_Nr = 0 + k;
		AML[Messung_Nr].m_Zeilendichte =
			(3.0 * AML_OLD[Messung_Nr].m_Zeilendichte
			 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Zeilendichte
			 + AML_OLD[Messung_Nr
			 + Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
		//cerr<<"Ecke 1\n";
		//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr+Anzahl_Linien<<" "
		//    <<Messung_Nr+Anzahl_Linien*Anzahl_Hoehen<<"\n";

		AML[Messung_Nr].m_Fehler_Zeilendichten =
			(3.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
		// ECKE 2
		Messung_Nr = (Anzahl_Hoehen - 1) * Anzahl_Linien + k;
		AML[Messung_Nr].m_Zeilendichte =
			(3.0 * AML_OLD[Messung_Nr].m_Zeilendichte
			 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Zeilendichte
			 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
		//cerr<<"Ecke 2\n";
		//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr-Anzahl_Linien<<" "
		//    <<Messung_Nr+Anzahl_Linien*Anzahl_Hoehen<<"\n";
		AML[Messung_Nr].m_Fehler_Zeilendichten =
			(3.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
		// ECKE 3
		Messung_Nr = Anzahl_Messungen - Anzahl_Linien * Anzahl_Hoehen + k;
		AML[Messung_Nr].m_Zeilendichte =
			(3.0 * AML_OLD[Messung_Nr].m_Zeilendichte
			 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Zeilendichte
			 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
		//cerr<<"Ecke 3\n";
		//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr+Anzahl_Linien<<" "
		//    <<Messung_Nr-Anzahl_Linien*Anzahl_Hoehen<<"\n";
		AML[Messung_Nr].m_Fehler_Zeilendichten =
			(3.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
		// ECKE 4
		Messung_Nr = Anzahl_Messungen - Anzahl_Linien + k;
		AML[Messung_Nr].m_Zeilendichte =
			(3.0 * AML_OLD[Messung_Nr].m_Zeilendichte
			 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Zeilendichte
			 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
		//cerr<<"Ecke 4\n";
		//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr-Anzahl_Linien<<" "
		//    <<Messung_Nr-Anzahl_Linien*Anzahl_Hoehen<<"\n";

		AML[Messung_Nr].m_Fehler_Zeilendichten =
			(3.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Fehler_Zeilendichten
			 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
		// 4 KANTEN ///////////////////////
		for (int u = 1; u < (Anzahl_Hoehen - 1); u++) {
			Messung_Nr = u * Anzahl_Linien + k;
			AML[Messung_Nr].m_Zeilendichte =
				(2.0 * AML_OLD[Messung_Nr].m_Zeilendichte
				 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Zeilendichte
				 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Zeilendichte
				 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
			//cerr<<"Kante 1\n";
			//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr +Anzahl_Linien
			//  <<" "<<Messung_Nr-Anzahl_Linien<<" "
			//  <<Messung_Nr+Anzahl_Linien*Anzahl_Hoehen<<"\n";
			AML[Messung_Nr].m_Fehler_Zeilendichten =
				(2.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr
				 + Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
			Messung_Nr = Anzahl_Messungen - Anzahl_Hoehen * Anzahl_Linien
						+ u * Anzahl_Linien + k;
			AML[Messung_Nr].m_Zeilendichte =
				(2.0 * AML_OLD[Messung_Nr].m_Zeilendichte
				 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Zeilendichte
				 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Zeilendichte
				 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
			//cerr<<"Kante 2\n";
			//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr +Anzahl_Linien
			//  <<" "<<Messung_Nr-Anzahl_Linien<<" "
			//  <<Messung_Nr-Anzahl_Linien*Anzahl_Hoehen<<"\n";
			AML[Messung_Nr].m_Fehler_Zeilendichten =
				(2.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
		} //ende for u
		for (int u = Anzahl_Hoehen * Anzahl_Linien;
				u < (Anzahl_Messungen - Anzahl_Hoehen * Anzahl_Linien);
				u += Anzahl_Hoehen * Anzahl_Linien) {
			Messung_Nr = u + k;
			AML[Messung_Nr].m_Zeilendichte =
				(2.0 * AML_OLD[Messung_Nr].m_Zeilendichte
				 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Zeilendichte
				 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte
				 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
			//cerr<<"Kante 3\n";
			//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr +Anzahl_Linien
			//  <<" "<<Messung_Nr-Anzahl_Linien*Anzahl_Hoehen<<" "
			//  <<Messung_Nr+Anzahl_Linien*Anzahl_Hoehen<<"\n";
			AML[Messung_Nr].m_Fehler_Zeilendichten =
				(2.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;

			Messung_Nr = u + (Anzahl_Hoehen - 1) * Anzahl_Linien + k;
			AML[Messung_Nr].m_Zeilendichte =
				(2.0 * AML_OLD[Messung_Nr].m_Zeilendichte
				 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Zeilendichte
				 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte
				 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte) / 5.0;
			//cerr<<"Kante 4\n";
			//cerr<<" Beitrag aus: "<<Messung_Nr<<" "<<Messung_Nr -Anzahl_Linien
			//  <<" "<<Messung_Nr-Anzahl_Linien*Anzahl_Hoehen<<" "
			//  <<Messung_Nr+Anzahl_Linien*Anzahl_Hoehen<<"\n";
			AML[Messung_Nr].m_Fehler_Zeilendichten =
				(2.0 * AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten
				 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten) / 5.0;
		} //ende for u

		// Innere Elemente
		for (int u = 1; u < (Anzahl_Hoehen - 1); u++) {
			for (int w = Anzahl_Hoehen * Anzahl_Linien;
					w < (Anzahl_Messungen - Anzahl_Hoehen * Anzahl_Linien);
					w += Anzahl_Hoehen * Anzahl_Linien) {
				Messung_Nr = w + u * Anzahl_Linien + k;
				AML[Messung_Nr].m_Zeilendichte =
					(AML_OLD[Messung_Nr].m_Zeilendichte
					 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Zeilendichte
					 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Zeilendichte
					 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Zeilendichte
					 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen ].m_Zeilendichte) / 5.0;
				//cerr<<"u= "<< u<<" w= "<<w<<"\n";
				//cerr<<" Beitrag aus: "<<Messung_Nr<<" "
				// <<Messung_Nr +Anzahl_Linien<<" "<<Messung_Nr -Anzahl_Linien
				// <<" "<<Messung_Nr-Anzahl_Linien*Anzahl_Hoehen<<" "
				// <<Messung_Nr+Anzahl_Linien*Anzahl_Hoehen<<"\n";
				AML[Messung_Nr].m_Fehler_Zeilendichten =
					(AML_OLD[Messung_Nr].m_Fehler_Zeilendichten
					 + AML_OLD[Messung_Nr + Anzahl_Linien].m_Fehler_Zeilendichten
					 + AML_OLD[Messung_Nr - Anzahl_Linien].m_Fehler_Zeilendichten
					 + AML_OLD[Messung_Nr + Anzahl_Linien * Anzahl_Hoehen].m_Fehler_Zeilendichten
					 + AML_OLD[Messung_Nr - Anzahl_Linien * Anzahl_Hoehen ].m_Fehler_Zeilendichten)
					/ 5.0;
			}//ende for w
		}//Ende for u
	}// ende for k

	return 0;
}
