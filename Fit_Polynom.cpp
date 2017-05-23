/*
 * Fit_Polynom.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 29.10.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include"Fit_Polynom.h"
#include<iostream>

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
	void dgetrs_(char *, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
}

////////////////////////////////////////////////////////////////////////////////
// Funktionsstart Fit_Polynom
////////////////////////////////////////////////////////////////////////////////
void Fit_Polynom(double *x, double *y, int Startindex, int Endindex, double x0,
		int Polynomgrad, double *Parametervektor)
{
	// Fit eines Polynoms a0+a1x+a2x^2.......+anx^n    mit n=Polynomgrad
	// Dabei wird zunächst ein lineares Gleichungssystem aufgestellt
	// und dann gelöst  (Beispiel in Matlab für Polynom 4ten Grades,
	// siehe unten....und Anleitung)


	double *LHS_Parameter_Zeile_eins;
	double *LHS_Parameter_letzte_Spalte;
	double *RHS;
	double *LHS;
	double *LHS_Parameter_Zeile_eins_Inkrement;
	double *LHS_Parameter_letzte_Spalte_Inkrement;
	double *RHS_Inkrement;
	//Speicher reservieren
	LHS_Parameter_Zeile_eins = new double[Polynomgrad + 1];
	LHS_Parameter_letzte_Spalte = new double[Polynomgrad + 1];
	RHS = new double[Polynomgrad + 1];
	LHS = new double[(Polynomgrad + 1) * (Polynomgrad + 1)];
	LHS_Parameter_Zeile_eins_Inkrement = new double[Polynomgrad + 1];
	LHS_Parameter_letzte_Spalte_Inkrement = new double[Polynomgrad + 1];
	RHS_Inkrement = new double[Polynomgrad + 1];
	// Alle Elemente mit 0 starten lassen
	for (int i = 0; i < Polynomgrad + 1; i++) {
		LHS_Parameter_Zeile_eins[i] = 0;
		LHS_Parameter_letzte_Spalte[i] = 0;
		RHS[i] = 0;
	}
	for (int i = 0; i < (Polynomgrad + 1) * (Polynomgrad + 1); i++) {
		LHS[i] = 0;
	}
	// Parameter aufbauen
	for (int Index = Startindex; Index <= Endindex; Index++) {
		double h = x[Index] - x0;
		//cout<<"x[Index]: "<<x[Index]<<" \n";
		//cout<<"x0: "<<x0<<" \n";
		//cout<<"h: "<<h<<" \n";
		LHS_Parameter_Zeile_eins_Inkrement[0] = 1;
		RHS_Inkrement[0]                           = y[Index];
		for (int Par_Zeile = 1; Par_Zeile < Polynomgrad + 1; Par_Zeile++) {
			LHS_Parameter_Zeile_eins_Inkrement[Par_Zeile]
				= LHS_Parameter_Zeile_eins_Inkrement[Par_Zeile - 1] * h;
			RHS_Inkrement[Par_Zeile] = RHS_Inkrement[Par_Zeile - 1] * h;
		}
		LHS_Parameter_letzte_Spalte_Inkrement[0]
			= LHS_Parameter_Zeile_eins_Inkrement[Polynomgrad];
		for (int Par_Spalte = 1; Par_Spalte < Polynomgrad + 1; Par_Spalte++) {
			LHS_Parameter_letzte_Spalte_Inkrement[Par_Spalte]
				= LHS_Parameter_letzte_Spalte_Inkrement[Par_Spalte - 1] * h;
		}
		// Inkrementationen durchführen
		for (int i = 0; i < Polynomgrad + 1; i++) {
			LHS_Parameter_Zeile_eins[i] = LHS_Parameter_Zeile_eins[i]
				+ LHS_Parameter_Zeile_eins_Inkrement[i];
			LHS_Parameter_letzte_Spalte[i] = LHS_Parameter_letzte_Spalte[i]
				+ LHS_Parameter_letzte_Spalte_Inkrement[i];
			RHS[i] = RHS[i] + RHS_Inkrement[i];
		}
	}
	// Parameter sind bekannt. Nun heißt es Matrix zusammenbauen
	// Der Solver ist ein Fortransolver und braucht die Matrix Spaltenweise
	// erste Zeile der Matrix Eintragen
	int Zeilennummer, Spaltennummer;
	Spaltennummer = 0;
	//cout<<"erste Zeile\n";
	for (Zeilennummer = 0; Zeilennummer < Polynomgrad + 1; Zeilennummer++) {
		LHS[Spaltennummer + Zeilennummer * (Polynomgrad + 1)]
			= LHS_Parameter_Zeile_eins[Zeilennummer];
		//cout<<"Element_LHS: "<<Spaltennummer+Zeilennummer*(Polynomgrad+1)<<"\n";
	}
	Zeilennummer = Polynomgrad;
	//cout<<"erste Spalte\n";
	for (Spaltennummer = 0; Spaltennummer < Polynomgrad + 1; Spaltennummer++) {
		LHS[Spaltennummer + Zeilennummer * (Polynomgrad + 1)]
			= LHS_Parameter_letzte_Spalte[Spaltennummer];
		//cout<<"Element_LHS: "<<Spaltennummer+Zeilennummer*(Polynomgrad+1)<<"\n";
	}
	// die erste Spalte und die letzte Zeile der Matrix sind nun besetzt
	// ab der 2ten Spalte die Spalte bis zur vorletzten Zeile durchgehen und
	// das Element mit einer Spalte weniger und einer Zeile mehr zuordnen (also
	// das links unten)
	for (Spaltennummer = 1; Spaltennummer < Polynomgrad + 1; Spaltennummer++) {
		for (Zeilennummer = 0; Zeilennummer < Polynomgrad; Zeilennummer++) {
			LHS[Spaltennummer + Zeilennummer * (Polynomgrad + 1)]
				= LHS[(Spaltennummer - 1) + (Zeilennummer + 1) * (Polynomgrad + 1)];
			//cout<<"Element_LHS: "<<Spaltennummer+Zeilennummer*(Polynomgrad+1)<<"\n";
		}
	}
	//RHS und LHS sind fertig
	//LAPACK Solver Routine vorbereiten
	int N = Polynomgrad + 1; //<---------- Feldgröße Speed propto N^3,
							 //LHS ist quadratisch, N ist Anz. der Gitterpunkte
	int *IPIV;  //array mit der Pivotisierungsmatrix sollte so groß wie N sein, alle Elemente 0
	IPIV = new int[N];
	int NRHS = 1; //Spalten von RHS 1 nehmen, um keine c/Fortran Verwirrungen zu provozieren
	int LDA = N;
	int LDB = N;
	int INFO;
	//Aufruf der LU-Zerlegung
	dgesv_(&N, &NRHS, LHS, &LDA, IPIV, RHS, &LDB, &INFO);
	// Lösungen liegen nun in RHS vor
	// Speicher für Parametervektor,
	// wurde vor Funktionsaufruf (Fit_Polynom) reserviert
	for (int i = 0; i < Polynomgrad + 1; i++) {
		Parametervektor[i] = RHS[i];
	}
	//Speicher freigeben
	delete[] LHS;
	delete[] LHS_Parameter_Zeile_eins;
	delete[] LHS_Parameter_letzte_Spalte;
	delete[] RHS;
	delete[] IPIV;
	delete[] LHS_Parameter_letzte_Spalte_Inkrement;
	delete[] RHS_Inkrement;
	delete[] LHS_Parameter_Zeile_eins_Inkrement;
	// Das Fitten eines Polynoms an einen Datensatz von x,y Punkten
	// Geschieht über das minimieren der Summe der Quadrate der
	// Abweichungen(Residuen) der Fitfunktion von den Messpunkten
	// Sum( (f-y)^2)=min
	// Leitet man ab erhält man für jeden linearen Parameter....
	// (das Wort linear dick unterstreichen)
	// Sum (f'f-f'y)=0  mit f'=df/da_i Ableitung nach linearem Parameter
	// in f' steckt a_i nichtmehr drin...
	// sodass man f'y auf die andere Seite zieht,
	// alle Gleichungen f'y ergeben dann die Rechte Seite von LHS x=RHS;
	// jede Summe f'f lässt sich als linearkombination
	// der parameter a_i darstellen, sodass man so LHS x
	// mit LHS Matrix und x Vektor der Parameter benutzen kann
	// für ein Polynom a0+a1 x+a2x^2... ergibt df/da_i gleich x^i
	//
	// konkretes Beispiel //
	// Zu illustrativen Zwecken das Beispiel für Polynomgrad=4 in Matlab
	//a0=0; b0=0; c0=0; d0=0; e0=0;
	//                                     e1=0;
	//                                     e2=0;
	//                                     e3=0;
	//                                     e4=0;
	// a0 b0 c0 d0 e0              //Beachte die Diagonalen r.o.-l.u.
	// b0 c0 d0 e0 e1
	// c0 d0 e0 e1 e2
	// d0 e0 e1 e2 e3
	// e0 e1 e2 e3 e4
	//
	//f0=0; f1=0; f2=0; f3=0;f4=0;
	//
	//for i=1:lang
	//    h=x(i)-mu;
	//    a0=a0+1;      // bestimmen der parameter...
	//                  // bei jedem weiteren *h...das schreit nach Schleife
	//    b0=b0+h;      //LHS
	//    c0=c0+h*h;
	//    d0=d0+h*h*h;
	//    e0=e0+h*h*h*h;
	//    e1=e1+h*h*h*h*h;
	//    e2=e2+h*h*h*h*h*h;
	//    e3=e3+h*h*h*h*h*h*h;
	//    e4=e4+h*h*h*h*h*h*h*h;
	//
	//    f0=f0+y(i);                         //RHS
	//    f1=f1+y(i)*h;
	//    f2=f2+y(i)*h*h;
	//    f3=f3+y(i)*h*h*h;
	//    f4=f4+y(i)*h*h*h*h;
	//end
	//dim=5;
	//LHS=zeros(dim,dim);
	// RHS=zeros(dim,1);
	// LHS(1,1)=a0; LHS(1,2)=b0; LHS(1,3)=c0; LHS(1,4)=d0; LHS(1,5)=e0;
	// LHS(2,1)=b0; LHS(2,2)=c0; LHS(2,3)=d0; LHS(2,4)=e0; LHS(2,5)=e1;
	// LHS(3,1)=c0; LHS(3,2)=d0; LHS(3,3)=e0; LHS(3,4)=e1; LHS(3,5)=e2;
	// LHS(4,1)=d0; LHS(4,2)=e0; LHS(4,3)=e1; LHS(4,4)=e2; LHS(4,5)=e3;
	// LHS(5,1)=e0; LHS(5,2)=e1; LHS(5,3)=e2; LHS(5,4)=e3; LHS(5,5)=e4;
	// RHS(1)=f0; RHS(2)=f1; RHS(3)=f2; RHS(4)=f3; RHS(5)=f4;
	// % Gleichung lösen
	// A=LHS\RHS                   // Lösung der Gleichung
	// % Lösungen übergeben
	// a0=A(1);    // a1=A(2);    // a2=A(3);    // a3=A(4);    // a4=A(5);
}
////////////////////////////////////////////////////////////////////////////////
// Ende Fit_Polynom
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Funktionsstart x_zu_Minimum_von_y_in_Intervall
////////////////////////////////////////////////////////////////////////////////
void x_zu_Minimum_von_y_in_Intervall(double *x, double *y, int Startindex,
									int Endindex, double &x_min, double &y_min,
									int &Indexmin)
{
	int Min_index = Startindex;
	for (int i = Startindex; i <= Endindex; i++) {
		if (y[i] < y[Min_index])
			Min_index = i;
	}
	x_min = x[Min_index];
	y_min = y[Min_index];
	Indexmin = Min_index;
}
////////////////////////////////////////////////////////////////////////////////
// Ende x_zu_Minimum_von_y_in_Intervall
////////////////////////////////////////////////////////////////////////////////
