/*
 * Retrievalfehler_Abschaetzung.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 16.09.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include "Retrievalfehler_Abschaetzung.h"
#include "Konfiguration.h"
#include "MPL_Matrix.h"

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B,
				int *LDB, int *INFO);
}

void Retrievalfehler_Abschaetzung(MPL_Matrix &S_x,
								 MPL_Matrix &S_x_meas,
								 MPL_Matrix &Averaging_Kernel_Matrix,
								 const MPL_Matrix &S_apriori,
								 const MPL_Matrix &S_y,
								 MPL_Matrix &S_Breite,
								 MPL_Matrix &S_Hoehe,
								 const double &Lambda_Breite,
								 const double &Lambda_Hoehe,
								 MPL_Matrix &AMF,
								 const Konfiguration &Konf)
{
	//TODO Auch hier kann man das Gleichungssystem
	//mit ATLAS/LAPACK FUNKTIONEN LÖSEN

	// Die Formeln für die Matrizen findet man in Marcos Arbeit
	MPL_Matrix AMF_trans = AMF.transponiert();
	MPL_Matrix S_Hoehe_trans = S_Hoehe.transponiert();
	MPL_Matrix S_Breite_trans = S_Breite.transponiert();

	// S_y here is equal to S_y^-1 in ususal retrieval equations,
	// as is S_apriori (~ S_a^-1)
	S_x = (AMF_trans * (S_y * AMF)); // hier noch invers, also noch invertieren
	S_x += S_apriori;
	Matrix_Invertieren(S_x);
	//Nun noch die Averaging Kernel Matrix bestimmen
	MPL_Matrix D = S_x * AMF_trans * S_y;
	Averaging_Kernel_Matrix = D * AMF;
	MPL_Matrix S_y_inv = S_y;
	Matrix_Invertieren(S_y_inv);
	S_x_meas = D * S_y_inv * D.transponiert();
}
///////////////////////////////////////////////////////
// Funktionsstart Matrix_Invertieren
///////////////////////////////////////////////////////
void Matrix_Invertieren(MPL_Matrix &M)
{
	// Die Funktion wandelt die Matrix M in ihre eigene Inverse um
	// Dafür muss das gehen (M quadratisch und nicht singulär)
	int N = M.m_Zeilenzahl; // M ist quadratisch !!!!
	MPL_Matrix dummy(N, N);
	MPL_Matrix Inverse = dummy.unity(); // Inverse als Einheitsmatrix der RHS initialisieren
	int *IPIV = new int[N]; // Lösung des Gleichungssystems M*Inverse=RHS vorbereiten
	int NRHS = N;
	int LDA = N;
	int LDB = N;
	int INFO;
	//Solveraufruf
	dgesv_(&N, &NRHS, M.m_Elemente, &LDA, IPIV, Inverse.m_Elemente, &LDB, &INFO);
	// Die Inverse von Fortran in C++ -> Inverse transponieren
	M = Inverse;  // Ergebnis in M deponieren
	delete[] IPIV;
}
///////////////////////////////////////////////////////
// ENDE Matrix_Invertieren
///////////////////////////////////////////////////////
