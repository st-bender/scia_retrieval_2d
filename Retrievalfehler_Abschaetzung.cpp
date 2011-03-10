/*
 * Retrievalfehler_Abschaetzung.cpp
 *
 *  Created on: 16.09.2010
 *      Author: martin
 */

#include "Retrievalfehler_Abschaetzung.h"
#include "Konfiguration.h"
#include "MPL_Matrix.h"

using namespace std;

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO);
	void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
				double *ALPHA, double *A, int *LDA, double *B, int *LDB, double *BETA, double *C, int *LDC);
}

int Retrievalfehler_Abschaetzung(MPL_Matrix &S_x,
								 MPL_Matrix &Averaging_Kernel_Matrix,
								 const MPL_Matrix &S_apriori,
								 const MPL_Matrix &S_y,
								 MPL_Matrix S_Breite,
								 MPL_Matrix S_Hoehe,
								 MPL_Matrix S_letzte_Hoehe,
								 const double &Lambda_Breite,
								 const double &Lambda_Hoehe,
								 MPL_Matrix AMF,
								 const Konfiguration &Konf)
{
	//TODO Auch hier kann man das Gleichungssystem mit ATLAS/LAPACK FUNKTIONEN LÖSEN



	// Die Formeln für die Matrizzen findet man in Marcos Arbeit
	MPL_Matrix S_x_invers;
	MPL_Matrix AMF_trans;
	AMF_trans = AMF.transponiert();
	MPL_Matrix S_Hoehe_trans;
	S_Hoehe_trans = S_Hoehe.transponiert();
	MPL_Matrix S_Breite_trans;
	S_Breite_trans = S_Breite.transponiert();
	MPL_Matrix S_letzte_Hoehe_trans;
	S_letzte_Hoehe_trans = S_Breite.transponiert();


	S_x = AMF_trans * (S_y * AMF) // hier noch invers, also noch invertieren
		  + S_apriori
		  + Lambda_Hoehe * (S_Hoehe_trans * S_Hoehe)
		  + Lambda_Breite * (S_Breite_trans * S_Breite)
		  + S_letzte_Hoehe_trans * S_letzte_Hoehe;
	Matrix_Invertieren(S_x);
	//  cout<<S_x_invers.m_Zeilenzahl<<"\t"<<S_x_invers.m_Spaltenzahl<<"\n";
	//  Die Matrix sollte quadratisch sein
	//  INVERSION
	//   int dim=S_x_invers.m_Zeilenzahl;
	//   MPL_Matrix Inversionsgleichungssystem(dim,2*dim);
	//   Inversionsgleichungssystem.Null_Initialisierung();
	//   S_x_Matrix als linke Blockmatrizze
	//   for(int zeile=0;zeile<dim;zeile++)
	//  {
	//    Inversionsgleichungssystem(zeile,dim+zeile)=1;
	//    for(int spalte=0;spalte<dim;spalte++)
	//    {
	//        Inversionsgleichungssystem(zeile,spalte)=S_x_invers(zeile,spalte);
	//    }//for spalte
	//  }// for zeile
	//  // Diagonalisieren
	//    int IERR;
	//  IERR=Inversionsgleichungssystem.Gausselimination_mit_Teilpivotisierung_ohne_Skalenfaktor();
	//  if(IERR==1)
	//  {
	//      cout<<"Inversion zur Bestimmung der Fehlermatrix gescheitert \n";
	//    return 1;
	//  }
	//Inversionsgleichungssystem.in_Datei_speichern("/tmp/mlangowski/0/Inversionsgleichungssystem.txt");
	// Rechte Block Matrix übergeben
	//for(int zeile=0;zeile<dim;zeile++)
	//{
	//    for(int spalte=0;spalte<dim;spalte++)
	//    {
	//        S_x(zeile,spalte)=Inversionsgleichungssystem(zeile,dim+spalte);
	//    }//ende spalte
	//}//ende zeile
	//Nun noch die Averaging Kernel Matrix bestimmen
	Averaging_Kernel_Matrix = S_x * (AMF_trans * (S_y * AMF));
	return 0;
}
///////////////////////////////////////////////////////
// Funktionsstart Matrix_Invertieren
///////////////////////////////////////////////////////
void Matrix_Invertieren(MPL_Matrix &M)
{
	// Die Funktion wandelt die Matrix M in ihre eigene Inverse um
	// Dafür muss das gehen (M quadratisch und nicht singulär)
	int N = M.m_Zeilenzahl; // M ist quadratisch !!!!
	MPL_Matrix Inverse(N, N);                                   // Inverse als Einheitsmatrix der RHS initialisieren
	int NN = N * N; // spare die Multiplikation
	for (int i = 0; i < NN; i++) {
		Inverse.m_Elemente[i] = 0;
	}
	for (int i = 0; i < N; i++) {
		Inverse.m_Elemente[i + N * i] = 1;
	}
	int *IPIV;
	IPIV = new int[N];              // Lösung des Gleichungssystems M*Inverse=RHS vorbereiten
	for (int i = 0; i < N; i++)  {
		IPIV[i] = 0;
	}
	int NRHS = N;
	int LDA = N;
	int LDB = N;
	int INFO;
	// M in Fortran Matrix umwandeln->M_transponieren
	M = M.transponiert();
	dgesv_(&N, &NRHS, M.m_Elemente, &LDA, IPIV, Inverse.m_Elemente, &LDB, &INFO);  //Solveraufruf
	// Die Inverse von Fortran in C++ -> Inverse transponieren
	M = Inverse.transponiert();                                  // Ergebnis in M deponieren
	delete[] IPIV;
}
///////////////////////////////////////////////////////
// ENDE Matrix_Invertieren
///////////////////////////////////////////////////////
