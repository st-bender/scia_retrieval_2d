/*
 * Matrix.h
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 28.05.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
/*
 * Changelog:
 * ==========
 *
 * 6.10.2010
 *  + und - hatten Matrixdimensionen vertauscht.....autsch
 *
 * 5.10.2010
 * Nutzen von ATLAS für Matrixmultiplikationen implementiert
 *
 * 29.09.2010
 * in datei speichern  schreiben in outfile nicht in cout
 *
 *  Last update: 23.09.2010
 *  *= und += waren falsch implementiert
 *  (wieso das erst jetzt auffällt...hatte odch alles getestet)
 *
 *  16.09.2010
 *  gausselimination mit Pivotisierung:
 *    das Zeilentauschen der RHS am Ende hinzugefügt
 *  gausselimination mit Pivotisierung
 *  ohne Pivotisierung ist glaub ich falsch und
 *    mit pivotisierung also erstmal auskommentieren
 *
 *  Created on: 28.05.2010
 *      Author: martin
 *      Matrixdefinitionen gibt es viele, MPL steht dann für meinen Namen,
 *      das führt dann hoffentlich nicht zu doppelt definierten Matrixtypen,
 *      was Mit Matrix alleine ziemlich sicher passieren würde.
 *
 *     Es ist durchaus sinnvoll Vektoren als Zeilen und Spaltenvektoren also
 *     als einzeilige oder einspaltige Matrix zu nutzen, weil dadurch für beide
 *     die selbe Algebra verwendet werden kann
 */
/*
 *  Zu den Rückgabewerten MPLMatrix,
 *    also Rückgabe als Variable ist immer möglich, aber etwas langsamer
 *  Rückgabe als Referenz MPLMatrix&,
 *    oder als MPLMatrix* brauchen die Rückgabevariable vorher
 *
 *  //Die Elemente dieser Matrix sind Zeilenweise angeordnet
 *    d.h. Elementzahl(i,j)=i+j*Spaltenzahl
 */
#ifndef MPLMATRIX_HH_
#define MPLMATRIX_HH_

#include<iostream>
#include <cstdio>
#include<fstream>
#include <algorithm>
#include <vector>
#include "gzstream.h"
#include "netcdf.h"

#ifdef HAVE_HDF5
#include "hdf5.h"
#endif /* HAVE_HDF5 */

extern "C" {
	void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
				double *ALPHA, const double *A, int *LDA, const double *B, int *LDB,
				double *BETA, double *C, int *LDC);
}
//////////////////////////////////////////////////////
// START KLASSENDEKLARATION
//////////////////////////////////////////////////////
class MPL_Matrix
{
public:
	//Konstruktoren //////////////////////////////////
	MPL_Matrix() : transposed(false), m_Zeilenzahl(0), m_Spaltenzahl(0),
		m_Elementanzahl(0), m_Elemente(0) {}
	MPL_Matrix(int Zeilenzahl, int Spaltenzahl, double value = 0.0) :
		transposed(false),
		m_Zeilenzahl(Zeilenzahl), m_Spaltenzahl(Spaltenzahl),
		m_Elementanzahl(Zeilenzahl * Spaltenzahl)
		{
			m_Elemente = new double[m_Elementanzahl];
			std::fill_n(m_Elemente, m_Elementanzahl, value);
		}
	MPL_Matrix(const MPL_Matrix &rhs);
	// Hier ist nochmehr denkbar z.b. einheitzsmatrix 0 matrix usw
	/////////////////////////////////////////////////////////
	//Destruktor///
	~MPL_Matrix();

	////////////////////
	//Überladene Operatoren/////////////////////
	// Diese werden inline definiert. damit das quasi in dieser Datei geschieht
	// wird am ende Eine .inl Funktion included
	// Zuweisungen
	MPL_Matrix &operator = (const MPL_Matrix &rhs);   // Zuweisung
	MPL_Matrix &operator += (const MPL_Matrix &rhs);  // Matrixaddition
	MPL_Matrix &operator -= (const MPL_Matrix &rhs);  // Matrixsubtraktion
	MPL_Matrix &operator *= (double rhs);  // Skalare Multiplikation
	MPL_Matrix &operator /= (double rhs);  // Skalare Division
	// Beachte...es gibt keine Division, da nicht jede Matrix eine Inverse hat

	//() Überladung -> Direkter Zugriff auf Das Elemente Array
	double &operator()(int Zeile, int Spalte) const;
	double &operator()(int Elementindex) const;

	// binary operators
	MPL_Matrix operator * (const MPL_Matrix &rhs) const;   //matmul
	MPL_Matrix operator + (const MPL_Matrix &rhs) const;   //matadd
	MPL_Matrix operator - (const MPL_Matrix &rhs) const;   //matsub
	MPL_Matrix operator * (const double &rhs) const;       //skalare Mult
	MPL_Matrix operator / (const double &rhs) const;       //skalare Div

	friend MPL_Matrix operator * (const double &lhs, const MPL_Matrix &rhs);

	bool operator == (const MPL_Matrix &rhs) const;
	////////////////////////////////////////////////////////////////

	// Methoden
	MPL_Matrix get_Zeile(int Zeilennummer); // gibt eine Zeile als Spaltenvektor aus
	MPL_Matrix get_Spalte(int Spaltennummer); // gibt eine Spalte als Spaltenvektor aus

	void transpose();
	MPL_Matrix transponiert() const; //transponierte Matrix
	MPL_Matrix transponiert_full() const; //transponierte Matrix
	MPL_Matrix transponiert_full2() const; //transponierte Matrix
//    MPLMatrix  invertiert();
//    //inverse Matrix, falls existent...
//      existiert nur bei quadratischen, nicht singulären Matrizen
	// transponieren, Inverse Matrix / Gauss, LU, Cholesky usw SVD
	void Zeile_Tauschen(int Zeile_a, int Zeile_b);
	void Zeile_Multiplizieren(int Zeile, double Faktor);
	void Vielfaches_einer_Zeile_addieren(int Summenzeile, int Additionszeile, double Faktor);
	MPL_Matrix row_diff();
	MPL_Matrix unity() const;
	double trace();
	//int simple_Gaussdiagonalisierung(); siehe ganz oben
	int Gausselimination_mit_Teilpivotisierung_ohne_Skalenfaktor();
	void in_Datei_speichern(std::string Dateiname, double precision = 0) const;
	int save_to_netcdf(std::string Dateiname, bool pack = false) const;
	int save_to_hdf5(std::string Dateiname, bool pack = false) const;

	//Membervariablen
	bool transposed;
	int m_Zeilenzahl;
	int m_Spaltenzahl;
	int m_Elementanzahl;
	double *m_Elemente;
};
//////////////////////////////////////////////////////
// ENDE KLASSENDEKLARATION
//////////////////////////////////////////////////////
//Deklaration fertig...nun Operatoren als inline Funktionen überladen
//Konstruktoren //////////////////////////////////
/////////////////////////////////////////////////////////
// Methodenstart MPL_Matrix Konstruktor
/////////////////////////////////////////////////////////
inline MPL_Matrix::MPL_Matrix(const MPL_Matrix &rhs)
{
	m_Elemente = 0;
	*this = rhs;
}
/////////////////////////////////////////////////////////
// Methodenstart Destruktor
/////////////////////////////////////////////////////////
inline MPL_Matrix::~MPL_Matrix()
{
	if (m_Elemente != 0) {
		delete[] m_Elemente;
		m_Elemente = 0;
	}
}
//Überladene Operatoren/////////////////////
//Diese werden inline definiert.
// damit das quasi in dieser Datei geschieht wird am ende
// Eine .inl Funktion included
// Zuweisungen
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// Zuweisung
inline MPL_Matrix &MPL_Matrix::operator = (const MPL_Matrix &rhs)
{
	if (this == &rhs)
		return *this;
	this->transposed = rhs.transposed;
	this->m_Zeilenzahl = rhs.m_Zeilenzahl;
	this->m_Spaltenzahl = rhs.m_Spaltenzahl;
	this->m_Elementanzahl = rhs.m_Elementanzahl;
	if (m_Elemente != 0) {
		delete[] m_Elemente;
		m_Elemente = 0;
	}
	m_Elemente = new (std::nothrow) double[m_Elementanzahl];
	if (!m_Elemente) {
		std::cout << "out of memory: cannot allocate "
			<< m_Elementanzahl << " doubles." << std::endl;
		exit(1);
	}
	std::copy_n(rhs.m_Elemente, m_Elementanzahl, m_Elemente);
	return *this;
}
// ende operator =

/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// Matrixaddition
inline MPL_Matrix &MPL_Matrix::operator +=(const MPL_Matrix &rhs)
{
	//Zeilen und Spaltenzahl muss gleich sein
	if ((m_Spaltenzahl != rhs.m_Spaltenzahl)
			&& (m_Zeilenzahl != rhs.m_Zeilenzahl)) {
		std::cerr << "Addition nicht möglich, ungleiche Matrixdimensionen"
				  << std::endl;
		return *this; //einfach garnix gemacht
	}
	// Wenn spaltenzahl und Zeilenzahl gleich,
	// dann ist auch die reihenfolge der elemente gleich
	std::transform(m_Elemente, m_Elemente + m_Elementanzahl, rhs.m_Elemente,
			m_Elemente, std::plus<double>());
	return *this;
} //ende +=
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// Matrixsubtraktion
inline MPL_Matrix &MPL_Matrix::operator -= (const MPL_Matrix &rhs)
{
	//Zeilen und Spaltenzahl muss gleich sein
	if ((m_Spaltenzahl != rhs.m_Spaltenzahl)
			&& (m_Zeilenzahl != rhs.m_Zeilenzahl)) {
		std::cerr << "Addition nicht möglich, ungleiche Matrixdimensionen"
				  << std::endl;
		return *this; //einfach garnix gemacht
	}
	// Wenn spaltenzahl und Zeilenzahl gleich,
	// dann ist auch die reihenfolge der elemente gleich
	std::transform(m_Elemente, m_Elemente + m_Elementanzahl, rhs.m_Elemente,
			m_Elemente, std::minus<double>());
	return *this;
} // ende -=
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// Skalare Multiplikation
inline MPL_Matrix &MPL_Matrix::operator *= (double rhs)
{
	std::transform(m_Elemente, m_Elemente + m_Elementanzahl, m_Elemente,
			std::bind2nd(std::multiplies<double>(), rhs));
	return *this;
}// ende *=
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// Skalare Division
inline MPL_Matrix &MPL_Matrix::operator /= (double rhs)
{
	if (rhs == 0) {
		std::cerr << "Division durch 0 wird nicht durchgeführt!" << std::endl;
		return *this;
	}
	std::transform(m_Elemente, m_Elemente + m_Elementanzahl, m_Elemente,
			std::bind2nd(std::divides<double>(), rhs));
	return *this;
}// ende /=

//() Überladung -> Direkter Zugriff auf Das Elemente Array
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline double &MPL_Matrix::operator()(int Zeile, int Spalte) const
{
	//A(1,2)=b;
	int idx = transposed ? Spalte * m_Zeilenzahl + Zeile
						 : Zeile * m_Spaltenzahl + Spalte;
	if ((uint)idx < (uint)m_Elementanzahl)
		return m_Elemente[idx];
	else {
		std::cerr << "Achtung!!! Zugriff auf Elemente ausserhalb der Matrix"
				  << std::endl;
		return m_Elemente[0];//auch schlecht, aber wenigstens nicht ausserhalb
	}
}
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline double &MPL_Matrix::operator()(int Elementindex) const
{
	if ((uint)Elementindex < (uint)m_Elementanzahl)
		return this->m_Elemente[Elementindex];
	else {
		std::cerr << "Achtung!!! Zugriff auf Elemente ausserhalb der Matrix"
				  << std::endl;
		return m_Elemente[0];//auch schlecht, aber wenigstens nicht ausserhalb
	}
}

// binary operators
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline MPL_Matrix MPL_Matrix::operator * (const MPL_Matrix &rhs) const
{
	//Zunächst prüfen, ob Multiplikation möglich ist
	if (this->m_Spaltenzahl != rhs.m_Zeilenzahl) {
		std::cerr << "*= Wrong Matrix Multiplication A*B, "
			 << "coloums number of A != rows number of B" << std::endl;
		std::cerr << "returning nonsense!!!!" << std::endl;
		return *this;
	}
	// gemm vorbereiten  (das ist ein Routine aus ATLAS)
	// Set memory order according to the transposed flag:
	// true = Fortran order, no transpose in dgemm_(),
	// false = C order, use transpose in dgemm_())
	char TRANSA = transposed ? 'n' : 't';
	char TRANSB = rhs.transposed ? 'n' : 't';
	int M = this->m_Zeilenzahl;
	int N = rhs.m_Spaltenzahl;
	int K = this->m_Spaltenzahl;
	double ALPHA = 1.0;
	int LDA = transposed ? M : K;
	int LDB = rhs.transposed ? K : N;
	double BETA = 0.0;
	int LDC = M;
	MPL_Matrix Produkt(M, N);
	// Matrixmultiplikation durchführen
	dgemm_(&TRANSA, &TRANSB, &M, &N, &K,
		   &ALPHA, m_Elemente, &LDA, rhs.m_Elemente, &LDB, &BETA,
		   Produkt.m_Elemente, &LDC);
	// We need to keep track of the element ordering in memory.
	// dgemm_() yields Fortran order which is transposed to C.
	// All subsequent multiplications will take care of that.
	Produkt.transposed = true;
	return Produkt;
}
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline MPL_Matrix MPL_Matrix::operator +(const MPL_Matrix &rhs) const
{

	//Zeilen und Spaltenzahl muss gleich sein
	if ((m_Spaltenzahl != rhs.m_Spaltenzahl)
			&& (m_Zeilenzahl != rhs.m_Zeilenzahl)) {
		std::cerr << "Addition nicht möglich, ungleiche Matrixdimensionen"
				  << std::endl;
		return *this; //einfach garnix gemacht
	}
	MPL_Matrix Summe(rhs.m_Zeilenzahl, rhs.m_Spaltenzahl);
	// Wenn spaltenzahl und Zeilenzahl gleich,
	// dann ist auch die reihenfolge der elemente gleich
	std::transform(m_Elemente, m_Elemente + m_Elementanzahl, rhs.m_Elemente,
			Summe.m_Elemente, std::plus<double>());
	return Summe;
}
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline MPL_Matrix MPL_Matrix::operator - (const MPL_Matrix &rhs) const
{
	//Zeilen und Spaltenzahl muss gleich sein
	if ((m_Spaltenzahl != rhs.m_Spaltenzahl)
			&& (m_Zeilenzahl != rhs.m_Zeilenzahl)) {
		std::cerr << "Subtraktion nicht möglich, ungleiche Matrixdimensionen"
				  << std::endl;
		return *this; //einfach garnix gemacht
	}
	MPL_Matrix Differenz(rhs.m_Zeilenzahl, rhs.m_Spaltenzahl);
	// Wenn spaltenzahl und Zeilenzahl gleich,
	// dann ist auch die reihenfolge der elemente gleich
	std::transform(m_Elemente, m_Elemente + m_Elementanzahl, rhs.m_Elemente,
			Differenz.m_Elemente, std::minus<double>());
	return Differenz;
}//ende -
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline MPL_Matrix MPL_Matrix::operator * (const double &rhs) const
{
	MPL_Matrix Skalare_Multiplikation;
	Skalare_Multiplikation = *this;
	std::transform(Skalare_Multiplikation.m_Elemente,
			Skalare_Multiplikation.m_Elemente + Skalare_Multiplikation.m_Elementanzahl,
			Skalare_Multiplikation.m_Elemente,
			std::bind2nd(std::multiplies<double>(), rhs));
	return Skalare_Multiplikation;
}// ende skalare Mult
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline MPL_Matrix MPL_Matrix::operator / (const double &rhs) const
{
	if (rhs == 0) {
		// ==0 ist bei double ziemlich undefiniert...
		// evtl nochmal grenzen einbaun
		std::cerr << "Division durch 0 wird nivht durchgeführt" << std::endl;
		return *this;
	}
	MPL_Matrix Skalare_Division;
	Skalare_Division = *this;
	std::transform(Skalare_Division.m_Elemente,
			Skalare_Division.m_Elemente + Skalare_Division.m_Elementanzahl,
			Skalare_Division.m_Elemente,
			std::bind2nd(std::divides<double>(), rhs));
	return Skalare_Division;
}//Ende skalare Division
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline MPL_Matrix operator * (const double &lhs, const MPL_Matrix &rhs)
{
	MPL_Matrix skalares_Produkt;
	skalares_Produkt = rhs;
	skalares_Produkt *= lhs; //nutze andere Überladene Operatoren
	return skalares_Produkt;
}// Ende *
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline bool MPL_Matrix::operator == (const MPL_Matrix &rhs) const
{
	//2 Matrizen sind gleich, wenn sie die gleichen Dimensionen haben,
	//und alle ihre Elemente gleich sind
	// erstmal Dimensionen prüfen
	if ((this->m_Spaltenzahl != rhs.m_Spaltenzahl)
			|| (this->m_Zeilenzahl != rhs.m_Zeilenzahl)) {
		return false;
	}
	//Gleichheit der Elemente...sowie ein ungleiches gefunden ->falsch
	for (int i = 0; i < m_Elementanzahl; i++) {
		if ((m_Elemente[i]) != (rhs.m_Elemente[i])) {
			return false;
		}
	}
	//wenns jetzt nicht foalsch war, muss es wohl richtig sein
	return true;
}

//Methoden
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// gibt eine Zeile   als Spaltenvektor aus
inline MPL_Matrix MPL_Matrix::get_Zeile(int Zeilennummer)
{
	MPL_Matrix aus(this->m_Spaltenzahl, 1);
	for (int i = 0; i < m_Spaltenzahl; i++) {
		aus(i, 0) = this->m_Elemente[i + Zeilennummer * this->m_Spaltenzahl];
	}
	return aus;
}
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
// gibt eine Spalte als Spaltenvektor aus
inline MPL_Matrix MPL_Matrix::get_Spalte(int Spaltennummer)
{
	MPL_Matrix aus(this->m_Zeilenzahl, 1);
	for (int i = 0; i < m_Zeilenzahl; i++) {
		aus(i, 0) = this->m_Elemente[Spaltennummer + i * this->m_Spaltenzahl];
	}
	return aus;
}
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline void MPL_Matrix::transpose()
{
	std::swap(m_Spaltenzahl, m_Zeilenzahl);
	transposed = !transposed;
}
inline MPL_Matrix MPL_Matrix::transponiert() const //transponierte Matrix
{
	MPL_Matrix Transponierte(*this);
	std::swap(Transponierte.m_Spaltenzahl, Transponierte.m_Zeilenzahl);
	Transponierte.transposed = !transposed;
	return Transponierte;
}
inline MPL_Matrix MPL_Matrix::transponiert_full() const
{
	//Zeilen und Spalten tauschen
	MPL_Matrix Transponierte(m_Spaltenzahl, m_Zeilenzahl);
	for (int i = 0; i < m_Spaltenzahl; i++)
		for (int j = 0; j < m_Zeilenzahl; j++)
			Transponierte.m_Elemente[j + i * m_Zeilenzahl]
				= m_Elemente[i + j * m_Spaltenzahl];
	return Transponierte;
}
inline MPL_Matrix MPL_Matrix::transponiert_full2() const
{
	MPL_Matrix Transponierte(m_Spaltenzahl, m_Zeilenzahl);
	for (int n = 0; n < m_Elementanzahl; n++) {
		int i = n / m_Zeilenzahl;
		int j = n % m_Zeilenzahl;
		Transponierte.m_Elemente[n] = m_Elemente[j * m_Spaltenzahl + i];
	}
	return Transponierte;
}
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline void MPL_Matrix::Zeile_Tauschen(int Zeile_a, int Zeile_b)
{
	if ((Zeile_a >= this->m_Zeilenzahl) || (Zeile_b >= this->m_Zeilenzahl)) {
		std::cerr << "Zeilen können nicht getauscht werden, "
			 << "weil eine Zeile nicht existiert!" << std::endl;
		return;
	}
	double dreieck;
	for (int i = 0; i < m_Spaltenzahl; i++) {
		dreieck = m_Elemente[i + Zeile_a * m_Spaltenzahl];
		m_Elemente[i + Zeile_a * m_Spaltenzahl]
			= m_Elemente[i + Zeile_b * m_Spaltenzahl];
		m_Elemente[i + Zeile_b * m_Spaltenzahl] = dreieck;
	}
}// Ende Zeile_Tauschen
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline void MPL_Matrix::Zeile_Multiplizieren(int Zeile, double Faktor)
{
	if (Zeile >= this->m_Zeilenzahl) {
		std::cerr << "Fehler in Zeile_Multiplizieren...Zeile existiert nicht"
				  << std::endl;
		return;
	}
	for (int i = 0; i < this->m_Spaltenzahl; i++) {
		this->m_Elemente[i + Zeile * m_Spaltenzahl] *= Faktor;
	}
}// Ende Zeile_Multiplizieren
/////////////////////////////////////////////////////////
// Methodenstart
/////////////////////////////////////////////////////////
inline void MPL_Matrix::Vielfaches_einer_Zeile_addieren(int Summenzeile,
		int Additionszeile, double Faktor)
{
	if ((Summenzeile >= this->m_Zeilenzahl)
			|| (Additionszeile >= this->m_Zeilenzahl)) {
		std::cerr << "Fehler in MPL_Matrix::Vielfaches_einer_Zeile_addieren..."
			 << " Zeile existiert nicht" << std::endl;
		return;
	}
	for (int i = 0; i < this->m_Spaltenzahl; i++) {
		this->m_Elemente[i + Summenzeile * m_Spaltenzahl] +=
			Faktor * this->m_Elemente[i + Additionszeile * m_Spaltenzahl];
	}
}// Ende Vielfaches_einer_Zeile_addieren

inline MPL_Matrix MPL_Matrix::row_diff()
{
	int m = m_Zeilenzahl, n = m_Spaltenzahl;
	MPL_Matrix diff(m - 1, n);

	for (int i = 0; i < m - 1; i++)
		for (int j = 0; j < n; j++)
			diff.m_Elemente[i * n + j] =
				m_Elemente[(i + 1) * n + j] - m_Elemente[i * n + j];

	return diff;
}
inline MPL_Matrix MPL_Matrix::unity() const
{
	int m = m_Zeilenzahl;
	MPL_Matrix E(m, m);

	for (int i = 0; i < m; i++)
		E.m_Elemente[i * m + i] = 1.;

	return E;
}

inline double MPL_Matrix::trace()
{
	if (m_Zeilenzahl != m_Spaltenzahl) {
		std::cerr << "Matrix is non-square." << std::endl;
		return 0.;
	}

	double trace = 0.;
	for (int i = 0; i < m_Zeilenzahl; ++i)
		trace += m_Elemente[i * m_Zeilenzahl + i];

	return trace;
}

////////////////////////////////////////////////////////////////////////////////
// Methodenstart Gausselimination_mit_Teilpivotisierung_ohne_Skalenfaktor
////////////////////////////////////////////////////////////////////////////////
// Die Funktion ist langsamer als Matlab...für 200*200 Matrizen ok,
// für 4000*4000 nicht (5min gegen 8 sekunden)
// TODO Funktion durch Lapackfunktion ersetzen oder Matlabcode einbinden
// oder selbst assembler schreiben...
// (sind ja nur 60 Zeilen code (so * 10 könnt man vll noch hinkriegen)
// TODO nachiteration...(das geht besser mit LU-Zerlegung)
inline int MPL_Matrix::Gausselimination_mit_Teilpivotisierung_ohne_Skalenfaktor()
{
	////////////////////////////////////////////////////////////////////////////
	// Die Matrix ist nicht quadratisch das Gleichungssystem
	// sieht so aus: LHS x= RHS
	// LHS ist eine gegebene Matrix, RHS ist ein gegebener Vektor
	// mehrere Zeilen der RHS werden als nebeneinander stehende Vektoren gewertet
	// Die Matrix hat die Form LHSy1y2 etc
	// Die Blockmatrix LHS wird durch Gausselimination in eine Einheitsmatrix
	// überführt und aus den Spalten y1,y2 usw
	// werden die Lösungen x1,x2, usw...d.h. LHS muss eine nxn Matrix sein und x
	// sowie RHS müssen nx1 Vektoren sein(Spaltenvektoren)
	// ACHTUNG DER LÖSUNGSVEKTOR WIRD SCHON RICHTIG DIAGONALISIERT...
	// nicht nochmal machen
	////////////////////////////////////////////////////////////////////////////
	int n = this->m_Zeilenzahl;
	//Feld in der die Reihenfolge der abgearbeiteten Zeilen stehn,
	//so spart man sich Zeilentauschs
	int *Matrixzeilen;
	Matrixzeilen = new int[n];
	for (int i = 0; i < n; i++) {
		Matrixzeilen[i] = -1;   //-1 initialisierung ...0 ist schlecht
	}

	// Schleife über alle Spalten
	for (int i = 0; i < n; i++) {
		// Suche größtes Element der Spalte
		for (int k = 0; k < n; k++) {
			// zunächst checken, ob aus dieser Zeile schon einmal
			// das größte Element stammte
			bool ignorieren = false;
			for (int j = 0; j < n; j++) {
				if (Matrixzeilen[j] == k) {
					ignorieren = true;
					continue;//j
				}
			} //for j
			if (ignorieren == true) {
				continue;     //k
			}
			//erstmögliche Zeile, die nicht ignoriert wird erstmal übernehmen
			if (Matrixzeilen[i] == -1) {
				Matrixzeilen[i] = k;
			}
			//Spaltenelemente miteinander vergleichen (Betragsmäßiger Vergleich)
			if (m_Elemente[i + Matrixzeilen[i] * m_Spaltenzahl]
					  * m_Elemente[i + Matrixzeilen[i]*m_Spaltenzahl]
					< m_Elemente[i + k * m_Spaltenzahl]
					  * m_Elemente[i + k * m_Spaltenzahl]) {
				Matrixzeilen[i] = k;
			}
		}//for k
		// Testen ob gefundenes Element 0 ist, dann Matrix singulär
		double epsilon = 1E-15;
		if ((m_Elemente[i + Matrixzeilen[i] * m_Spaltenzahl] + epsilon > 0)
				&& (m_Elemente[i + Matrixzeilen[i] * m_Spaltenzahl] - epsilon < 0)) {
			std::cerr << "Matrix ist singulär, Gausselimination scheitert"
					  << std::endl;
			//Speicher Freigeben
			delete[] Matrixzeilen;
			return 1;
		}


		// Diagonale 1 erzeugen..(division der Zeile durch das führende Element)
		double Faktor_1 = 1.0 / m_Elemente[i + Matrixzeilen[i] * m_Spaltenzahl];
//        std::cerr<<"Faktor_1: "<<Faktor_1<<"\n";
		this->Zeile_Multiplizieren(Matrixzeilen[i], Faktor_1);
//        //////////////////////////////////////////////////////////////////////
//        // zum Test Matrixzeilen und Matrix ausgeben
//        std::cerr<<"Matrixzeilen: ";
//        for(int b=0;b<n;b++)
//        {std::cerr<< Matrixzeilen[b]<<"\t";}
//        std::cerr<<"\n";
//        // Matrix ausgeben
//        {
//            for (int f=0;f<m_Zeilenzahl;f++)
//            {
//                for (int g=0;g<m_Spaltenzahl;g++)
//                {
//                    std::cerr<<m_Elemente[g+f*m_Spaltenzahl];
//                    std::cerr<<"  \t";
//                }//ende for i
//                std::cerr<<"\n";
//            }//ende for j
//            std::cerr<<"\n";
//        }//Ende Matrix_Ausgeben(MPLMatrix M)
//        //////////////////////////////////////////////////////////////////////
		// Die auserwählte Zeile von den restlichen Zeilen so abziehen,
		// dass diese in der gerade betrachteten Spalte 0 sind
		for (int k = 0; k < n; k++) {
			if (Matrixzeilen[i] == k) {
				continue; // Die Zeile soll ja die 1 haben und die anderen die 0
			}
			double Faktor_2 = -1.0 * (m_Elemente[i + k * m_Spaltenzahl]
					/ m_Elemente[i + Matrixzeilen[i] * m_Spaltenzahl]);
			this->Vielfaches_einer_Zeile_addieren(k, Matrixzeilen[i], Faktor_2);
		}// for k
	}// for i

	// Letzter Schritt...die Matrix auf echte Diagonalform bringen...
	// (sonst Lösungen in falscher Reihenfolge)
	//...nicht komplett nötig:
	//   Es reicht die RHS in die richtige Reihenfolge zu bringen
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (Matrixzeilen[j] == i) {
				//Zeilentauschen
				for (int Spalte = m_Zeilenzahl; Spalte < m_Spaltenzahl; Spalte++) {
					double d_dreieck;

					d_dreieck = m_Elemente[Spalte + Matrixzeilen[i] * m_Spaltenzahl];
					m_Elemente[Spalte + Matrixzeilen[i]*m_Spaltenzahl]
						= m_Elemente[Spalte + Matrixzeilen[j] * m_Spaltenzahl];
					m_Elemente[Spalte + Matrixzeilen[j]*m_Spaltenzahl] = d_dreieck;
				}
				//Buchhaltung updaten
				int dreieck = Matrixzeilen[j];
				Matrixzeilen[j] = Matrixzeilen[i];
				Matrixzeilen[i] = dreieck;
				break;
			}
		}//for j
	}//for i
	//Speicher Freigeben
	delete[] Matrixzeilen;
	return 0;
}
////////////////////////////////////////////////////////////////////////////////
// Ende Gausselimination_mit_Teilpivotisierung_ohne_Skalenfaktor
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Methodenstart in_Datei_speichern
////////////////////////////////////////////////////////////////////////////////
inline void MPL_Matrix::in_Datei_speichern(std::string Dateiname, double precision) const
{
	ogzstream outfile;
	MPL_Matrix A{transposed ? this->transponiert_full() : *this};
	outfile.open(Dateiname.c_str());
	if (precision != 0)
		outfile.precision(precision);
	for (int i = 0; i < A.m_Zeilenzahl; i++) {
		for (int j = 0; j < A.m_Spaltenzahl; j++) {
			outfile << m_Elemente[j + i * A.m_Spaltenzahl];
			if (j < (A.m_Spaltenzahl - 1)) {
				outfile << "\t";
			}
		}
		outfile << std::endl;
	}
	outfile.close();
	outfile.clear();
}
////////////////////////////////////////////////////////////////////////////////
// ENDE in_Datei_speichern
////////////////////////////////////////////////////////////////////////////////

/* Alternative method to store the matrix contents as a (compressed) netcdf4.
 * It should be portable and saves a few bytes on disk space.
 * If requested (pack == true), the data will be further packed using the
 * integer representation of the data (here using 32 bit unsigned int) as
 * described in http://nco.sourceforge.net/nco.html#Packed-data . */
inline int MPL_Matrix::save_to_netcdf(std::string Dateiname, bool pack) const
{
	int ncid, dimidx, dimidy, varid, ret;
	int dimids[2];
	int shuffle = NC_SHUFFLE;
	int deflate = 1;
	int deflate_level = 9;
	MPL_Matrix A{transposed ? this->transponiert_full() : *this};
	ret = nc_create(Dateiname.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid);
	if (ret) return ret;
	ret = nc_def_dim(ncid, "rows", A.m_Zeilenzahl, &dimidx);
	if (ret) return ret;
	ret = nc_def_dim(ncid, "cols", A.m_Spaltenzahl, &dimidy);
	if (ret) return ret;
	dimids[0] = dimidx;
	dimids[1] = dimidy;
	if (pack) {
		/* Manually convert to packed 32 bit integer data, see
		 * http://nco.sourceforge.net/nco.html#Packed-data
		 * for details. */
		unsigned pack_ndrv = 4294967293U; // pack ndrv for int32
		std::vector<unsigned> idata;
		auto minmax =
			std::minmax_element(m_Elemente, m_Elemente + m_Elementanzahl);
		double scale_factor = (*minmax.second - *minmax.first) / pack_ndrv;
		double add_offset = 0.5 * (*minmax.second + *minmax.first);
		std::transform(A.m_Elemente, A.m_Elemente + A.m_Elementanzahl,
				std::back_inserter(idata),
				[=](double upk) { return (upk - add_offset) / scale_factor; });
		ret = nc_def_var(ncid, "data", NC_UINT, 2, dimids, &varid);
		if (ret) return ret;
		ret = nc_put_att_double(ncid, varid, "scale_factor", NC_DOUBLE, 1, &scale_factor);
		if (ret) return ret;
		ret = nc_put_att_double(ncid, varid, "add_offset", NC_DOUBLE, 1, &add_offset);
		if (ret) return ret;
		ret = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level);
		if (ret) return ret;
		ret = nc_enddef(ncid);
		if (ret) return ret;
		ret = nc_put_var_uint(ncid, varid, idata.data());
		if (ret) return ret;
	} else {
		/* We set the netcdf variable to single precisions although the
		 * actual data are double precision. According to the libnetcdf
		 * documentation, they are converted automatically. Although we lose
		 * precision this way, it should still be fine for all practical
		 * purposes. It also saves some disk space. */
		ret = nc_def_var(ncid, "data", NC_FLOAT, 2, dimids, &varid);
		if (ret) return ret;
		ret = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level);
		if (ret) return ret;
		ret = nc_enddef(ncid);
		if (ret) return ret;
		/* We still have to tell libnetcdf that the memory block contains
		 * double precision floats. */
		ret = nc_put_var_double(ncid, varid, A.m_Elemente);
		if (ret) return ret;
	}
	ret = nc_close(ncid);
	return ret;
}

#ifdef HAVE_HDF5
/* Alternative method to store the matrix contents as a (compressed) hdf5.
 * It should be portable and saves a few bytes on disk space.
 * If requested (pack == true), the data will be further packed using the
 * integer representation of the data (here using 32 bit unsigned int) as
 * described in http://nco.sourceforge.net/nco.html#Packed-data .
 * In contrast to the netcdf4 variant above, we don't need to create
 * dimensions first and canjust save the data as is. However, we have to
 * define chunksizes by hand for the compression to work and the files seem to
 * be slightly larger than their netcdf4 counterparts. */
inline int MPL_Matrix::save_to_hdf5(std::string Dateiname, bool pack) const
{
	hid_t file_id, plist_id, dataset_id, dataspace_id;  /* identifiers */
	hsize_t dims[2], cdims[2] = { 128, 128 };
	herr_t ret;

	MPL_Matrix A{transposed ? this->transponiert_full() : *this};
	dims[0] = A.m_Zeilenzahl;
	dims[1] = A.m_Spaltenzahl;

	/* Open an existing file. */
	file_id = H5Fcreate(Dateiname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
			H5P_DEFAULT);

	dataspace_id = H5Screate_simple(2, dims, NULL);
	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	ret = H5Pset_chunk(plist_id, 2, cdims);
	if (ret) return ret;
	ret = H5Pset_shuffle(plist_id);
	if (ret) return ret;
	ret = H5Pset_deflate(plist_id, 9);
	if (ret) return ret;

	if (pack) {
		/* Manually convert to packed 32 bit integer data, see
		 * http://nco.sourceforge.net/nco.html#Packed-data
		 * for details. */
		unsigned pack_ndrv = 4294967293U; // pack ndrv for int32
		std::vector<unsigned> idata;
		auto minmax =
			std::minmax_element(m_Elemente, m_Elemente + m_Elementanzahl);
		double scale_factor = (*minmax.second - *minmax.first) / pack_ndrv;
		double add_offset = 0.5 * (*minmax.second + *minmax.first);
		hid_t attribute_id;

		std::transform(A.m_Elemente, A.m_Elemente + A.m_Elementanzahl,
				std::back_inserter(idata),
				[=](double upk) { return (upk - add_offset) / scale_factor; });

		dataset_id = H5Dcreate2(file_id, "data", H5T_STD_I32LE, dataspace_id,
				H5P_DEFAULT, plist_id, H5P_DEFAULT);

		attribute_id = H5Acreate2(dataset_id, "scale_factor", H5T_IEEE_F32LE,
				dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &scale_factor);
		if (ret) return ret;
		attribute_id = H5Acreate2(dataset_id, "add_offset", H5T_IEEE_F32LE,
				dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		ret = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &add_offset);
		if (ret) return ret;
		ret = H5Aclose(attribute_id);
		if (ret) return ret;
		/* Write the dataset. */
		ret = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
				H5P_DEFAULT, idata.data());
		if (ret) return ret;
	} else {
		/* We set the hdf5 variable to single precisions although the
		 * actual data are double precision. According to the hdf5
		 * documentation, they are converted automatically? Although we lose
		 * precision this way, it should still be fine for all practical
		 * purposes. It also saves some disk space. */
		dataset_id = H5Dcreate2(file_id, "data", H5T_IEEE_F32LE, dataspace_id,
				H5P_DEFAULT, plist_id, H5P_DEFAULT);
		/* Write the dataset. */
		ret = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
				H5P_DEFAULT, A.m_Elemente);
		if (ret) return ret;
	}
	ret = H5Dclose(dataset_id);
	if (ret) return ret;
	ret = H5Pclose(plist_id);
	if (ret) return ret;
	ret = H5Sclose(dataspace_id);
	if (ret) return ret;
	ret = H5Fclose(file_id);
	return ret;
}
#else /* HAVE_HDF5 */
/* Dummy method if hdf5 is not available on the system */
inline int MPL_Matrix::save_to_hdf5(std::string Dateiname, bool pack) const
{
	std::cerr << "Saving to hdf5 is not available in this build, "
			<< "use ascii or netcdf instead." << std::endl;
	return -1;
}
#endif /* HAVE_HDF5 */

#endif /* MPLMATRIX_HH_ */
