/*
 * MPL_Vektor.h
 *
 *  28.09.2010
 *  skalare division...da stand * statt /
 *
 *
 *  Stand 9.09.2010
 *  neu Betrag()
 *        normieren()
 *
 *
 *
 *
 *  Created on: 15.06.2010
 *      Author: martin
 */
// Es ist doch Sinnvoll neben der Matrixklasse eine mathematische Vektorklasse
// einzuführen, da Vektoren einer viel einfacheren
// Arithmetik unterliegen

#ifndef MPL_VEKTOR_HH_
#define MPL_VEKTOR_HH_

#include <cmath>
#include <iostream>

class MPL_Vektor
{
public:
	// Konstruktoren ////////
	MPL_Vektor();
	MPL_Vektor(int Elementzahl);
	MPL_Vektor(const MPL_Vektor &rhs);  //copyconstructor
	// Destruktor /////////////
	~MPL_Vektor();
	// Überladene Operatoren //////
	MPL_Vektor &operator =(const MPL_Vektor &rhs);     //= Zuweisung
	MPL_Vektor &operator +=(const MPL_Vektor &rhs);   // Vektoraddition
	MPL_Vektor &operator -=(const MPL_Vektor &rhs);    // Vektorsubtraktion

	MPL_Vektor &operator *=(const double &rhs);  // skalare Multiplikation
	MPL_Vektor &operator /=(const double &rhs);          // skalare Division

	// Zugriff auf Elemente
	double &operator()(int Element);

	//binäre Operationen
	MPL_Vektor operator +(const MPL_Vektor &rhs);  //Vektoraddition
	MPL_Vektor operator -(const MPL_Vektor &rhs);   //Vektorsubtraktion
	MPL_Vektor operator *(const double &rhs);         // skalare Multiplikation
	MPL_Vektor operator /(const double &rhs);         // skalare Division
	double operator *(const MPL_Vektor &rhs);         //Skalarprodukt

	friend MPL_Vektor operator * (const double &lhs, const MPL_Vektor &rhs);

	bool operator == (const MPL_Vektor &rhs) const;
	bool operator != (const MPL_Vektor &rhs) const;

	// Methoden
	void Elementzahl_festlegen(int E);
	void Null_Initialisierung();
	double Betrag_ausgeben();
	void Normieren();
	//MPL_Vektor Kreuzprodukt(const MPL_Vektor& rhs);
	// nur für 3er Vektoren sinnvoll
	//Membervariablen
	int m_Elementanzahl;
	double *m_Elemente;
};

//Implementation der Inlinefunktionen

// Konstruktoren ////////
inline MPL_Vektor::MPL_Vektor()
{
	m_Elementanzahl = 0;
	m_Elemente = 0;
}// Ende defaultkonstruktor
inline MPL_Vektor::MPL_Vektor(int Elementzahl)
{
	m_Elementanzahl = Elementzahl;
	m_Elemente = new double[Elementzahl];
	for (int i = 0; i < Elementzahl; i++) {
		m_Elemente[i] = 0;
	}
} //Ende MPL_Vektor::MPL_Vektor(int Elementzahl)

inline MPL_Vektor::MPL_Vektor(const MPL_Vektor &rhs)  //copyconstructor
{
	m_Elementanzahl = 0;
	m_Elemente = 0;
	*this = rhs;
}
// Destruktor /////////////
inline MPL_Vektor::~MPL_Vektor()
{
	if (m_Elemente != 0) {
		delete[] m_Elemente;
	}
}// Ende Destruktor

// Überladene Operatoren //////
//= Zuweisung
inline MPL_Vektor &MPL_Vektor::operator =(const MPL_Vektor &rhs)
{
	if (this == &rhs)
		return *this;
	m_Elementanzahl = rhs.m_Elementanzahl;
	if (m_Elemente != 0) {
		delete[] m_Elemente;
	}
	m_Elemente = new double[this->m_Elementanzahl];
	for (int i = 0; i < m_Elementanzahl; i++) {
		m_Elemente[i] = rhs.m_Elemente[i];
	}
	return *this;
}// Ende =

// Vektoraddition
inline MPL_Vektor &MPL_Vektor::operator +=(const MPL_Vektor &rhs)
{
	if (this->m_Elementanzahl != rhs.m_Elementanzahl) {
		std::cerr << "Vektoraddition verschieden langer Vektoren ist unsinnig!!!!"
				  << std::endl;
		return *this;
	}
	for (int i = 0; i < m_Elementanzahl; i++) {
		this->m_Elemente[i] += rhs.m_Elemente[i];
	}
	return *this;
}// Ende +=

// Vektorsubtraktion
inline MPL_Vektor &MPL_Vektor::operator -=(const MPL_Vektor &rhs)
{
	if (this->m_Elementanzahl != rhs.m_Elementanzahl) {
		std::cerr << "Vektorsubtraktion verschieden langer Vektoren ist unsinnig!!!!"
				  << std::endl;
		return *this;
	}
	for (int i = 0; i < m_Elementanzahl; i++) {
		this->m_Elemente[i] -= rhs.m_Elemente[i];
	}
	return *this;
}// Ende -=

// skalare Multiplikation
inline MPL_Vektor &MPL_Vektor::operator *=(const double &rhs)
{
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		m_Elemente[i] *= rhs;
	}
	return *this;
}//Ende *=

// skalare Division
inline MPL_Vektor &MPL_Vektor::operator /=(const double &rhs)
{
	if (rhs == 0) { // Achtung ...eigentlich in epsilonumgebung betrachten
		std::cerr << "Achtung Division durch 0 bei skalarer Vektordivision"
				  << std::endl;
		return *this;
	}

	for (int i = 0; i < this->m_Elementanzahl; i++) {
		m_Elemente[i] /= rhs;
	}
	return *this;
}//Ende /=

// Zugriff auf Elemente
// Die braucht man aus obskuren Gründen 2mal
inline double &MPL_Vektor::operator()(int Element)
{
	if ((Element >= 0) && (Element < m_Elementanzahl))
		return this->m_Elemente[Element];
	else {
		std::cerr << "Achtung!!! Zugriff auf Elemente ausserhalb des Vektors"
				  << std::endl;
		return m_Elemente[0];//auch schlecht, aber wenigstens nicht ausserhalb
	}
}

//binäre Operationen
//Vektoraddition
inline MPL_Vektor MPL_Vektor::operator +(const MPL_Vektor &rhs)
{
	if (this->m_Elementanzahl != rhs.m_Elementanzahl) {
		std::cerr << "Vektoraddition verschieden langer Vektoren ist unsinnig!!!!"
				  << std::endl;
		return *this;
	}
	MPL_Vektor aus(m_Elementanzahl);
	for (int i = 0; i < m_Elementanzahl; i++) {
		aus.m_Elemente[i] = this->m_Elemente[i] + rhs.m_Elemente[i];
	}
	return aus;
}// Ende +
//Vektorsubtraktion
inline MPL_Vektor MPL_Vektor::operator -(const MPL_Vektor &rhs)
{
	if (this->m_Elementanzahl != rhs.m_Elementanzahl) {
		std::cerr << "Vektorsubtraktion verschieden langer Vektoren ist unsinnig!!!!"
				  << std::endl;
		return *this;
	}
	MPL_Vektor aus(m_Elementanzahl);
	for (int i = 0; i < m_Elementanzahl; i++) {
		aus.m_Elemente[i] = this->m_Elemente[i] - rhs.m_Elemente[i];
	}
	return aus;
}//Ende -

//Skalarprodukt
inline double MPL_Vektor::operator *(const MPL_Vektor &rhs)
{
	if (this->m_Elementanzahl != rhs.m_Elementanzahl) {
		std::cerr << "Skalarprodukt verschieden langer Vektoren ist unsinnig!!!!"
				  << std::endl;
		return 0;
	}
	double aus = 0;
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		aus += this->m_Elemente[i] * rhs.m_Elemente[i];
	}
	return aus;
}// Ende Skalarprodukt
// skalare Division
inline MPL_Vektor MPL_Vektor::operator /(const double &rhs)
{
	if (rhs == 0) {
		std::cerr << "Achtung Division durch 0!!!!!" << std::endl;
		return *this;
	}
	MPL_Vektor aus(m_Elementanzahl);
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		aus.m_Elemente[i] = this->m_Elemente[i] / rhs;
	}
	return aus;

}
// skalare Multiplikation von skalar rechts
inline MPL_Vektor MPL_Vektor::operator *(const double &rhs)
{
	MPL_Vektor aus(m_Elementanzahl);
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		aus.m_Elemente[i] = this->m_Elemente[i] * rhs;
	}
	return aus;
}

// skalare Multiplikation von skalar links
inline MPL_Vektor operator * (const double &lhs, const MPL_Vektor &rhs)
{
	MPL_Vektor aus(rhs.m_Elementanzahl);
	for (int i = 0; i < rhs.m_Elementanzahl; i++) {
		aus.m_Elemente[i] = rhs.m_Elemente[i] * lhs;
	}
	return aus;
}

inline bool MPL_Vektor::operator == (const MPL_Vektor &rhs) const
{
	//Zwei Vektoren sind gleich, wenn alle ihre Elemente gleich sind
	if (this->m_Elementanzahl != rhs.m_Elementanzahl)
		return false;
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		if (this->m_Elemente[i] != rhs.m_Elemente[i])
			return false;
	}
	return true;
}// Ende ==
inline bool MPL_Vektor::operator != (const MPL_Vektor &rhs) const
{
	if (*this == rhs)
		return false;
	return true;
}// Ende !=

// Methoden

inline void MPL_Vektor::Elementzahl_festlegen(int E)
{
	// Falls Vektor reduziert wird, die ersten Einträge behalten
	double *altes_Feld = 0;
	int altes_Feld_Elementzahl = m_Elementanzahl;
	//altes Feld für Dreieckstauch zwischenkopieren
	if (m_Elementanzahl != 0) {
		altes_Feld = new double[m_Elementanzahl];
		for (int i = 0; i < m_Elementanzahl; i++) {
			altes_Feld[i] = m_Elemente[i];
		}
		//Altes Feld löschen
		delete[] m_Elemente;
		m_Elemente = 0;        //Zeile eigentlich überflüssig
		m_Elementanzahl = 0; //Zeile eigentlich überflüssig
	}
	//neuen Speicher allokieren
	m_Elementanzahl = E;
	m_Elemente = new double[m_Elementanzahl];
	//neues Feld mit altem so weit es geht von vorne auffüllen...
	//falls Rest, mit 0 auffüllen
	int Grenze;
	if (m_Elementanzahl < altes_Feld_Elementzahl) {
		Grenze = m_Elementanzahl;
	} else {
		Grenze = altes_Feld_Elementzahl;
	}
	for (int i = 0; i < Grenze; i++) {
		m_Elemente[i] = altes_Feld[i];
	}
	for (int i = Grenze; i < m_Elementanzahl; i++) {
		m_Elemente[i] = 0;
	}
	//altes Feld löschen
	if (altes_Feld != 0) {
		delete[] altes_Feld;
	}
}

inline void MPL_Vektor::Null_Initialisierung()
{
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		this->m_Elemente[i] = 0;

	}
}// Ende Nullinitialisierung
inline double MPL_Vektor::Betrag_ausgeben()
{
	double d_Betrag = 0;
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		d_Betrag += m_Elemente[i] * m_Elemente[i];
	}
	d_Betrag = sqrt(d_Betrag);
	return d_Betrag;
}

inline void MPL_Vektor::Normieren()
{
	double d_Betrag = Betrag_ausgeben();
	if (d_Betrag < 0.000001) {
		std::cerr << "Es wird versucht einen Vektor zu normieren, "
			 << "dessen Betrag <0.000001 ist. Dies wird nicht gemacht"
			 << std::endl;
		return; // Nullvektor abbbruch
	}
	for (int i = 0; i < this->m_Elementanzahl; i++) {
		m_Elemente[i] /= d_Betrag;
	}
}

#endif /* MPL_VEKTOR_HH_ */
