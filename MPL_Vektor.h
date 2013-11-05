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
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <iterator>

class MPL_Vektor
{
public:
	// Konstruktoren ////////
	MPL_Vektor() : m_Elemente(0) {}
	explicit MPL_Vektor(unsigned long Elementzahl);
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

	// Methoden
	void Elementzahl_festlegen(int E);
	void Null_Initialisierung();
	double Betrag_ausgeben() const;
	void Normieren();
	//MPL_Vektor Kreuzprodukt(const MPL_Vektor& rhs);
	// nur für 3er Vektoren sinnvoll
private:
	//Membervariablen
	std::vector<double> m_Elemente;
};

//Implementation der Inlinefunktionen

inline MPL_Vektor::MPL_Vektor(unsigned long Elementzahl) :
	m_Elemente(Elementzahl)
{ } //Ende MPL_Vektor::MPL_Vektor(int Elementzahl)

// Überladene Operatoren //////
//= Zuweisung
inline MPL_Vektor &MPL_Vektor::operator =(const MPL_Vektor &rhs)
{
	if (this == &rhs)
		return *this;

	m_Elemente.resize(rhs.m_Elemente.size());
	std::copy(rhs.m_Elemente.begin(), rhs.m_Elemente.end(),
			m_Elemente.begin());

	return *this;
}// Ende =

// Vektoraddition
inline MPL_Vektor &MPL_Vektor::operator +=(const MPL_Vektor &rhs)
{
	std::transform(m_Elemente.begin(), m_Elemente.end(),
			rhs.m_Elemente.begin(), m_Elemente.begin(),
			std::plus<double>());

	return *this;
}// Ende +=

// Vektorsubtraktion
inline MPL_Vektor &MPL_Vektor::operator -=(const MPL_Vektor &rhs)
{
	std::transform(m_Elemente.begin(), m_Elemente.end(),
			rhs.m_Elemente.begin(), m_Elemente.begin(),
			std::minus<double>());

	return *this;
}// Ende -=

// skalare Multiplikation
inline MPL_Vektor &MPL_Vektor::operator *=(const double &rhs)
{
	std::transform(m_Elemente.begin(), m_Elemente.end(),
			m_Elemente.begin(),
			std::bind2nd(std::multiplies<double>(), rhs));

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

	std::transform(m_Elemente.begin(), m_Elemente.end(),
			m_Elemente.begin(),
			std::bind2nd(std::divides<double>(), rhs));

	return *this;
}//Ende /=

// Zugriff auf Elemente
// Die braucht man aus obskuren Gründen 2mal
inline double &MPL_Vektor::operator()(int Element)
{
	return m_Elemente.at(Element);
}

//binäre Operationen
//Vektoraddition
inline MPL_Vektor MPL_Vektor::operator +(const MPL_Vektor &rhs)
{
	MPL_Vektor aus(m_Elemente.size());

	std::transform(m_Elemente.begin(), m_Elemente.end(),
			rhs.m_Elemente.begin(), aus.m_Elemente.begin(),
			std::plus<double>());

	return aus;
}// Ende +
//Vektorsubtraktion
inline MPL_Vektor MPL_Vektor::operator -(const MPL_Vektor &rhs)
{
	MPL_Vektor aus(m_Elemente.size());

	std::transform(m_Elemente.begin(), m_Elemente.end(),
			rhs.m_Elemente.begin(), aus.m_Elemente.begin(),
			std::minus<double>());

	return aus;
}//Ende -

//Skalarprodukt
inline double MPL_Vektor::operator *(const MPL_Vektor &rhs)
{
	return std::inner_product(m_Elemente.begin(), m_Elemente.end(),
			rhs.m_Elemente.begin(), 0.0);
}// Ende Skalarprodukt
// skalare Division
inline MPL_Vektor MPL_Vektor::operator /(const double &rhs)
{
	if (rhs == 0) {
		std::cerr << "Achtung Division durch 0!!!!!" << std::endl;
		return *this;
	}

	MPL_Vektor aus(m_Elemente.size());

	std::transform(m_Elemente.begin(), m_Elemente.end(),
			aus.m_Elemente.begin(),
			std::bind2nd(std::divides<double>(), rhs));

	return aus;
}
// skalare Multiplikation von skalar rechts
inline MPL_Vektor MPL_Vektor::operator *(const double &rhs)
{
	MPL_Vektor aus(m_Elemente.size());

	std::transform(m_Elemente.begin(), m_Elemente.end(),
			aus.m_Elemente.begin(),
			std::bind2nd(std::multiplies<double>(), rhs));

	return aus;
}

// skalare Multiplikation von skalar links
inline MPL_Vektor operator * (const double &lhs, const MPL_Vektor &rhs)
{
	MPL_Vektor aus(rhs);

	std::transform(rhs.m_Elemente.begin(), rhs.m_Elemente.end(),
			aus.m_Elemente.begin(),
			std::bind2nd(std::multiplies<double>(), lhs));

	return aus;
}

inline bool MPL_Vektor::operator == (const MPL_Vektor &rhs) const
{
	//Zwei Vektoren sind gleich, wenn alle ihre Elemente gleich sind
	if (this->m_Elemente.size() != rhs.m_Elemente.size())
		return false;

	return std::equal(m_Elemente.begin(), m_Elemente.end(),
			rhs.m_Elemente.begin());
}// Ende ==

// Methoden

inline void MPL_Vektor::Elementzahl_festlegen(int E)
{
	m_Elemente.resize(E);
}

inline void MPL_Vektor::Null_Initialisierung()
{
	std::fill(m_Elemente.begin(), m_Elemente.end(), 0.0);
}// Ende Nullinitialisierung

inline double MPL_Vektor::Betrag_ausgeben() const
{
	return std::sqrt(std::inner_product(m_Elemente.begin(),
				m_Elemente.end(), m_Elemente.begin(), 0.0));
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

	std::transform(m_Elemente.begin(), m_Elemente.end(),
			m_Elemente.begin(),
			std::bind2nd(std::divides<double>(), d_Betrag));
}

#endif /* MPL_VEKTOR_HH_ */
