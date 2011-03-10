/*
 * Gitterpunkt.cpp
 *
 *  Created on: 27.05.2010
 *      Author: martin
 */
#include "Gitterpunkt.h"

// Konstruktoren
Gitterpunkt::Gitterpunkt()
{
	//leer
	this->m_vorderer_Durchstosspunkt.Elementzahl_festlegen(3);
	this->m_hinterer_Durchstosspunkt.Elementzahl_festlegen(3);

}
Gitterpunkt::Gitterpunkt(const Gitterpunkt &rhs)
{
	*this = rhs; //Siehe operator =
}
// Überladene Operatoren
Gitterpunkt &Gitterpunkt::operator =(const Gitterpunkt &rhs)
{
	if (this == &rhs)
		return *this;
	//Nur statische Variablen
	this->m_Breite = rhs.m_Breite;
	this->m_Hoehe = rhs.m_Hoehe;
	this->m_Index_Nord_Nachbar              = rhs.m_Index_Nord_Nachbar;
	this->m_Index_Sued_Nachbar             = rhs.m_Index_Sued_Nachbar;
	this->m_Index_oberer_Nachbar           = rhs.m_Index_oberer_Nachbar;
	this->m_Index_unterer_Nachbar          = rhs.m_Index_unterer_Nachbar;
	this->m_Index_oberer_Nord_Nachbar   = rhs.m_Index_oberer_Nord_Nachbar ;
	this->m_Index_unterer_Nord_Nachbar  = rhs.m_Index_unterer_Nord_Nachbar;
	this->m_Index_oberer_Sued_Nachbar  = rhs.m_Index_oberer_Sued_Nachbar ;
	this->m_Index_unterer_Sued_Nachbar = rhs.m_Index_unterer_Sued_Nachbar;
	this->m_Max_Breite = rhs.m_Max_Breite;
	this->m_Max_Hoehe = rhs.m_Max_Hoehe;
	this->m_Min_Breite = rhs.m_Min_Breite;
	this->m_Min_Hoehe = rhs.m_Min_Hoehe;
	this->m_eigener_Index = rhs.m_eigener_Index;
	// Durchstoßpunkthoehen
	this->m_vorderer_Durchstosspunkt = rhs.m_vorderer_Durchstosspunkt;
	this->m_hinterer_Durchstosspunkt = rhs.m_hinterer_Durchstosspunkt;


	//this->m_SZA=rhs.m_SZA; // TODO Winkel in Raytracing...hier später löschen
	return *this;
}
//Methoden
bool Gitterpunkt::Punkt_in_Gitterpunkt(double Lat, double Hoehe)
{
	if ((Lat <= m_Max_Breite) && (Lat > m_Min_Breite) && (Hoehe <= m_Max_Hoehe) && (Hoehe > m_Min_Hoehe)) {
		//Punkt im Gitter
		return true;
	}
	return false;
}
