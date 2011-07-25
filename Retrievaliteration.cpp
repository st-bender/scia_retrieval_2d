/*
 * Retrievaliteration.cpp
 *
 *  Created on: 14.09.2010
 *      Author: martin
 */

#include"Retrievaliteration.h"

#include"MPL_Matrix.h"
#include "Konfiguration.h"
#include <cmath>

using namespace std;

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B,
			int *LDB, int *INFO);
	void dgetrs_(char *, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
			double *B, int *LDB, int *INFO);
}

int Retrievaliteration(MPL_Matrix &Dichten,
					   MPL_Matrix &Dichten_apriori,
					   MPL_Matrix &Saeulendichten,
					   MPL_Matrix &S_apriori,
					   MPL_Matrix &S_y,
					   MPL_Matrix &S_Breite,
					   MPL_Matrix &S_Hoehe,
					   const double &Lambda_Breite,
					   const double &Lambda_Hoehe,
					   MPL_Matrix &AMF,
					   Konfiguration &Konf)
{
	// Die Sache mit dem Apriori ist noch glaub ich noch nicht ganz sauber...
	// bei 0 apriori kein problem beim nachiterieren

	// Kritik an alter Version;
	// Die Iteration ist schlecht.... apriori sollte nicht neu gesetzt werden
	// Lieber eine iteration mit dem Rest, der noch da ist
	// Das wird nun gemacht
	// Debugging Output aus alter Routine wurde entfernt...
	// kann von da aber hier rein kopiert werden

	//linke Seite der Gleichung(anzuwenden auf Dichten, um RHS zu erhalten)
	MPL_Matrix LHS;
	MPL_Matrix AMF_trans = AMF.transponiert();
	MPL_Matrix S_Breite_trans = S_Breite.transponiert();
	MPL_Matrix S_Hoehe_trans = S_Hoehe.transponiert();

	// Solange man die lambdas für die constraints nicht ändern will,
	// sieht die LHS immer gleich aus
	LHS = (AMF_trans * (S_y * AMF));
	LHS += (S_apriori);
	LHS += (Lambda_Breite * (S_Breite_trans * S_Breite));  // Breitenglattung
	LHS += (Lambda_Hoehe * (S_Hoehe_trans * S_Hoehe)); // Hoehenglattung
//    cout<<"LHS: "<<LHS.m_Zeilenzahl<<"\t"<<LHS.m_Spaltenzahl<<"\n";
	////////////////////////////////////////////////////////////////////////////
	// TODO Der Absatz muss neu geschrieben werden, weil stimmt nichtmehr
	// ich halte mich erstmal an Marcos Routine, bis auf die Diagonalisierung
	// der Matrix...das mach ich mit Gauss statt cholesky, da auch negative
	// Dichten als Lösungen prinzipiell erstmal zugelassen werden sollen (das
	// sagt einem ja dann auch was über die Genauigkeit des Retrievals aus,
	// ausserdem überschätzt man sonst die Dichte)

	// Die Iteration in Marcos Programm läuft folgendermaßen ab
	// In der Hauptschleife(Iterationsschleife)
	// wird der Solver der Normal_Gleichung aufgerufen
	// danach werden die Residuen berechnet und es wird überprüft, ob diese zum
	// vorherigen Schritt kleiner geworden sind falls nicht, wird abgebrochen
	// falls die Residuen kleiner geworden sind, werden alle positiven Einträge
	// von X als neues apriori übergeben( die anderen bleiben gleich)

	//Der Solver löst die Gleichung mehrmals, falls IERR nicht 1 ist, d.h. die
	//LHS nicht positiv definit, (dann wird Lambda_apriori erhöht, was Unsinn
	//ist, da das nur in die RHS einfließt, die keine Rolle bei der
	//feststellung der Definitheit der LHS hat.
	//Das heißt, falls es kene Lösung gibt, braucht das Programm auch noch
	//länger, da ja das Flag so nie 1 werden kann Diesen Teil lass ich also weg
	//und mach nur die Iteration
	//Nach meinem letzten Vortrag hieß es, dass eine derartige Anpassung des
	//aprior auch nicht dichter an der Lösung liegt....  Evtl ist das apriori
	//eher sowas wie der Iterationsanfang

	//Als Solver für das Programm werden die ATLAS/LAPACK Routinen für die
	//LU-Zerlegung einer Matrix, und das Rückeinsetzen mit hilfe dieser
	//Zerlegung verwendet.
	//Diese implementation gehört zu den schnellsten ihrer art.
	//Es handelt sich dabei um 2 Schritte...die LU-Zerlegung selbst und die
	//Rückeinsetzung Die LU-Zerlegung braucht Zeit, das Rückeinsetzen erfolgt
	//dagegen fast instantan für positiv definite Matrizen ist
	//Choleskyzerlegung genauer, das schließt aber negative Säulendichten aus,
	//also wird LU-benutzt
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// LU Zerlegung der LHS
	////////////////////////////////////////////////////////////////////////////
	MPL_Matrix Dichten_alt = Dichten;
	MPL_Matrix RHS(Dichten.m_Zeilenzahl, 1);
	MPL_Matrix Saeulendichten_neu(Saeulendichten.m_Zeilenzahl, 1);
	MPL_Matrix Saeulendichten_rest(Saeulendichten.m_Zeilenzahl, 1);
	//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";
	for (int i = 0; i < RHS.m_Elementanzahl; i++) {
		RHS(i) = 1;    //Dummy Rechte Seite erstellen
	}
	// Zunächst die LU Zerlegung der LHS durchführen mit dummy RHS,
	// in der Iteration dann nur das Rückeinsetzen nutzen
	// START LU ZERLEGUNG
	// VORBEREITEN
	//Fortran Matrizen sind zu C++ Matrizen transponiert
	MPL_Matrix A = LHS.transponiert();
	// Man kann auch die MPL_Matrix nach Fortran Nomenklatur anpassen, aber
	// transponieren ist nicht zeitaufwändig
	// Feldgröße Speed propto N^3, LHS ist quadratisch, N ist Anzahl
	// der Gitterpunkte
	int N = LHS.m_Zeilenzahl;
	//array mit der Pivotisierungsmatrix sollte so groß wie N sein,
	//alle Elemente 0
	int *IPIV;
	IPIV = new int[N];
	// ------ RHS oben definiert
	//Spalten von RHS 1 nehmen, um keine c/Fortran Verwirrungen zu provozieren
	int NRHS = 1;
	int LDA = N;
	int LDB = N;
	int INFO;
	//Anzahl sollte die Integer grenzen nicht überschreiten,
	//aber danbn sollte der Aufbau von LHS schon stören
	//int Anzahl=N*N;
	char textflag = 'N'; // "No transpose"; //fürs Rückeinsetzen
	// AUFRUF   A ist LHS.transponiert und B ist RHS
	dgesv_(&N, &NRHS, A.m_Elemente, &LDA, IPIV, RHS.m_Elemente, &LDB, &INFO);
	// TODO im Prinizp muss man RHS, falls es eine Matrix ist,
	// jetzt auch transponieren
	// ENDE LU ZERLEGUNG
	//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";
	////////////////////////////////////////////////////////////////////////////
	// ENDE LU Zerlegung der LHS
	////////////////////////////////////////////////////////////////////////////

	double Residual, Residual_1, residual_prev = 0.001;
	MPL_Matrix Mat_Residual;
	Residual = 0;
	Residual_1 = 0;
	int Itmax = Konf.m_Max_Zahl_Iterationen;
	double Threshold = Konf.m_Convergence_Treshold;
	//Threshold=1E-5;  // für 12 Hoehen gut
	//Threshold=1E-2;
	//cout<<"Itmax: "<<Itmax<<"\n";
	//MPL_Matrix RHS_Teil1;
	//RHS_Teil1=AMF_trans*(S_y*Saeulendichten);
	//Die Schleife ist echt schnell...das sollten höchstens 10 sekunden sein
	///////////////////////////////////////
	// erster Schritt
	///////////////////////////////////////
	//RHS sollte ein Spaltenvektor sein
	RHS = AMF_trans * (S_y * Saeulendichten) + S_apriori * Dichten_apriori;
	// Lösungen durch Rückeinsetzen finden
	dgetrs_(&textflag, &N, &NRHS, A.m_Elemente, &LDA, IPIV, RHS.m_Elemente,
			&LDB, &INFO);
	Dichten = RHS;
	///////////////////////////////////////
	// ENDE erster Schritt
	///////////////////////////////////////

	///////////////////////////////////////
	// Nachiteration
	///////////////////////////////////////
	for (int Iterationsschritt = 0; Iterationsschritt < Itmax; Iterationsschritt++) {
		Saeulendichten_neu = AMF * Dichten;
		Saeulendichten_rest = Saeulendichten - Saeulendichten_neu;
		MPL_Matrix Dichten_apriori_rest = Dichten_apriori - Dichten;
		Mat_Residual = (Saeulendichten_rest.transponiert() * S_y * Saeulendichten_rest)
			+ Dichten_apriori_rest.transponiert() * S_apriori * Dichten_apriori_rest;

		// Anfangsresiduum bestimmen
		Residual = sqrt(Mat_Residual(0));
		if ((Iterationsschritt == 0) || (Iterationsschritt == Itmax - 1)) {
			cerr << "Iterationsschritt:" << Iterationsschritt << "\t"
				 << "Residual: " << Residual << "\n";
		}
		if (Iterationsschritt == 0) {
			// erstes Residuum als ungefähre Fehlerabschätzung
			Residual_1 = Residual;
		}
		if (Residual < Threshold * Residual_1) {
			cerr << "Konvergenz bei Iterationsschritt: " << Iterationsschritt << endl;
			cerr << "Residual: " << Residual << endl;
			break;
		}
		if (abs(residual_prev - Residual) / residual_prev < Threshold) {
			cerr << "Konvergenz bei Iterationsschritt: " << Iterationsschritt << endl;
			cerr << "Residual: " << Residual << endl;
			break;
		}
		residual_prev = Residual;
		//RHS sollte ein Spaltenvektor sein
		RHS = AMF_trans * (S_y * Saeulendichten_rest)
			  + S_apriori * Dichten_apriori_rest;
		// Lösungen durch Rückeinsetzen finden
		dgetrs_(&textflag, &N, &NRHS, A.m_Elemente, &LDA, IPIV, RHS.m_Elemente,
				&LDB, &INFO);
		Dichten += RHS; // inkrementierung
	}//for Iterationsschritt
	// keine Probleme während Iteration aufgetreten
	//Dichten.in_Datei_speichern("/tmp/mlangowski/0/Dichten_nach_Iteration.txt");
	// dynamischen Kram löschen
	delete[] IPIV;
	return 0;
}



int Retrievaliteration_old(MPL_Matrix &Dichten,
						   MPL_Matrix &Dichten_apriori,
						   MPL_Matrix &Saeulendichten,
						   MPL_Matrix &S_apriori,
						   MPL_Matrix &S_y,
						   MPL_Matrix &S_Breite,
						   MPL_Matrix &S_Hoehe,
						   MPL_Matrix &S_letzte_Hoehe,
						   const double &Lambda_Breite,
						   const double &Lambda_Hoehe,
						   MPL_Matrix &AMF,
						   Konfiguration &Konf)
{
	// Die Iteration ist schlecht.... apriori sollte nicht neu gesetzt werden
	// Lieber eine iteration mit dem Rest, der noch da ist

	//linke Seite der Gleichung(anzuwenden auf Dichten, um RHS zu erhalten)
	MPL_Matrix LHS;
	MPL_Matrix AMF_trans = AMF.transponiert();
	MPL_Matrix S_Breite_trans = S_Breite.transponiert();
	MPL_Matrix S_Hoehe_trans = S_Hoehe.transponiert();
	MPL_Matrix S_letzte_Hoehe_trans = S_letzte_Hoehe.transponiert();

//    cout<<"Dim S_y: "<<S_y.m_Zeilenzahl<<"\t"<<S_y.m_Spaltenzahl<<"\n";
//    cout<<"AMF: "<<AMF.m_Zeilenzahl<<"\t"<<AMF.m_Spaltenzahl<<"\n";
//    cout<<"AMF_trans: "<<AMF_trans.m_Zeilenzahl<<"\t"<<AMF_trans.m_Spaltenzahl<<"\n";
	LHS = (AMF_trans * (S_y * AMF));
	LHS += (S_apriori);
	LHS += (Lambda_Breite * (S_Breite_trans * S_Breite));  // Breitenglattung
	LHS += (Lambda_Hoehe * (S_Hoehe_trans * S_Hoehe)); // Hoehenglattung
	LHS += S_letzte_Hoehe_trans * S_letzte_Hoehe;  // letzte Hoehe auf 0 zwingen
//    cout<<"LHS: "<<LHS.m_Zeilenzahl<<"\t"<<LHS.m_Spaltenzahl<<"\n";
	////////////////////////////////////////////////////////////////////////////
	// ich halte mich erstmal an Marcos Routine, bis auf die Diagonalisierung
	// der Matrix...das mach ich mit Gauss statt cholesky, da auch negative
	// Dichten als Lösungen prinzipiell erstmal zugelassen werden sollen(das
	// sagt einem ja dann auch was über die Genauigkeit des Retrievals aus,
	// ausserdem überschätzt man sonst die Dichte)

	// Die Iteration in Marcos Programm läuft folgendermaßen ab
	// In der Hauptschleife(Iterationsschleife)
	// wird der Solver der Normal_Gleichung aufgerufen
	// danach werden die Residuen berechnet und es wird überprüft, ob diese zum
	// vorherigen Schritt kleiner geworden sind falls nicht, wird abgebrochen
	// falls die Residuen kleiner geworden sind, werden alle positiven Einträge
	// von X als neues apriori übergeben( die anderen bleiben gleich)

	//Der Solver löst die Gleichung mehrmals, falls IERR nicht 1 ist, d.h. die
	//LHS nicht positiv definit, (dann wird Lambda_apriori erhöht, was Unsinn
	//ist, da das nur in die RHS einfließt, die keine Rolle bei der
	//feststellung der Definitheit der LHS hat.
	//Das heißt, falls es kene Lösung gibt, braucht das Programm auch noch
	//länger, da ja das Flag so nie 1 werden kann Diesen Teil lass ich also weg
	//und mach nur die Iteration
	//Nach meinem letzten Vortrag hieß es, dass eine derartige Anpassung des
	//aprior auch nicht dichter an der Lösung liegt....  Evtl ist das apriori
	//eher sowas wie der Iterationsanfang

	//Als Solver für das Programm werden die ATLAS/LAPACK Routinen für die
	//LU-Zerlegung einer Matrix, und das Rückeinsetzen mit hilfe dieser
	//Zerlegung verwendet.
	//Diese implementation gehört zu den schnellsten ihrer art.
	//Es handelt sich dabei um 2 Schritte...die LU-Zerlegung selbst und die
	//Rückeinsetzung Die LU-Zerlegung braucht Zeit, das Rückeinsetzen erfolgt
	//dagegen fast instantan für positiv definite Matrizen ist
	//Choleskyzerlegung genauer, das schließt aber negative Säulendichten aus,
	//also wird LU-benutzt
	////////////////////////////////////////////////////////////////////////////
	int Itmax = Konf.m_Max_Zahl_Iterationen;
	double Threshold = Konf.m_Convergence_Treshold;
	//Threshold=1E-5;  // für 12 Hoehen gut
	MPL_Matrix Dichten_alt = Dichten;
	MPL_Matrix RHS(Dichten.m_Zeilenzahl, 1);
	//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";
	for (int i = 0; i < RHS.m_Elementanzahl; i++) {
		RHS(i) = 1;    //Dummy Rechte Seite erstellen
	}
	// Zunächst die LU Zerlegung der LHS durchführen mit dummy RHS,
	// in der Iteration dann nur das Rückeinsetzen nutzen
	// START LU ZERLEGUNG
	// VORBEREITEN
	//Fortran Matrizen sind zu C++ Matrizen transponiert
	MPL_Matrix A = LHS.transponiert();
	// Man kann auch die MPL_Matrix nach Fortran Nomenklatur anpassen, aber
	// transponieren ist nicht zeitaufwändig
	// Feldgröße Speed propto N^3, LHS ist quadratisch,
	// N ist Anzahl der Gitterpunkte
	int N = LHS.m_Zeilenzahl;
	//array mit der Pivotisierungsmatrix sollte so groß wie N sein,
	//alle Elemente 0
	int *IPIV;
	IPIV = new int[N];
	// ------ RHS oben definiert
	//Spalten von RHS 1 nehmen, um keine c/Fortran Verwirrungen zu provozieren
	int NRHS = 1;
	int LDA = N;
	int LDB = N;
	int INFO;
	//Anzahl sollte die Integer grenzen nicht überschreiten,
	//aber danbn sollte der Aufbau von LHS schon stören
	//int Anzahl=N*N;
	char textflag = 'N'; // "No transpose"; //fürs Rückeinsetzen
	// AUFRUF   A ist LHS.transponiert und B ist RHS
	dgesv_(&N, &NRHS, A.m_Elemente, &LDA, IPIV, RHS.m_Elemente, &LDB, &INFO);
	// ENDE LU ZERLEGUNG
	//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";

	//cout<<"Itmax: "<<Itmax<<"\n";
	double Residual, Residual_1, residual_prev = 0.;
	Residual = 0;
	Residual_1 = 0;
	MPL_Matrix RHS_Teil1;
	RHS_Teil1 = AMF_trans * (S_y * Saeulendichten);

	//Die Schleife ist echt schnell...das sollten höchstens 10 sekunden sein

	for (int Iterationsschritt = 0; Iterationsschritt < Itmax; Iterationsschritt++) {

		//Auf der RHS ändert sich die Apriori-Lösung...d.h. neu berechnen
		//RHS sollte ein Spaltenvektor sein
		RHS = RHS_Teil1 + S_apriori * Dichten_apriori;
		//for(int i=0; i<RHS.m_Zeilenzahl;i++)
		//{
		//    cout<<RHS(i)<<"\n";
		//}
		//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";
		//cout<<"RHS.m_Spaltenzahl: "<<RHS.m_Spaltenzahl<<"\n";
		// Lösungen durch Rückeinsetzen finden
		dgetrs_(&textflag, &N, &NRHS, A.m_Elemente, &LDA, IPIV, RHS.m_Elemente,
				&LDB, &INFO);
		// RHS enthält nun die Lösungen
		Dichten_alt = Dichten;
		Dichten = RHS;
		//for(int i=0; i<RHS.m_Zeilenzahl;i++)
		//{
		//    cout<<RHS(i)<<"\n";
		//}
		MPL_Matrix Differenz = Dichten - Dichten_alt; //Spaltenvektor
		//cout<<"Differenzenvektor\n";
		//cout<<"Dichten.m_Zeilenzahl: "<<Dichten.m_Zeilenzahl<<"\n";
		//cout<<"Dichten_alt.m_Zeilenzahl: "<<Dichten_alt.m_Zeilenzahl<<"\n";
		//cout<<"Dichten.m_Spaltenzahl: "<<Dichten.m_Spaltenzahl<<"\n";
		//cout<<"Dichten_alt.m_Spaltenzahl: "<<Dichten_alt.m_Spaltenzahl<<"\n";
		//skalar als 1x1 Matrix getarnt
		Differenz = Differenz.transponiert() * Differenz;

		double Summe_der_quadratischen_Abweichungen = Differenz(0); //
		//cout<<"Summe_der_quadratischen_Abweichungen: "
		//  <<Summe_der_quadratischen_Abweichungen<<"\n";
		//das ist jetzt sowas wie die Standardabweichung....
		//der Nenner ist nur zur Normierung da
		//(spielt das noch irgendwo ne Rolle...zumindest nicht im Hauptprogramm)
		// TODO herausfinden
		Residual = sqrt(Summe_der_quadratischen_Abweichungen);
		// /sqrt((double)LHS.m_Zeilenzahl-1);
		//cout<<"Iterationsschritt: "<<Iterationsschritt<<"\n";
		//cout<<"Residual: "<<Residual<<"\n";
		if ((Iterationsschritt == 0) || (Iterationsschritt == Itmax - 1)) {
			cerr << "Iterationsschritt:" << Iterationsschritt << "\t"
				 << "Residual: " << Residual << "\n";
		}
		if (Iterationsschritt == 0) {
			// erstes Residuum als ungefähre Fehlerabschätzung
			Residual_1 = Residual;
		}
		//cout<<"Residual_1: "<<Residual_1<<"\n";
		if (Residual < Threshold * Residual_1) {
			//cout<<"Residual: "<<Residual<<"\n";
			//cout<<"Residual_1: "<<Residual_1<<"\n";
			//cout<<"Threshold: "<<Threshold<<"\n";
			//cout<<"Threshold*Residual_1: "<<Threshold*Residual_1<<"\n";
			//evtl konvergenzflag setzen
			cerr << "Konvergenz bei Iterationsschritt: " << Iterationsschritt << "\n";
			break;
			// achtung mit break und continue in for-schleifen
			// (vor allem mit continue->(i++; continue;)
		}
		if (abs(residual_prev - Residual) / Residual < Threshold) {
			cerr << "Konvergenz bei Iterationsschritt: " << Iterationsschritt << endl;
			break;
		}
		residual_prev = Residual;

		Dichten_apriori = Dichten;
		if (Iterationsschritt == (Itmax - 1)) {
			cout << "Residual: " << Residual << "\n";
			cout << "Threshold*Residual_1: " << Threshold *Residual_1 << "\n";
			cout << "Iteration konvergiert nicht, oder zu langsam\n";
			// dynamischen Kram löschen
			Dichten.in_Datei_speichern("/tmp/mlangowski/0/Dichten_nach_Iteration.txt");
			delete[] IPIV;
			return 1; //Fehlerflag 1
		}
	}//for Iterationsschritt
	// keine Probleme während Iteration aufgetreten
	Dichten.in_Datei_speichern("/tmp/mlangowski/0/Dichten_nach_Iteration.txt");
	// dynamischen Kram löschen
	delete[] IPIV;
	return 0;
}


