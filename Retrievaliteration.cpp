/*
 * Retrievaliteration.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 * Copyright (c) 2010-2011 Martin Langowski
 *
 * Initial version created on: 14.09.2010
 *      Author: Martin Langowski
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */

#include"Retrievaliteration.h"

#include"MPL_Matrix.h"
#include "Konfiguration.h"
#include <cmath>

using std::cout;
using std::endl;

extern "C" {
	void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B,
			int *LDB, int *INFO);
	void dgetrs_(char *, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
			double *B, int *LDB, int *INFO);
}

/* Construct the a priori covariance matrix using the a priori values.
 * A common approach seems to be to use 1000%, i.e. set factor to 10. */
MPL_Matrix relative_Sapriori(MPL_Matrix &Dichten_apriori,
		MPL_Matrix &S_apriori, double factor)
{
	MPL_Matrix S_a_i{S_apriori};

	for (int i = 0; i < Dichten_apriori.m_Elementanzahl; ++i) {
		double x_a_i = Dichten_apriori(i);
		if (x_a_i != 0.)
			// relative (co)variance (squared weights)
			S_a_i(i, i) = 1. / (factor*factor * x_a_i*x_a_i);
	}

	return S_a_i;
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
	MPL_Matrix LHS, R, Rl, Ra;
	MPL_Matrix AMF_trans_S_y = AMF.transponiert() * S_y;

	// Solange man die lambdas für die constraints nicht ändern will,
	// sieht die LHS immer gleich aus
	Rl = (Lambda_Breite * (S_Breite.transponiert() * S_Breite));  // Breitenglattung
	Ra = (Lambda_Hoehe * (S_Hoehe.transponiert() * S_Hoehe)); // Hoehenglattung
	R  = Rl + Ra;
	// S_y here is equal to S_y^-1 in ususal retrieval equations,
	// as is S_apriori (~ S_a^-1)
	LHS = (AMF_trans_S_y * AMF);
	LHS += R;
	LHS += S_apriori;
	// LHS += Rl + Ra;
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
	double beta_inv = (Dichten_apriori.transponiert() *
			S_apriori * Dichten_apriori)(0, 0);
	if (beta_inv != 0 && Konf.NO_apriori_scale == -2)
		LHS -= 1 / beta_inv * (S_apriori * Dichten_apriori) *
				(Dichten_apriori.transponiert() * S_apriori);

	////////////////////////////////////////////////////////////////////////////
	// LU Zerlegung der LHS
	////////////////////////////////////////////////////////////////////////////
	MPL_Matrix RHS(Dichten.m_Zeilenzahl, 1, 1.);
	MPL_Matrix Saeulendichten_rest(Saeulendichten.m_Zeilenzahl, 1);
	//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";
	// Zunächst die LU Zerlegung der LHS durchführen mit dummy RHS,
	// in der Iteration dann nur das Rückeinsetzen nutzen
	// START LU ZERLEGUNG
	// VORBEREITEN
	int N = LHS.m_Zeilenzahl;
	//array mit der Pivotisierungsmatrix sollte so groß wie N sein,
	//alle Elemente 0
	int *IPIV = new int[N];
	// ------ RHS oben definiert
	//Spalten von RHS 1 nehmen, um keine c/Fortran Verwirrungen zu provozieren
	int NRHS = 1;
	int LDA = N;
	int LDB = N;
	int INFO;
	//Anzahl sollte die Integer grenzen nicht überschreiten,
	//aber danbn sollte der Aufbau von LHS schon stören
	//int Anzahl=N*N;
	char textflag = 'T'; // "transpose"; //fürs Rückeinsetzen

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
	/* according to [Funke 2005, Steck & von Clarmann 2001],
	 * the iterative solution is
	 * x_{i+1} = x_i + (K^T S_y^{-1} K + R)^{-1}
	 *             \times [K^T S_y^{-1} (y - F(x_i)) + R (x_a - x_i)]
	 * where we have an arbitrary choice for x_0.
	 * In the case of x_0 = x_a, (x_a - x_0) vanishes, but we have
	 * to add x_a to the first solution below.
	 * Instead here we use x_0 = 0 and since F is linear
	 * (F(x) = K*x in our case), F(0) = K*0 = 0.
	 * Therefore, y - F(x_0) = y for calculating x_1.
	 */
	RHS = AMF_trans_S_y * (Saeulendichten);
	if (beta_inv == 0 || Konf.NO_apriori_scale != -2)
		RHS += S_apriori * Dichten_apriori;
#ifdef DEBUG_RETRIEVAL_MATRICES
	LHS.in_Datei_speichern("/tmp/LHS1.dat.gz");
	RHS.in_Datei_speichern("/tmp/RHS1.dat.gz");
	AMF.in_Datei_speichern("/tmp/AMF1.dat.gz");
	Saeulendichten.in_Datei_speichern("/tmp/SDN1.dat.gz");
	Dichten_apriori.in_Datei_speichern("/tmp/DAP1.dat.gz");
	S_y.in_Datei_speichern("/tmp/SY1.dat.gz");
	S_apriori.in_Datei_speichern("/tmp/SAP1.dat.gz");
	R.in_Datei_speichern("/tmp/R1.dat.gz");
#endif /* DEBUG_RETRIEVAL_MATRICES */
	// Lösungen durch LU Zerlegung finden
	// LU Komponenten für später in A abgelegt
	dgesv_(&N, &NRHS, LHS.m_Elemente, &LDA, IPIV, RHS.m_Elemente, &LDB, &INFO);
	//cout<<"RHS.m_Zeilenzahl: "<<RHS.m_Zeilenzahl<<"\n";
	////////////////////////////////////////////////////////////////////////////
	// ENDE LU Zerlegung der LHS
	////////////////////////////////////////////////////////////////////////////
	/* This should be Dichten_apriori + RHS,
	 * if x_0 = x_a was used above. */
	Dichten = RHS;
	///////////////////////////////////////
	// ENDE erster Schritt
	///////////////////////////////////////
	if (beta_inv != 0 && Konf.NO_apriori_scale == -2) {
		double alpha = 1. / beta_inv *
				(Dichten_apriori.transponiert() * S_apriori * Dichten)(0, 0);
		std::cout << "# apriori fit factor = " << alpha << std::endl;
		Dichten_apriori *= alpha;
	}

	///////////////////////////////////////
	// Nachiteration
	///////////////////////////////////////
	for (int Iterationsschritt = 0; Iterationsschritt < Itmax; Iterationsschritt++) {
		Saeulendichten_rest = Saeulendichten - AMF * Dichten;
		MPL_Matrix Dichten_apriori_rest = Dichten_apriori - Dichten;
		Mat_Residual = 0.5 * (Saeulendichten_rest.transponiert() * S_y * Saeulendichten_rest
			+ Dichten_apriori_rest.transponiert() * S_apriori * Dichten_apriori_rest
			+ Dichten.transponiert() * R * Dichten);

		// Anfangsresiduum bestimmen
		Residual = Mat_Residual(0);
		if ((Iterationsschritt == 0) || (Iterationsschritt == Itmax - 1)) {
			cout << "Iterationsschritt:" << Iterationsschritt << "\t"
				 << "Residual: " << Residual << "\n";
		}
		if (Iterationsschritt == 0) {
			// erstes Residuum als ungefähre Fehlerabschätzung
			Residual_1 = Residual;
		}
		cout << "Iterationsschritt:" << Iterationsschritt << "\t"
			 << "residual_prev: " << residual_prev << "\t"
			 << "Residual_1: " << Residual_1 << "\t"
			 << "Residual: " << Residual << "\n";
		if (Residual < Threshold * Residual_1) {
			cout << "Konvergenz 1 bei Iterationsschritt: " << Iterationsschritt << endl;
			cout << "Residual: " << Residual << endl;
			break;
		}
		if (abs(residual_prev - Residual) < Threshold) {
			cout << "Konvergenz 2 bei Iterationsschritt: " << Iterationsschritt << endl;
			cout << "Residual: " << Residual << endl;
			break;
		}
		residual_prev = Residual;
		//RHS sollte ein Spaltenvektor sein
		RHS = AMF_trans_S_y * Saeulendichten_rest
			  + S_apriori * Dichten_apriori_rest
			  - R * Dichten;
		// Lösungen durch Rückeinsetzen finden
		dgetrs_(&textflag, &N, &NRHS, LHS.m_Elemente, &LDA, IPIV, RHS.m_Elemente,
				&LDB, &INFO);
		Dichten += RHS; // inkrementierung
	}//for Iterationsschritt
	// keine Probleme während Iteration aufgetreten
	//Dichten.in_Datei_speichern("/tmp/mlangowski/0/Dichten_nach_Iteration.txt");
	// dynamischen Kram löschen
	delete[] IPIV;
	S_apriori += R; // update to speed up AKM calculation
	return 0;
}

#ifdef HAVE_EIGEN3
#include <Eigen/Dense>
/*
 * Retrieval variant using Eigen3 (http://eigen.tuxfamily.org).
 * This provides an easy interface to matrix operations and should also be
 * robust and hopefully a little more modern compared to the self-made
 * MPL_Matrix class.
 */
int Retrievaliteration_Eigen(MPL_Matrix &Dichten,
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
	using Eigen::Dynamic;
	using Eigen::LDLT;
	using Eigen::Map;
	using Eigen::Matrix;
	using Eigen::MatrixXd;
	using Eigen::PartialPivLU;
	using Eigen::RowMajor;

	bool fit_apriori = (Konf.NO_apriori_scale == -2);
	int grad_add_rows = 1 ? fit_apriori : 0;
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		S_alt{S_Hoehe.m_Elemente, S_Hoehe.m_Zeilenzahl, S_Hoehe.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		S_lat{S_Breite.m_Elemente, S_Breite.m_Zeilenzahl, S_Breite.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		K_m{AMF.m_Elemente, AMF.m_Zeilenzahl, AMF.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		Sy_m{S_y.m_Elemente, S_y.m_Zeilenzahl, S_y.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		Sa_m{S_apriori.m_Elemente, S_apriori.m_Zeilenzahl, S_apriori.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		x_m{Dichten.m_Elemente, Dichten.m_Zeilenzahl, Dichten.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		xa_m{Dichten_apriori.m_Elemente, Dichten_apriori.m_Zeilenzahl, Dichten_apriori.m_Spaltenzahl};
	Map<Matrix<double, Dynamic, Dynamic, RowMajor> >
		y_m{Saeulendichten.m_Elemente, Saeulendichten.m_Zeilenzahl, Saeulendichten.m_Spaltenzahl};

	MatrixXd R_m{Lambda_Hoehe * S_alt.transpose() * S_alt +
			Lambda_Breite * S_lat.transpose() * S_lat};
	MatrixXd K_T_Sy{K_m.transpose() * Sy_m};
	MatrixXd LHS_m{K_T_Sy * K_m + Sa_m + R_m};
	MatrixXd RHS_m{K_T_Sy * y_m + Sa_m * xa_m};

	Matrix<double, Dynamic, Dynamic, RowMajor>
		nHess(LHS_m.rows() + grad_add_rows, LHS_m.cols() + grad_add_rows);
	Matrix<double, Dynamic, Dynamic, RowMajor>
		gradf(RHS_m.rows() + grad_add_rows, RHS_m.cols());
	Matrix<double, Dynamic, Dynamic, RowMajor> dx, nlp_m;
	if (fit_apriori) {
		nHess << LHS_m, -1. * Sa_m * xa_m,
			  -1. * xa_m.transpose() * Sa_m, xa_m.transpose() * Sa_m * xa_m;
	} else
		nHess = LHS_m;
	//PartialPivLU<MatrixXd> nHess_fact;
	LDLT<MatrixXd> nHess_fact; // LDLT gives slightly higher accuracy
	nHess_fact.compute(nHess);
	double alpha = 1.;
	double nlp_old = 1e25, nlp = 1e24;
	int it_cnt = 0;
	while (std::abs(nlp_old - nlp) > Konf.m_Convergence_Treshold
			&& it_cnt < Konf.m_Max_Zahl_Iterationen) {
		nlp_old = nlp;
		if (fit_apriori)
			gradf << K_T_Sy * (y_m - K_m * x_m) - Sa_m * (x_m - alpha * xa_m) - R_m * x_m,
					xa_m.transpose() * Sa_m * (x_m - alpha * xa_m);
		else
			gradf = K_T_Sy * (y_m - K_m * x_m) - Sa_m * (x_m - alpha * xa_m) - R_m * x_m;
		dx = nHess_fact.solve(gradf);
		x_m += dx.topRows(dx.rows() - grad_add_rows);
		if (fit_apriori) {
			double dalpha = dx.bottomRows(grad_add_rows)(0, 0);
			alpha += dalpha;
			std::cout << "# apriori fit factor adj = " << dalpha << std::endl;
			std::cout << "# apriori fit factor = " << alpha << std::endl;
		}
		// negative log posterior probability
		nlp_m = 0.5 * ((y_m - K_m * x_m).transpose() * Sy_m * (y_m - K_m * x_m)
				+ (x_m - alpha * xa_m).transpose() * Sa_m * (x_m - alpha * xa_m)
				+ x_m.transpose() * R_m * x_m);
		nlp = nlp_m(0, 0);
		std::cout << "iteration nr. " << it_cnt;
		std::cout << " negative log posterior: " << nlp << std::endl;
		++it_cnt;
	}
	Sa_m += R_m; // update to speed up AKM calculation
	return 0;
}
#endif /* HAVE_EIGEN3 */


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
			cout << "Iterationsschritt:" << Iterationsschritt << "\t"
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
			cout << "Konvergenz bei Iterationsschritt: " << Iterationsschritt << "\n";
			cout << "Residual: " << Residual << endl;
			break;
			// achtung mit break und continue in for-schleifen
			// (vor allem mit continue->(i++; continue;)
		}
		if (abs(residual_prev - Residual) / Residual < Threshold) {
			cout << "Konvergenz bei Iterationsschritt: " << Iterationsschritt << endl;
			cout << "Residual: " << Residual << endl;
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


