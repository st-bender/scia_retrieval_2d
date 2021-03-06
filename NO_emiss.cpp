/*
 * NO_emiss.cpp
 *
 * Copyright (c) 2011-2017 Stefan Bender
 *
 * Initial version created on: 20.04.2011
 *      Author: Stefan Bender
 *
 * This file is part of scia_retrieval_2d
 *
 * scia_retrieval_2d is free software: you can redistribute it or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 2.
 * See accompanying COPYING.GPL2 file or http://www.gnu.org/licenses/gpl-2.0.html.
 */
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

#include "constants.h"
#include "NO_emiss.h"
#include "Sonnenspektrum.h"
#include "Messung_Limb.h"
#include "Glaetten.h"

double B_v(int v)
{
	return 1.7049
		- 0.01754 * (v + 0.5)
		- 1.304e-5 * (v + 0.5) * (v + 0.5)
		- 2.48e-7 * (v + 0.5) * (v + 0.5) * (v + 0.5)
		- 3.86e-8 * (v + 0.5) * (v + 0.5) * (v + 0.5) * (v + 0.5);
}
double D_v(int v)
{
	return 5.4542e-6
		+ 1.515e-8 * (v + 0.5)
		+ 6.99e-10 * (v + 0.5) * (v + 0.5);
}
double A_v(int v)
{
	return 123.258
		- 0.2363 * (v + 0.5)
		- 3.582e-3 * (v + 0.5) * (v + 0.5)
		- 2.938e-4 * (v + 0.5) * (v + 0.5) * (v + 0.5);
}
double Evibl(int v)
{
	return NO_const::w_l * (v + 0.5)
		- NO_const::w_Xl * (v + 0.5) * (v + 0.5)
		+ NO_const::w_Yl * (v + 0.5) * (v + 0.5) * (v + 0.5)
		+ NO_const::w_Zl * (v + 0.5) * (v + 0.5) * (v + 0.5) * (v + 0.5);
}
double Evibu(int v)
{
	return NO_const::w_u * (v + 0.5)
		- NO_const::w_Xu * (v + 0.5) * (v + 0.5)
		+ NO_const::w_Yu * (v + 0.5) * (v + 0.5) * (v + 0.5);
}
double F0(double B, double D, double Y, double lam, double j)
{
	return B * ((j + 0.5) * (j + 0.5) - lam * lam
			- 0.5 * std::sqrt(4. * (j + 0.5) * (j + 0.5)
				+ Y * (Y - 4.) * lam * lam)) - D * j * j * j * j;
}
double F1(double B, double D, double Y, double lam, double j)
{
	return B * ((j + 0.5) * (j + 0.5) - lam * lam
			+ 0.5 * std::sqrt(4. * (j + 0.5) * (j + 0.5)
				+ Y * (Y - 4.) * lam * lam))
		- D * (j + 1) * (j + 1) * (j + 1) * (j + 1);
}

// default constructor
NO_emiss::NO_emiss(int vu, int vl, int vl_abs, double T) :
	v_u(vu), v_l(vl), v_l_abs(vl_abs), NJ(40), Temp(T)
{
	alloc_memory();
	set_constants();
	populate_Fs();
	calc_lines_emiss_absorp();
	set_Hoenl_London_abs();
	set_Hoenl_London_emiss();
}

void NO_emiss::alloc_memory()
{
	// allocate
	MPL_Matrix Fl(2, NJ + 3);
	MPL_Matrix Fl_abs(2, NJ + 3);
	MPL_Matrix Fu(2, NJ + 3);
	MPL_Matrix NjtoN(2, NJ + 1);

	MPL_Matrix xl_K(12, NJ + 1);
	MPL_Matrix xl_K_abs(12, NJ + 1);
	MPL_Matrix lam_K(12, NJ + 1);
	MPL_Matrix lam_K_abs(12, NJ + 1);

	MPL_Matrix vfHL_j(12, NJ + 1);
	MPL_Matrix vfHL_k(12, NJ + 1);

	MPL_Matrix sol(12, NJ + 1);

	MPL_Matrix fos(NO_const::l_vu + 1, NO_const::l_vl + 1);
	MPL_Matrix ffc(NO_const::l_vu + 1, NO_const::l_vl + 1);
	MPL_Matrix fa(NO_const::l_vu + 1, NO_const::l_vl + 1);
	MPL_Matrix flam(NO_const::l_vu + 1, NO_const::l_vl + 1);

	MPL_Matrix exc(2, NJ + 1);
	MPL_Matrix exc_pqr(3, NJ + 1);

	MPL_Matrix vfHL_emiss(12, NJ + 1);
	MPL_Matrix vfHL_emiss_K(12, NJ + 1);
	MPL_Matrix quantjup(12, NJ + 1);

	MPL_Matrix gam_j(12, NJ + 1);

	// copy
	F_l = Fl;
	F_l_abs = Fl_abs;
	F_u = Fu;
	NJ_to_N = NjtoN;

	xlines_K = xl_K;
	xlines_K_abs = xl_K_abs;
	lambda_K = lam_K;
	lambda_K_abs = lam_K_abs;

	vf_HL_J = vfHL_j;
	vf_HL_K = vfHL_k;

	solar = sol;

	f_osc = fos;
	f_FC = ffc;
	f_A = fa;
	f_lam = flam;

	excit = exc;
	excit_pqr = exc_pqr;

	vf_HL_emiss = vfHL_emiss;
	vf_HL_emiss_K = vfHL_emiss_K;
	quant_j_up = quantjup;

	gamma_j = gam_j;
}

void NO_emiss::set_constants()
{
	// rotational parameters of the upper state, from Eparvier and Barth
	B_Vu = 1.9965 - 0.01915 * (v_u + 0.5);
	D_Vu = 5.4e-6 + 3.4056e-8 * (v_u + 0.5);
	gam_u = -2.765e-3;
	// rotational parameters of the lower state, from Eparvier and Barth
	// lower state of the emission
	B_Vl = B_v(v_l);
	D_Vl = D_v(v_l);
	A_Vl = A_v(v_l);
	Y_l = A_Vl / B_Vl;
	// lower state of the absorption
	B_Vl_abs = B_v(v_l_abs);
	D_Vl_abs = D_v(v_l_abs);
	A_Vl_abs = A_v(v_l_abs);
	Y_l_abs = A_Vl_abs / B_Vl_abs;
	lambda_l = 1.0;
	// lower state mean energy, no fine-splitting
	E_l = A_Vl * lambda_l * 0.5;
	E_l_abs = A_Vl_abs * lambda_l * 0.5;
	// total electric energy
	E_tot = NO_const::E_u - E_l;
	E_tot_abs = NO_const::E_u - E_l_abs;
	// vibrational term energy
	// lower state of the emission
	E_vib_l = Evibl(v_l);
	//lower state of the absorption
	E_vib_l_abs = Evibl(v_l_abs);

	E_vib_u = Evibu(v_u);

	E_vib = E_vib_u - E_vib_l;
	W_vib = E_tot + E_vib;

	E_vib_abs = E_vib_u - E_vib_l_abs;
	W_vib_abs = E_tot_abs + E_vib_abs;

	phi = 0.5; // degeneracy of rotational states
	f_boltz = 100. / (phys::k2 * Temp);
	part_el = 1.0 + std::exp(-f_boltz * NO_const::E_u);

	// Rotational and total partitioning sum
	sum_j = 0.;
	for (int i = 0; i <= NJ + 2; i++) {
		double j = i + 0.5;
		sum_j += (2. * j + 1.) * std::exp(-f_boltz * B_Vl_abs * j * (j + 0.5));
	}
	// total partitioning sum considering vibrational ground state
	// and Lambda splitting (Hund's case a)
	sum_j *= part_el * (std::exp(-f_boltz * E_l_abs)
			+ std::exp(+f_boltz * E_l_abs));
}
void NO_emiss::populate_Fs()
{
	for (int i = 0; i <= NJ + 2; i++) {
		double j = i + 0.5;
		double flow = B_Vl_abs * j * (j + 0.5);
		F_l(0, i) = F0(B_Vl, D_Vl, Y_l, lambda_l, j);
		F_l(1, i) = F1(B_Vl, D_Vl, Y_l, lambda_l, j);
		F_l_abs(0, i) = F0(B_Vl_abs, D_Vl_abs, Y_l_abs, lambda_l, j);
		F_l_abs(1, i) = F1(B_Vl_abs, D_Vl_abs, Y_l_abs, lambda_l, j);

		F_u(0, i) = B_Vu * i * (i + 1)
			- D_Vu * i * i * (i + 1) * (i + 1)
			+ 0.5 * gam_u * i;
		F_u(1, i) = B_Vu * i * (i + 1)
			- D_Vu * i * i * (i + 1) * (i + 1)
			- 0.5 * gam_u * (i + 1);

		if (i <= NJ) {
			NJ_to_N(0, i) = phi * (2. * j + 1.) / sum_j
				* std::exp(-f_boltz * (flow - E_l_abs));
			if (i > 0 && i < NJ)
				NJ_to_N(1, i + 1) = phi * (2. * j + 1.) / sum_j
					* std::exp(-f_boltz * (flow + E_l_abs));
		}
	}
}

void NO_emiss::calc_lines_emiss_absorp()
{
	int i, j, l;
	double E_rot, E_rot_abs;

	for (i = 0; i <= NJ; i++) {
		double k_l = i;
		double j_l = 0., j_u = 0., k_u = 0.;
		int i_l, i_u;
		for (l = 0; l <= 11; l++) {
			if (l == 0 || l == 1 || l == 2 || l == 6 || l == 8 || l == 10)
				j_l = k_l + 0.5;
			else
				j_l = k_l - 0.5;

			if (l == 0 || l == 3 || l == 6 || l == 9) j_u = j_l - 1.;
			if (l == 1 || l == 4 || l == 7 || l == 10) j_u = j_l;
			if (l == 2 || l == 5 || l == 8 || l == 11) j_u = j_l + 1.;

			if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
				k_u = j_u - 0.5;
			else
				k_u = j_u + 0.5;

			if (k_u >= 0 && j_u >= 0.5 && j_l >= 0.5) {
				if (k_u - j_u + 0.5 == 0.) i_u = 0;
				else i_u = 1;
				if (k_l - j_l + 0.5 == 0.) i_l = 0;
				else i_l = 1;

				if (k_u >= 0 && j_l >= 0.5) {
					E_rot = F_u(i_u, k_u) - F_l(i_l, j_l - 0.5);
					E_rot_abs = F_u(i_u, k_u) - F_l_abs(i_l, j_l - 0.5);
				} else
					E_rot = E_rot_abs = 0.;

				if (E_rot_abs != 0.) {
					if (k_l > 1 || (k_l == 1 && j_l == 1.5) || (k_l == 0)) {
						xlines_K_abs(l, k_l) = E_tot_abs + E_vib_abs + E_rot_abs;
						lambda_K_abs(l, k_l) = 1.e8 / xlines_K_abs(l, k_l);
					}
				}
				if (E_rot != 0.) {
					if ((k_l > 1 || (k_l == 1 && j_l == 1.5)
							|| k_l == 0) && k_u <= NJ) {
						xlines_K(l, k_u) = E_tot + E_vib + E_rot;
						lambda_K(l, k_u) = 1.e8 / xlines_K(l, k_u);
					}
				}
			}
		}
		if (i == 0) quant_K_vec.push_back(0.5);
		else if (i == 1) quant_K_vec.push_back(1.5);
		else quant_K_vec.push_back(k_l);
	}

	// fix up the rest
	for (i = 0; i <= NJ; i++)
		for (j = 0; j <= 11; j++)
			if (xlines_K(j, i) != 0.)
				lambda_K(j, i) = 1.e8 / xlines_K(j, i);
}

// Absorption Hoenl-London factors from ground state
// from Earls, 1935, normalised to 4*(2J+1)
void NO_emiss::set_Hoenl_London_abs()
{
	int i;
	double j, u, d1, d2, d3;

	// j = 0.5
	// P1
	vf_HL_J(0, 0) = vf_HL_K(0, 0) = 4. / 6.;
	// Q1
	vf_HL_J(1, 0) = vf_HL_K(1, 0) = 4. / 3.;
	// R1
	vf_HL_J(2, 0) = vf_HL_K(2, 0) = 0.;
	// P2
	vf_HL_J(3, 0) = vf_HL_K(3, 1) = 0.;
	// Q2
	vf_HL_J(4, 0) = vf_HL_K(4, 1) = 0.;
	// R2
	vf_HL_J(5, 0) = vf_HL_K(5, 1) = 0.;
	// qP21
	vf_HL_J(6, 0) = vf_HL_K(6, 0) = 0.;
	// pQ12
	vf_HL_J(7, 0) = vf_HL_K(7, 1) = 4. / 3.;
	// sR21
	vf_HL_J(8, 0) = vf_HL_K(8, 0) = 0.;
	// oP12
	vf_HL_J(9, 0) = vf_HL_K(9, 1) = 4. / 6.;
	// rQ21
	vf_HL_J(10, 0) = vf_HL_K(10, 0) = 0.;
	// qR12
	vf_HL_J(11, 0) = vf_HL_K(11, 1) = 0.;

	// j > 0.5
	for (i = 1; i <= NJ; i++) {
		j = i + 0.5;
		u = 1. / std::sqrt(Y_l * Y_l - 4. * Y_l
				+ (2. * j + 1.) * (2. * j + 1.));
		d1 = 8. * j; // 32 * j / 4
		d2 = 8. * (j + 1.); // 32 * (j + 1) / 4
		d3 = 8. * j * (j + 1.); // 32 * j * (j + 1) / 4

		// P1
		vf_HL_J(0, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		vf_HL_K(0, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		// Q1
		vf_HL_J(1, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				+ u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		vf_HL_K(1, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				+ u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		// R1
		vf_HL_J(2, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		vf_HL_K(2, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		// P2
		vf_HL_J(3, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d2;
		if (i != NJ)
			vf_HL_K(3, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d2;
		// Q2
		vf_HL_J(4, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				+ u * (8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l)) / d3;
		if (i != NJ)
			vf_HL_K(4, i + 1) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
					+ u * (8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l)) / d3;
		// R2
		vf_HL_J(5, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
		if (i != NJ)
			vf_HL_K(5, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
		// qP21
		vf_HL_J(6, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		vf_HL_K(6, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		// pQ12
		vf_HL_J(7, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				- u * (8.*j*j*j + 12.*j*j - 2.*j + 1. - 2*Y_l)) / d3;
		if (i != NJ)
			vf_HL_K(7, i + 1) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
					- u * (8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l)) / d3;
		// sR21
		vf_HL_J(8, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		vf_HL_K(8, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		// oP12
		vf_HL_J(9, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d2;
		if (i != NJ)
			vf_HL_K(9, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d2;
		// rQ21
		vf_HL_J(10, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				- u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		vf_HL_K(10, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				- u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		// qR12
		vf_HL_J(11, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
		if (i != NJ)
			vf_HL_K(11, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
	}
}
// Hoenl-London factors of emission, relative to upper state K
// and upper state J ? anyways, these are normalised to 2*(2J+1)
void NO_emiss::set_Hoenl_London_emiss()
{
	int i;
	double j, u, d1, d2, d3;

	// j' = 0.5
	// Q11
	vf_HL_emiss(1, 0) = 2. / 3.;
	vf_HL_emiss_K(1, 0) = 2. / 3.;
	quant_j_up(1, 0) = 0.5;
	// R11
	vf_HL_emiss(2, 0) = 2. / 6.;
	vf_HL_emiss_K(2, 1) = 2. / 6.;
	quant_j_up(2, 1) = 0.5;
	// sR21
	vf_HL_emiss(8, 0) = 2. / 6.;
	vf_HL_emiss_K(8, 2) = 2. / 6.;
	quant_j_up(8, 2) = 0.5;
	// rQ21
	vf_HL_emiss(10, 0) = 2. / 3.;
	vf_HL_emiss_K(10, 1) = 2. / 3.;
	quant_j_up(10, 1) = 0.5;

	// j' > 0.5
	for (i = 1; i <= NJ; i++) {
		j = i + 0.5;
		u = 1. / std::sqrt(Y_l * Y_l - 4. * Y_l
				+ (2. * j + 1.) * (2. * j + 1.));
		d1 = 16. * j; // 32 * j / 2
		d2 = 16. * (j + 1.); // 32 * (j + 1) / 2
		d3 = 16. * j * (j + 1.); // 32 * j * (j + 1) / 2

		// P1
		vf_HL_emiss(0, i) = ((2.*j + 1.)*(2.*j + 1.)
				+ (2.*j + 1.)*u*(4.*j*j + 4.*j + 1. - 2.*Y_l))/d1;
		if (i > 0) {
			vf_HL_emiss_K(0, i - 1) = ((2.*j + 1.)*(2.*j + 1.)
					+ (2.*j + 1.)*u*(4.*j*j + 4.*j + 1. - 2.*Y_l))/d1;
			quant_j_up(0, i - 1) = j;
		}
		// Q1
		vf_HL_emiss(1, i) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
				+ u*(8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l))/d3;
		vf_HL_emiss_K(1, i) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
				+ u*(8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l))/d3;
		quant_j_up(1, i) = j;
		// R1
		vf_HL_emiss(2, i) = ((2.*j + 1.)*(2.*j + 1.)
				+ (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d2;
		if (i != NJ) {
			vf_HL_emiss_K(2, i + 1) = ((2.*j + 1.)*(2.*j + 1.)
					+ (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d2;
			quant_j_up(2, i + 1) = j;
		}
		// P2
		vf_HL_emiss(3, i) = ((2.*j + 1.)*(2.*j + 1.)
				+ (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d1;
		vf_HL_emiss_K(3, i) = ((2.*j + 1.)*(2.*j + 1.)
				+ (2.*j + 1.)*u*(4.*j*j + 4.*j - 7 + 2.*Y_l))/d1;
		quant_j_up(3, i) = j;
		// Q2
		vf_HL_emiss(4, i) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
				+ u*(8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l))/d3;
		if (i != NJ) {
			vf_HL_emiss_K(4, i + 1) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
					+ u*(8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l))/d3;
			quant_j_up(4, i + 1) = j;
		}
		// R2
		vf_HL_emiss(5, i) = ((2.*j + 1.)*(2.*j + 1.)
				+ (2.*j + 1.)*u*(4.*j*j + 4.*j + 1 - 2.*Y_l))/d2;
		if (i < NJ - 1) {
			vf_HL_emiss_K(5, i + 2) = ((2.*j + 1.)*(2.*j + 1.)
					+ (2.*j + 1.)*u*(4.*j*j + 4.*j + 1 - 2.*Y_l))/d2;
			quant_j_up(5, i + 2) = j;
		}
		// qP21
		vf_HL_emiss(6, i) = ((2.*j + 1.)*(2.*j + 1.)
				- (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d1;
		vf_HL_emiss_K(6, i) = ((2.*j + 1.)*(2.*j + 1.)
				- (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d1;
		quant_j_up(6, i) = j;
		// pQ12
		vf_HL_emiss(7, i) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
				- u*(8.*j*j*j + 12.*j*j - 2.*j - 7. + 2*Y_l))/d3;
		vf_HL_emiss_K(7, i) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
				- u*(8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l))/d3;
		quant_j_up(7, i) = j;
		// sR21
		vf_HL_emiss(8, i) = ((2.*j + 1.)*(2.*j + 1.)
				- (2.*j + 1.)*u*(4.*j*j + 4.*j + 1. - 2.*Y_l))/d2;
		if (i < NJ - 1) {
			vf_HL_emiss_K(8, i + 2) = ((2.*j + 1.)*(2.*j + 1.)
					- (2.*j + 1.)*u*(4.*j*j + 4.*j + 1. - 2.*Y_l))/d2;
			quant_j_up(8, i + 2) = j;
		}
		// oP12
		vf_HL_emiss(9, i) = ((2.*j + 1.)*(2.*j + 1.)
				- (2.*j + 1.)*u*(4.*j*j + 4.*j + 1. - 2.*Y_l))/d1;
		if (i > 0) {
			vf_HL_emiss_K(9, i - 1) = ((2.*j + 1.)*(2.*j + 1.)
					- (2.*j + 1.)*u*(4.*j*j + 4.*j + 1. - 2.*Y_l))/d1;
			quant_j_up(8, i - 1) = j;
		}
		// rQ21
		vf_HL_emiss(10, i) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
				- u*(8.*j*j*j + 12.*j*j - 2.*j + 1 - 2.*Y_l))/d3;
		if (i != NJ) {
			vf_HL_emiss_K(10, i + 1) = (2.*j + 1.)*((4.*j*j + 4.*j - 1.)
					- u*(8.*j*j*j + 12.*j*j - 2.*j + 1 - 2.*Y_l))/d3;
			quant_j_up(10, i + 1) = j;
		}
		// qR12
		vf_HL_emiss(11, i) = ((2.*j + 1.)*(2.*j + 1.)
				- (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d2;
		if (i != NJ) {
			vf_HL_emiss_K(11, i + 1) = ((2.*j + 1.)*(2.*j + 1.)
					- (2.*j + 1.)*u*(4.*j*j + 4.*j - 7. + 2.*Y_l))/d2;
			quant_j_up(11, i + 1) = j;
		}
	}
}

// Get solar spectral data
void NO_emiss::get_solar_data(Sonnenspektrum &sol_spec)
{
	int i, j;

	for (i = 0; i <= NJ; i++)
		for (j = 0; j <= 11; j++)
			if (lambda_K_abs(j, i) != 0.) {
				// translate from Å to nm and from ph/s/cm^2/nm to ph/s/cm^2/Å
				solar(j, i) = 0.1 * interpolate(sol_spec.m_Wellenlaengen,
						sol_spec.m_Intensitaeten, 0.1 * lambda_K_abs(j, i));
			}
}

// Get data from Luque et al.
int NO_emiss::read_luque_data_from_file(std::string filename)
{
	int i, j;
	int vu, vl;
	double d_dum;
	std::string s_dum;
	std::ifstream lfile(filename.c_str());

	if (!(lfile.is_open())) {
		std::cerr << "Failed to open Luque et al. data from `"
			 << filename.c_str() << "'." << std::endl;
		return 1;
	}

	// skip the first three lines
	for (i = 0; i < 3; i++) std::getline(lfile, s_dum);

	for (i = 0; i <= NO_const::l_vu; i++)
		for (j = 0; j <= NO_const::l_vl; j++) {
			std::stringstream ss;
			std::getline(lfile, s_dum);
			ss << s_dum;
			ss >> vu >> vl >> f_FC(i, j) >> d_dum
				  >> f_lam(i, j) >> f_osc(i, j) >> f_A(i, j);
		}

	return 0;
}

// the excitation of J': sum over all J
// Excitation of upper state J-levels: k_u = j_u - 0.5, k_u = j_u + 0.5
void NO_emiss::calc_excitation()
{
	int i, l, k_l, k_u, nidx;
	double j_l = 0., j_u = 0., nj_frac, sum1, sum2;
	double f_FC_tot;

	for (i = 0; i <= NJ; i++) {
		j_u = i + 0.5;
		sum1 = sum2 = 0.0;
		for (l = 0; l <= 11; l++) {
			if (l == 0 || l == 3 || l == 6 || l == 9) j_l = j_u + 1.;
			if (l == 1 || l == 4 || l == 7 || l == 10) j_l = j_u;
			if (l == 2 || l == 5 || l == 8 || l == 11) j_l = j_u - 1.;

			if (l == 0 || l == 1 || l == 2 || l == 6 || l == 8 || l == 10) {
				k_l = j_l - 0.5;
				nidx = 0;
			} else {
				k_l = j_l + 0.5;
				nidx = 1;
			}
			if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
				k_u = j_u - 0.5;
			else
				k_u = j_u + 0.5;

			if ((k_l >= 0) && (k_u >= 0) && (j_l >= 0.5) && (k_l <= NJ)) {
				if ((k_l >= 1) || (k_l == 1 && j_l == 1.5) || (k_l == 0)) {
					nj_frac = NJ_to_N(nidx, k_l);
					if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
						sum1 += lambda_K_abs(l, k_l) * lambda_K_abs(l, k_l)
							* solar(l, k_l) * vf_HL_K(l, k_l) / (2. * j_l + 1.)
							* nj_frac;
					else
						sum2 += lambda_K_abs(l, k_l) * lambda_K_abs(l, k_l)
							* solar(l, k_l) * vf_HL_K(l, k_l) / (2. * j_l + 1.)
							* nj_frac;
					excit_pqr(l % 3, i) += vf_HL_K(l, k_l) / (2. * j_l + 1.)
						* nj_frac;
				}
			}
		}
		sum1 *= phys::flux * f_osc(v_u, v_l_abs);
		sum2 *= phys::flux * f_osc(v_u, v_l_abs);
		excit(0, i) = sum1;
		if (i < NJ) excit(1, i + 1) = sum2;
	}

	f_FC_tot = 0.;
	for (i = 0; i <= NO_const::l_vl; i++) {
		double l_vib = 1.e8 / f_lam(v_u, i);
		f_FC_tot += f_FC(v_u, i) * l_vib*l_vib*l_vib;
	}
	W_vib = 1.e8 / f_lam(v_u, v_l);
	f_FC_v = f_FC(v_u, v_l) * W_vib*W_vib*W_vib / f_FC_tot;
}

void NO_emiss::calc_line_emissivities()
{
	int i, l, k_l, k_u;
	double j_l = 0., j_u = 0.;

	emiss_tot = 0.;
	for (i = 0; i <= NJ; i++) {
		k_u = i;
		for (l = 0; l <= 11; l++) {
			if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
				j_u = k_u + 0.5;
			else
				j_u = k_u - 0.5;

			if (l == 0 || l == 3 || l == 6 || l == 9) j_l = j_u + 1.;
			if (l == 1 || l == 4 || l == 7 || l == 10) j_l = j_u;
			if (l == 2 || l == 5 || l == 8 || l == 11) j_l = j_u - 1.;

			if (l == 0 || l == 1 || l == 2 || l == 6 || l == 8 || l == 10)
				k_l = j_l - 0.5;
			else
				k_l = j_l + 0.5;

			if ((k_l == 0) || (k_l == 1 && j_l == 1.5) || (k_l > 1)) {
				if (j_l >= 0.5 && j_u >= 0.5) {
					if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
						gamma_j(l, k_u) =
							f_FC_v * vf_HL_emiss_K(l, k_u) / (2.*j_u + 1.)
							* excit(0, j_u - 0.5);
					else
						gamma_j(l, k_u) =
							f_FC_v * vf_HL_emiss_K(l, k_u) / (2.*j_u + 1.)
							* excit(1, j_u + 0.5);
				}
				emiss_tot += gamma_j(l, k_u);
			}
		}
	}
}

/* line polarisations as per
 * Zare, J. Chem. Phys. 45, 4510 (1966)
 * doi:10.1063/1.1727531
 */
std::vector<double> NO_emiss::calc_polarisation()
{
	double j_u = 0.;
	std::vector<double> pol_j(12, 0.);

	for (int i = 0; i <= NJ; i++) {
		double excit_p = excit_pqr(0, i);
		double excit_q = excit_pqr(1, i);
		double excit_r = excit_pqr(2, i);
		for (size_t j = 0; j < 12; j++) {
			if (j == 0 || j == 1 || j == 2 || j == 7 || j == 9 || j == 11)
				j_u = i + 0.5;
			else
				j_u = i - 0.5;

			switch (j % 3) {
			case 0: /* P branch */
				pol_j.at(j) += - (2.*j_u - 1.) / (6.*j_u + 7.) * excit_q
							+ (2.*j_u*j_u - j_u) / (14.*j_u*j_u + 33.*j_u + 20.) * excit_p
							+ 1./7. * excit_r;
				break;
			case 1: /* Q branch */
				pol_j.at(j) +=  (4.*(j_u*j_u + j_u) - 3.) / (8.*(j_u*j_u + j_u) - 1.) * excit_q
							- (2.*j_u - 1.) / (6.*j_u + 7.) * excit_p
							- (2.*j_u + 3.) / (6.*j_u - 1.) * excit_r;
				break;
			case 2: /* R branch */
				pol_j.at(j) += - (2.*j_u + 3.) / (6.*j_u - 1.) * excit_q
							+ (2.*j_u*j_u + 5.*j_u + 3.) / (14.*j_u*j_u - 5.*j_u + 1.) * excit_r
							+ 1./7. * excit_p;
				break;
			}
		}
	}

	return pol_j;
}

/* SCIAMACHY polarisation correction from
 * "SCIAMACHY Level 0 to 1c Processing
 *  Algorithm Theoretical Basis Document"
 * Doc. No. ENV-ATB-DLR-SCIA-0041, pages 96--110
 * and from Liebing et al, Atmos. Meas. Tech. 6, 1503--1520, 2013
 *          doi:10.5194/amt-6-1503-2013
 * and from M.L.'s paper in preparation.
 */
void NO_emiss::pol_corr(double SZA, double rel_SAA,
		double mu2, double mu3)
{
	const double deg = M_PI / 180.0;
	double sin_th0 = std::sin(SZA * deg);
	double cos_phi = std::cos(rel_SAA * deg);
	double sin_phi = std::sin(rel_SAA * deg);
	double cos_Theta = sin_th0 * cos_phi;
	double cos_chi = sin_th0 * sin_phi / std::sqrt(1. - cos_Theta*cos_Theta);
	double chi;

	std::vector<double> Ps = calc_polarisation();

	if (SZA < 90.)
		chi = std::acos(cos_chi);
	else
		chi = std::acos(-cos_chi);

	emiss_tot = 0.;

	for (int i = 0; i < 12; i++) {
		double Q = Ps.at(i) * std::cos(2. * chi);
		double U = Ps.at(i) * std::sin(2. * chi);
		double f = 1. + mu2 * Q + mu3 * U;
		for (int j = 0; j <= NJ; j++) {
			gamma_j(i, j) *= f;
			emiss_tot += gamma_j(i, j);
		}
	}
}

// integral of the slit function, using Bronstein formula 87.
double int_slit_func(double fwhm, double x)
{
	const double cnorm_inv = fwhm*fwhm*fwhm * 0.25 * M_1_PI * M_SQRT1_2;
	const double a = 0.5 * fwhm;
	const double f1 = 0.25 / (M_SQRT2 * a*a*a);
	const double num = x*x + M_SQRT2*a*x + a*a;
	const double den = num - 2. * M_SQRT2*a*x;
	const double y = M_SQRT2 * x/a;

	return cnorm_inv * f1 * (std::log(num / den)
		+ 2. * (std::atan(y + 1.) + std::atan(y - 1.)));
}

void NO_emiss::scia_convolve(Messung &ml)
{
	int i, j;
	std::vector<double> x = ml.m_Wellenlaengen;
	std::vector<double>::iterator x_it;
	std::vector<double>::iterator spec_max;
	spec_scia_res.resize(x.size());

	for (i = 0; i <= NJ; i++) {
		for (j = 0; j < 12; j++) {
			double NO_wl = get_lambda_K(j, i);
			for (x_it = x.begin(); x_it != x.end() - 1; ++x_it) {
				double dl_i = std::abs(NO_wl - *x_it);
				if (dl_i < 1.5) {
					ptrdiff_t l = std::distance(x.begin(), x_it);
					// divide by 4*pi to get [gamma]/sr = ph/s/sr
					double NO_rad = get_gamma_j(j, i) * 0.25 * M_1_PI;
					double xdiff = *(x_it + 1) - *x_it;
					double x0i = *x_it - 0.5 * xdiff;
					double w = int_slit_func(0.22, NO_wl - x0i)
						- int_slit_func(0.22, NO_wl - x0i - xdiff);
					spec_scia_res.at(l) += w * NO_rad; // * xdiff / npix;
				}
			}
		}
	}
	// set the last element to zero, just to be sure.
	spec_scia_res.back() = 0.;
	spec_max = std::max_element(spec_scia_res.begin(), spec_scia_res.end());
	ptrdiff_t midx = std::distance(spec_scia_res.begin(), spec_max);
	scia_wl_at_max = x.at(midx);
	/* integrated band emission,
	 * corrects for the 1/(4.*pi) for the whole solid angle,
	 * and the d(lambda) = 0.11 might not be fully accurate
	 * but is close enough. */
	scia_band_emiss = 4. * M_PI * 0.11 *
		std::accumulate(spec_scia_res.begin(), spec_scia_res.end(), 0.);
}

void NO_emiss::print_lines_emiss_absorp()
{
	int i, j;

	std::cout << "Wavelength of absorption transition, by lower-state K"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		std::cout << quant_K_vec.at(i);
		for (j = 0; j <= 11; j++) {
			std::cout << "\t" << lambda_K_abs(j, i);
		}
		std::cout << std::endl;
	}

	std::cout << "Wavelength of emission transition, by lower-state K"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		std::cout << quant_K_vec.at(i);
		for (j = 0; j <= 11; j++) {
			std::cout << "\t" << lambda_K(j, i);
		}
		std::cout << std::endl;
	}
}

void NO_emiss::print_Hoenl_London_abs()
{
	int i, j;
	double sum;

	std::cout << "Absorption Hoenl-London factor, by ground-state K"
			  << std::endl;
	std::cout << " P1, Q1, R1, P2, Q2, R2, qP21, pQ12, sR21, oP12, rQ21, qR12"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		sum = 0.;
		std::cout << i << "\t" << quant_K_vec.at(i);
		for (j = 0; j <= 11; j++) {
			sum += vf_HL_K(j, i);
			std::cout << "\t" << vf_HL_K(j, i);
		}
		std::cout << "\t" << sum << std::endl;
	}
}
void NO_emiss::print_Hoenl_London_emiss()
{
	int i, j;
	double sum;

	std::cout << "Emission Hoenl-London factor, by upper-state K"
			  << std::endl;
	std::cout << " P1, Q1, R1, P2, Q2, R2, qP21, pQ12, sR21, oP12, rQ21, qR12"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		sum = 0.;
		std::cout << i << "\t" << quant_K_vec.at(i);
		for (j = 0; j <= 11; j++) {
			sum += vf_HL_emiss_K(j, i);
			std::cout << "\t" << vf_HL_emiss_K(j, i);
		}
		std::cout << "\t" << sum << std::endl;
	}
}
void NO_emiss::print_solar_data()
{
	int i, j;

	std::cout << "K_l, Solar radiation at absorption lines, "
			  << "NJ/N: K_l = J_l - 0.5, K_l = J_l + 0.5"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		std::cout << quant_K_vec.at(i);
		for (j = 0; j <= 11; j++) {
			std::cout << "\t" << solar(j, i);
		}
		std::cout << "\t" << NJ_to_N(0, i);
		std::cout << "\t" << NJ_to_N(1, i);
		std::cout << std::endl;
	}
}

void NO_emiss::print_excitation()
{
	int i;
	double j_u;

	std::cout << "Excitation of upper state J-levels, K_u = J_u - 0.5, K_u = J_u + 0.5"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		j_u = i + 0.5;
		std::cout << j_u;
		std::cout << "\t" << excit(0, j_u - 0.5);
		if (i < NJ) std::cout << "\t" << excit(1, j_u + 0.5);
		std::cout << std::endl;
	}
}

void NO_emiss::print_line_emissivities()
{
	int i, j;

	std::cout << "Gamma-factor at " << Temp << " K, by lower state K"
			  << std::endl;

	for (i = 0; i <= NJ; i++) {
		std::cout << i;
		for (j = 0; j < 12; j++) {
			std::cout << "\t" << gamma_j(j, i);
		}
		std::cout << std::endl;
	}
	std::cout << "band emission rate factor of the " << v_u << "-" << v_l
		<< " transition at " << Temp << " K, photons molec-1 s-1:" << std::endl;
	std::cout << emiss_tot << std::endl;
}

int NO_emiss::get_NJ()
{
	return NJ;
}
int NO_emiss::get_vu()
{
	return v_u;
}
int NO_emiss::get_vl()
{
	return v_l;
}
int NO_emiss::get_vl_abs()
{
	return v_l_abs;
}
double NO_emiss::get_lambda_K(int i, int j)
{
	// translate from Å to nm
	return 0.1 * lambda_K(i, j);
}
double NO_emiss::get_gamma_j(int i, int j)
{
	// convert from 1/Å to 1/nm
	return 10. * gamma_j(i, j);
}
double NO_emiss::get_spec_scia_res(int i)
{
	return spec_scia_res.at(i);
}
double NO_emiss::get_spec_scia_max()
{
	return *std::max_element(spec_scia_res.begin(), spec_scia_res.end());
}
double NO_emiss::get_scia_wl_at_max()
{
	return scia_wl_at_max;
}
double NO_emiss::get_wl_abs_median()
{
	// convert the lambdas (MPL_Matrix elements) to a std::vector
	std::vector<double> lambda_abs(lambda_K_abs.m_Elemente,
			lambda_K_abs.m_Elemente + (NJ + 1) * 12);
	// remove all entries equal to zero
	lambda_abs.erase(std::remove_if(lambda_abs.begin(), lambda_abs.end(),
				std::bind2nd(std::equal_to<double>(), 0.)),
			lambda_abs.end());
	// sort the vector from the shortest to the longest absorption wavelength
	std::sort(lambda_abs.begin(), lambda_abs.end());

	// return the median in nm
	return 0.1 * lambda_abs.at(lambda_abs.size() / 2);
}
double NO_emiss::get_wl_abs_vu_0()
{
	return 0.1 * f_lam(v_u, 0);
}
double NO_emiss::get_wl_emiss_vu_vl()
{
	return 0.1 * f_lam(v_u, v_l);
}
double NO_emiss::get_band_emiss()
{
	return emiss_tot;
}
double NO_emiss::get_scia_band_emiss()
{
	return scia_band_emiss;
}
