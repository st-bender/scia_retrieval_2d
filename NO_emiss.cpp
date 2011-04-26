/*
 * NO_emiss.cpp
 *
 *  Created on: 20-Apr-2011
 *      Author: bender-s
 */
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "constants.h"
#include "NO_emiss.h"

int line_count(std::string filename)
{
	std::ifstream file(filename.c_str());
	return std::count(std::istreambuf_iterator<char>(file),
			std::istreambuf_iterator<char>(), '\n') + 1;
}
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
NO_emiss::NO_emiss(int vu, int vl, int vl_abs, double T)
{
	v_u = vu;
	v_l = vl;
	v_l_abs = vl_abs;
	Temp = T;
	NJ = 37;

	alloc_memory();
	set_constants();
}

int NO_emiss::alloc_memory()
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

	vf_HL_emiss = vfHL_emiss;
	vf_HL_emiss_K = vfHL_emiss_K;
	quant_j_up = quantjup;

	gamma_j = gam_j;

	return 0;
}

int NO_emiss::set_constants()
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

	E_vib_u = NO_const::w_u * (v_u + 0.5)
		- NO_const::w_Xu * (v_u + 0.5) * (v_u + 0.5)
		+ NO_const::w_Yu * (v_u + 0.5) * (v_u + 0.5) * (v_u + 0.5);
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

	return 0;
}
int NO_emiss::populate_Fs()
{
	NJ_to_N.Null_Initialisierung();
	F_l.Null_Initialisierung();
	F_l_abs.Null_Initialisierung();
	F_u.Null_Initialisierung();

	for (int i = 0; i <= NJ + 2; i++) {
		double j = i + 0.5;
		double flow = B_Vl_abs * j * (j + 0.5);
		F_l(0, i) = F0(B_Vl, D_Vl, Y_l, lambda_l, j);
		F_l(1, i) = F1(B_Vl, D_Vl, Y_l, lambda_l, j);
		F_l_abs(0, i) = F0(B_Vl_abs, D_Vl_abs, Y_l_abs, lambda_l, j);
		F_l_abs(1, i) = F1(B_Vl_abs, D_Vl_abs, Y_l_abs, lambda_l, j);
		
		F_u(0, i) = B_Vu * i * (i + 1)
			- D_Vu * i * i * (i + 1.) * (i + 1)
			+ 0.5 * gam_u * i;
		F_u(1, i) = B_Vu * i * (i + 1)
			- D_Vu * i * i * (i + 1.) * (i + 1)
			- 0.5 * gam_u * (i + 1.);

		if (i <= NJ) {
			NJ_to_N(0, i) = phi * (2. * j + 1.) / sum_j
				* std::exp(-f_boltz * (flow - E_l_abs));
			if (i > 0 && i < NJ)
				NJ_to_N(1, i + 1) = phi * (2. * j + 1.) / sum_j
					* std::exp(-f_boltz * (flow + E_l_abs));
		}
	}

	return 0;
}

int NO_emiss::calc_lines_emiss_absorp()
{
	int i, j, l;
	double E_rot, E_rot_abs;

	xlines_K.Null_Initialisierung();
	xlines_K_abs.Null_Initialisierung();
	lambda_K.Null_Initialisierung();
	lambda_K_abs.Null_Initialisierung();

	for (i = 0; i <= NJ; i++) {
		double k_l = i;
		double j_l, j_u, k_u;
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

	return 0;
}

// Absorption Hoenl-London factors from ground state
// from Earls, 1935, normalised to 4*(2J+1)
int NO_emiss::set_Hoenl_London_abs()
{
	int i;
	double j, u, d1, d2, d3;

	vf_HL_J.Null_Initialisierung();
	vf_HL_K.Null_Initialisierung();

	// j = 0.5
	// P1
	vf_HL_J(0, 0) = vf_HL_K(0, 0) = 0.;
	// Q1
	vf_HL_J(1, 0) = vf_HL_K(1, 0) = 4. / 3.;
	// R1
	vf_HL_J(2, 0) = vf_HL_K(2, 0) = 4. / 6.;
	// P2
	vf_HL_J(3, 0) = vf_HL_K(3, 1) = 0.;
	// Q2
	vf_HL_J(4, 0) = vf_HL_K(4, 1) = 0.;
	// R2
	vf_HL_J(5, 0) = vf_HL_K(5, 1) = 0.;
	// qP21
	vf_HL_J(6, 0) = vf_HL_K(6, 0) = 0.;
	// pQ12
	vf_HL_J(7, 0) = vf_HL_K(7, 1) = 0.;
	// sR21
	vf_HL_J(8, 0) = vf_HL_K(8, 0) = 4. / 6.;
	// oP12
	vf_HL_J(9, 0) = vf_HL_K(9, 1) = 0.;
	// rQ21
	vf_HL_J(10, 0) = vf_HL_K(10, 0) = 4. / 3.;
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
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		vf_HL_K(0, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		// Q1
		vf_HL_J(1, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				+ u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		vf_HL_K(1, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				+ u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		// R1
		vf_HL_J(2, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		vf_HL_K(2, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		// P2
		vf_HL_J(3, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
		if (i != NJ)
			vf_HL_K(3, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					+ (2.*j + 1.) * u * (4.*j*j + 4.*j - 7 + 2.*Y_l)) / d1;
		// Q2
		vf_HL_J(4, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				+ u * (8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l)) / d3;
		if (i != NJ)
			vf_HL_K(4, i + 1) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
					+ u * (8.*j*j*j + 12.*j*j - 2.*j + 1. - 2.*Y_l)) / d3;
		// R2
		vf_HL_J(5, i) = ((2.*j + 1.) * (2.*j + 1.)
				+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1 - 2.*Y_l)) / d2;
		if (i != NJ)
			vf_HL_K(5, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					+ (2.*j + 1.) * u * (4.*j*j + 4.*j + 1 - 2.*Y_l)) / d2;
		// qP21
		vf_HL_J(6, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
		vf_HL_K(6, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d1;
		// pQ12
		vf_HL_J(7, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				- u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2*Y_l)) / d3;
		if (i != NJ)
			vf_HL_K(7, i + 1) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
					- u * (8.*j*j*j + 12.*j*j - 2.*j - 7. + 2.*Y_l)) / d3;
		// sR21
		vf_HL_J(8, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d2;
		vf_HL_K(8, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d2;
		// oP12
		vf_HL_J(9, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		if (i != NJ)
			vf_HL_K(9, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					- (2.*j + 1.) * u * (4.*j*j + 4.*j + 1. - 2.*Y_l)) / d1;
		// rQ21
		vf_HL_J(10, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				- u * (8.*j*j*j + 12.*j*j - 2.*j + 1 - 2.*Y_l)) / d3;
		vf_HL_K(10, i) = (2.*j + 1.) * ((4.*j*j + 4.*j - 1.)
				- u * (8.*j*j*j + 12.*j*j - 2.*j + 1 - 2.*Y_l)) / d3;
		// qR12
		vf_HL_J(11, i) = ((2.*j + 1.) * (2.*j + 1.)
				- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
		if (i != NJ)
			vf_HL_K(11, i + 1) = ((2.*j + 1.) * (2.*j + 1.)
					- (2.*j + 1.) * u * (4.*j*j + 4.*j - 7. + 2.*Y_l)) / d2;
	}

	return 0;
}
// Hoenl-London factors of emission, relative to upper state K
// and upper state J ? anyways, these are normalised to 2*(2J+1)
int NO_emiss::set_Hoenl_London_emiss()
{
	int i;
	double j, u, d1, d2, d3;

	vf_HL_emiss.Null_Initialisierung();
	vf_HL_emiss_K.Null_Initialisierung();

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

	return 0;
}

double interpolate(std::vector<double> &x, std::vector<double> &y, double x0)
{
	int i;
	std::vector<double>::iterator x_it;
	x_it = std::upper_bound(x.begin(), x.end(), x0);

	if (x_it == x.begin()) return y.at(0);
	if (x_it == x.end()) return *(y.end() - 1);

	i = distance(x.begin(), x_it) - 1;

	return y.at(i)
		+ (x0 - x.at(i)) * (y.at(i) - y.at(i + 1)) / (x.at(i) - x.at(i + 1));
}

// Get solar spectral data
int NO_emiss::get_solar_data(Sonnenspektrum &sol_spec)
{
	int i, j;

	solar.Null_Initialisierung();

	for (i = 0; i <= NJ; i++)
		for (j = 0; j <= 11; j++)
			if (lambda_K_abs(j, i) != 0.) {
				// translate from Å to nm
				solar(j, i) = interpolate(sol_spec.m_Wellenlaengen,
						sol_spec.m_Intensitaeten, 0.1 * lambda_K_abs(j, i));
			}

	return 0;
}

// Get data from Luque et al.
int NO_emiss::read_luque_data_from_file(std::string filename)
{
	int i, j;
	int vu, vl;
	double d_dum;
	std::string s_dum;
	std::ifstream lfile;

	lfile.open(filename.c_str());

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

	lfile.close();

	return 0;
}

// the excitation of J': sum over all J
// Excitation of upper state J-levels: k_u = j_u - 0.5, k_u = j_u + 0.5
int NO_emiss::calc_excitation()
{
	int i, l, k_l, k_u;
	double j, j_l, j_u, nj_frac, sum1, sum2;
	double f_FC_tot;

	excit.Null_Initialisierung();

	for (i = 0; i <= NJ; i++) {
		j = j_u = i + 0.5;
		sum1 = sum2 = 0.0;
		for (l = 0; l <= 11; l++) {
			if (l == 0 || l == 3 || l == 6 || l == 9) j_l = j_u + 1.;
			if (l == 1 || l == 4 || l == 7 || l == 10) j_l = j_u;
			if (l == 2 || l == 5 || l == 8 || l == 11) j_l = j_u - 1.;

			if (l == 0 || l == 1 || l == 2 || l == 6 || l == 8 || l == 10)
				k_l = j_l - 0.5;
			else
				k_l = j_l + 0.5;
			if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
				k_u = j_u - 0.5;
			else
				k_u = j_u + 0.5;

			if ((k_l >= 0) && (k_u >= 0) && (j_l >= 0.5) && (k_l <= NJ)) {
				if (l == 0 || l == 1 || l == 2 || l == 6 || l == 8 || l == 10)
					nj_frac = NJ_to_N(0, j_l - 0.5);
				else
					nj_frac = NJ_to_N(1, j_l + 0.5);

				if ((k_l >= 1) || (k_l == 1 && j_l == 1.5) || (k_l == 0)) {
					if (l == 0 || l == 1 || l == 2 || l == 7 || l == 9 || l == 11)
						sum1 += lambda_K_abs(l, k_l) * lambda_K_abs(l, k_l)
							* solar(l, k_l) * vf_HL_K(l, k_l) / (2. * j_l + 1.)
							* nj_frac;
					else
						sum2 += lambda_K_abs(l, k_l) * lambda_K_abs(l, k_l)
							* solar(l, k_l) * vf_HL_K(l, k_l) / (2. * j_l + 1.)
							* nj_frac;
					}
			}
		}
		sum1 *= phys::flux * f_osc(v_u, v_l);
		sum2 *= phys::flux * f_osc(v_u, v_l);
		excit(0, i) = sum1;
		if (i < NJ) excit(1, i + 1) = sum2;
	}

	f_FC_tot = 0.;
	for (i = 0; i <= NO_const::l_vl; i++) {
		double l_vib = 1.e8 / f_lam(v_u, i);
		f_FC_tot += f_FC(v_u, i) * l_vib*l_vib*l_vib*l_vib*l_vib;
	}
	W_vib = 1.e8 / f_lam(v_u, v_l);
	f_FC_v = f_FC(v_u, v_l) * W_vib*W_vib*W_vib*W_vib*W_vib / f_FC_tot;

	return 0;
}

int NO_emiss::calc_line_emissivities()
{
	int i, l, k_l, k_u;
	double j_l, j_u;

	gamma_j.Null_Initialisierung();

	emiss_tot = 0.;
	for (i = 0; i <= NJ - 2; i++) {
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

	return 0;
}
