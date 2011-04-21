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
NO_emiss::NO_emiss()
{
	v_u = 2;
	v_l = 4;
	v_l_abs = 0;
	Temp = 200.;
	NJ = 37;

	alloc_memory();
	set_constants();
}

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

	// copy
	F_l = Fl;
	F_l_abs = Fl_abs;
	F_u = Fu;
	NJ_to_N = NjtoN;

	xlines_K = xl_K;
	xlines_K_abs = xl_K_abs;
	lambda_K = lam_K;
	lambda_K_abs = lam_K_abs;

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
	f_boltz = 1.0 / (phys::k2 * Temp);
	part_el = 1.0 + std::exp(-f_boltz * NO_const::E_u);

	// Rotational and total partitioning sum
	sum_j = 0;
	for (int i = 0; i <= NJ + 2; i++) {
		double j = i + 0.5;
		sum_j += (2. * j + 1.) * std::exp(-f_boltz * B_Vl_abs * j * (j + 1.));
	}
	// total partitioning sum considering vibrational ground state
	// and Lambda splitting (Hund's case a)
	sum_j *= part_el * (std::exp(-f_boltz * E_l_abs)
			+ std::exp(+f_boltz * E_l_abs));

	return 0;
}
int NO_emiss::populate_Fs()
{
	for (int i = 0; i < NJ + 2; i++) {
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
				NJ_to_N(1, i) = phi * (2. * j + 1.) / sum_j
					* std::exp(-f_boltz * (flow + E_l_abs));
		}
	}

	return 0;
}

int NO_emiss::calc_lines_emiss_absorp()
{
	int i, l;
	double E_rot, E_rot_abs;

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

			if (k_u >=0 && j_u >= 0.5 && j_l >= 0.5) {
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

	return 0;
}
