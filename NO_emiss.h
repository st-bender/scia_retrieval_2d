/*
 * NO_emiss.h
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

#ifndef NO_EMISS_H_
#define NO_EMISS_H_

#include <vector>
#include "MPL_Matrix.h"

namespace NO_const {
	// all energy values in cm^{-1}
	const double E_u = 43965.7;
	// Vibrational parameters, from Hertzberg, 1950, and Huber and Hertzberg, 1976
	// same as in Eparvier and Barth
	const double w_l = 1904.12177;
	const double w_Xl = 14.09231;
	const double w_Yl = 0.011307;
	const double w_Zl = 3.315e-4;
	const double w_u = 2374.31;
	const double w_Xu = 16.106;
	const double w_Yu = -0.0465;
	// for data from Luque et al.
	const int l_vu = 3;
	const int l_vl = 10;
}

class NO_emiss {
	private:
	int v_u, v_l, v_l_abs; // vibrational quantum number of upper and lower state
	int NJ; // maximal J value
	double Temp; // temperature
	double E_u; // upper state electronical energy, from Eparvier and Barth
	// rotational parameters of the upper state, from Eparvier and Barth
	double B_Vu, D_Vu, gam_u;
	// rotational parameters of the lower state, from Eparvier and Barth
	double B_Vl, D_Vl, A_Vl, Y_l;
	// lower state of the absorption
	double B_Vl_abs, D_Vl_abs, A_Vl_abs, Y_l_abs;
	double lambda_l;
	// lower state mean energy, no fine-splitting
	double E_l, E_l_abs;
	// total electric energy
	double E_tot, E_tot_abs;
	// vibrational term energy
	double E_vib_l, E_vib_l_abs, E_vib_u;
	double E_vib, W_vib, E_vib_abs, W_vib_abs;
	// more variables...
	double f_boltz; // hc/kT
	// electronic partitioning function considering only the unsplit ground-state
	double part_el;
	double phi; // degeneracy of rotational states
	double sum_j; // total rotational-vibration-electric partitioning sum
	// vectors...
	MPL_Matrix F_l, F_l_abs, F_u;
	// partitioning of J-excited lower electronic state
	MPL_Matrix NJ_to_N;
	// line positions of emission and absorption transition
	MPL_Matrix xlines_J, xlines_K, lambda_K;
	MPL_Matrix xlines_K_abs, lambda_K_abs;
	// upper state quanta N, the upper state being strictly Hund's case b
	MPL_Matrix quant_N;
	MPL_Matrix quant_J; // lower state quante J=Lamda + Omega
	// lower state quante K=J+S, K=J for lowest quantum numbers because
	// ground state is intermediate between Hund's case a and b - meaning
	// that for low quantum numbers, S is allayed around the internuclear
	// axis (case a), for high quantum numbers, it is not (case b).
	MPL_Matrix quant_K;
	std::vector<double> quant_K_vec;
	// space for the Hoenl-London factors
	MPL_Matrix vf_HL_J, vf_HL_K;
	MPL_Matrix vf_HL_emiss, vf_HL_emiss_K, quant_j_up;
	// coefficient matrices for data from Luque et al.
	double f_FC_v;
	MPL_Matrix f_osc; // band oscillator strength
	MPL_Matrix f_FC;  // Franck-Condon factor
	MPL_Matrix f_A;   // Einsteins spontaneous emission probability
	MPL_Matrix f_lam; // band wavelength [nm]
	// the excitation of J'
	MPL_Matrix excit;
	MPL_Matrix excit_pqr;
	// line emissivities by upper state k
	double emiss_tot;
	MPL_Matrix gamma_j;
	double scia_wl_at_max;
	double scia_band_emiss;

	public:
	// solar spectrum
	MPL_Matrix solar;
	std::vector<double> spec_scia_res;
	explicit NO_emiss(int vu = 2, int vl = 4, int vl_abs = 0, double T = 200.);
	void alloc_memory();
	void set_constants();
	void populate_Fs();
	void calc_lines_emiss_absorp();
	void set_Hoenl_London_abs();
	void set_Hoenl_London_emiss();
	void get_solar_data(class Sonnenspektrum &sol_spec);
	int read_luque_data_from_file(std::string filename);
	void calc_excitation();
	void calc_line_emissivities();
	std::vector<double> calc_polarisation();
	void pol_corr(double SZA, double rel_SAA, double mu2 = 0.17, double mu3 = -0.2);
	void scia_convolve(class Messung &ml);
	void print_lines_emiss_absorp();
	void print_Hoenl_London_abs();
	void print_Hoenl_London_emiss();
	void print_excitation();
	void print_line_emissivities();
	void print_solar_data();
	int get_NJ();
	int get_vu();
	int get_vl();
	int get_vl_abs();
	double get_lambda_K(int i, int j);
	double get_gamma_j(int i, int j);
	double get_spec_scia_res(int i);
	double get_spec_scia_max();
	double get_scia_wl_at_max();
	double get_wl_abs_median();
	double get_wl_abs_vu_0();
	double get_wl_emiss_vu_vl();
	double get_band_emiss();
	double get_scia_band_emiss();
};

#endif /* NO_EMISS_H_ */
