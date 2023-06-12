//
// Created by finge on 6/5/2023.
//

#ifndef PULSERECONSTRUCTION_AUTOCORRELATION_H
#define PULSERECONSTRUCTION_AUTOCORRELATION_H

#include "physics_functions.h"
#include "sub_matrix_elements.h"
#include <cstring>
#include <iostream>
#include <fstream>

vector<double> Autocorrelation(double initial_energy, int N_free_states_l0, vector <vector<double>> two_photon_dipoles_l0,
                               vector<double> free_energies_l0, int N_free_states_l2,
                               vector <vector<double>> two_photon_dipoles_l2, vector<double> free_energies_l2,
                               int N_free_states_l1, vector<double> one_photon_dipoles_l1, vector<double> free_energies_l1,
                               int N_bound_states_l1, const vector<double>& two_photon_dipoles_l1, const vector<double>& bound_energies_l1,
                               int N_gaussians, vector <vector<double>> laser_list, int N_delays,
                               vector<double> correlation_delay);

#endif //PULSERECONSTRUCTION_AUTOCORRELATION_H
