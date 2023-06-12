//
// Created by finge on 6/5/2023.
//

#ifndef PULSERECONSTRUCTION_SUB_MATRIX_ELEMENTS_H
#define PULSERECONSTRUCTION_SUB_MATRIX_ELEMENTS_H

#include "physics_functions.h"
#include <iostream>

complex<double> TwoPhotonAmplitude(double initial_energy, int N_states, vector<double> middle_dipole, vector<double> middle_energies,
                                   vector<double> final_dipole, double final_energy, vector<double> laser_one_parameters,
                                   vector<double> laser_two_parameters, double position_one, double position_two);

complex<double> OnePhotonAmplitude(double initial_energy, double dipole, double final_energy, vector<double> laser_parameters,
                                   double position);

#endif //PULSERECONSTRUCTION_SUB_MATRIX_ELEMENTS_H
