//
// Created by finge on 6/5/2023.
//

#include "sub_matrix_elements.h"

complex<double> TwoPhotonAmplitude(double initial_energy, int N_states, vector<double> middle_dipole, vector<double> middle_energies,
                                   vector<double> final_dipole, double final_energy, vector<double> laser_one_parameters,
                                   vector<double> laser_two_parameters, double position_one, double position_two) {
    // Laser One Parameters
    double omega_one = laser_one_parameters.at(0);
    complex<double> period_one =
            1.0 / sqrt(1.0 / pow(laser_one_parameters.at(4),2) + 1i * laser_one_parameters.at(5));
    complex<double> amplitude_one = laser_one_parameters.at(2) + 1i * laser_one_parameters.at(3);

    // Laser Two Parameters
    double omega_two = laser_two_parameters.at(0);
    complex<double> period_two =
            1.0 / sqrt(1.0 / pow(laser_two_parameters.at(4),2) + 1i * laser_two_parameters.at(5));
    complex<double> amplitude_two = laser_two_parameters.at(2) + 1i * laser_two_parameters.at(3);

    // Compute Two-Photon Ionization Amplitude
    complex<double> ionization_amplitude = 0;
    for (int mid_state = 0; mid_state < N_states; mid_state++) {
        double middle_energy = middle_energies.at(mid_state);
        double delta_fm = final_energy - middle_energy - omega_two;
        double delta_mi = middle_energy - initial_energy - omega_one;

        complex<double> prefactor = -0.25 * M_PI * amplitude_two * amplitude_one * period_two * period_one;

        complex<double> laser_factor =
                -0.5 * (period_one * period_one * delta_mi * delta_mi + period_two * period_two * delta_fm * delta_fm) +
                1i * (position_two * delta_fm + position_one * delta_mi);
        complex<double> erf_factor = 1.0 / sqrt(2.0 * (period_one * period_one + period_two * period_two)) *
                                     (position_two - position_one +
                                      1i * (delta_fm * period_two * period_two - delta_mi * period_one * period_one));
        complex<double> laser = exp(1i * laser_factor.imag()) *
                                (exp(laser_factor.real()) + ErrorFunction(erf_factor, laser_factor.real()));

        ionization_amplitude += prefactor * final_dipole.at(mid_state) * middle_dipole.at(mid_state) * laser;
    }
    return ionization_amplitude;
}

complex<double> OnePhotonAmplitude(double initial_energy, double dipole, double final_energy, vector<double> laser_parameters,
                                   double position) {

    // Laser Parameters
    double omega = laser_parameters.at(0);
    complex<double> period = 1.0 / sqrt(1.0 / pow(laser_parameters.at(4),2) + 1i * laser_parameters.at(5));
    complex<double> amplitude = laser_parameters.at(2) + 1i * laser_parameters.at(3);

    // Compute One-Photon Ionization Amplitude
    complex<double> prefactor = -1i * sqrt(M_PI / 2) * amplitude * period;
    double delta_fi = final_energy - initial_energy - omega;
    complex<double> laser_factor = -1i * delta_fi * position - 0.5 * (delta_fi * delta_fi * period * period);

    return prefactor * dipole * exp(laser_factor);
}