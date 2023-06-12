//
// Created by finge on 6/5/2023.
//

#include "autocorrelation.h"

vector<double> Autocorrelation(double initial_energy, int N_free_states_l0, vector <vector<double>> two_photon_dipoles_l0,
                       vector<double> free_energies_l0, int N_free_states_l2,
                       vector <vector<double>> two_photon_dipoles_l2, vector<double> free_energies_l2,
                       int N_free_states_l1, vector<double> one_photon_dipoles_l1, vector<double> free_energies_l1,
                       int N_bound_states_l1, const vector<double>& two_photon_dipoles_l1, const vector<double>& bound_energies_l1,
                       int N_gaussians, vector <vector<double>> laser_list, int N_delays,
                       vector<double> correlation_delay) {
    vector<double> signal;
    signal.reserve(correlation_delay.size());
    for (int delay = 0; delay < N_delays; delay++) {
        double this_delay = correlation_delay.at(delay);

        // L=0 Final State
        double autocorrelation_signal = 0;
        for (int end_state = 0; end_state < N_free_states_l0; end_state++) {
            double final_energy = free_energies_l0.at(end_state);
            for (int n = 0; n < N_gaussians; n++) {
                vector<double> laser_one = laser_list.at(n);
                double omega_one = laser_one.at(0);
                for (int m = 0; m < N_gaussians; m++) {
                    vector<double> laser_two = laser_list.at(m);
                    double omega_two = laser_two.at(0);
                    // Final state L=0
                    if (final_energy < 3 * max(omega_one, omega_two) + initial_energy) {
                        double E_fi = final_energy - initial_energy;

                        complex<double> term_one =
                                TwoPhotonAmplitude(initial_energy, N_bound_states_l1, two_photon_dipoles_l1,
                                                   bound_energies_l1, two_photon_dipoles_l0.at(end_state), final_energy,
                                                   laser_one, laser_two, laser_one.at(1), laser_two.at(1)) *
                                exp(1i * E_fi * this_delay);
                        complex<double> term_two =
                                TwoPhotonAmplitude(initial_energy, N_bound_states_l1, two_photon_dipoles_l1,
                                                   bound_energies_l1, two_photon_dipoles_l0.at(end_state), final_energy,
                                                   laser_one, laser_two, laser_one.at(1) + this_delay,
                                                   laser_two.at(1)) * exp(1i * omega_one * this_delay);
                        complex<double> term_three =
                                TwoPhotonAmplitude(initial_energy, N_bound_states_l1, two_photon_dipoles_l1,
                                                   bound_energies_l1, two_photon_dipoles_l0.at(end_state), final_energy,
                                                   laser_one, laser_two, laser_one.at(1),
                                                   laser_two.at(1) + this_delay) * exp(1i * omega_two * this_delay);
                        autocorrelation_signal += norm(term_one + term_two + term_three);
                    }
                }
            }
        }

        // L=2 Final State
        for (int end_state = 0; end_state < N_free_states_l2; end_state++) {
            double final_energy = free_energies_l2.at(end_state);
            for (int n = 0; n < N_gaussians; n++) {
                vector<double> laser_one = laser_list.at(n);
                double omega_one = laser_one.at(0);
                for (int m = 0; m < N_gaussians; m++) {
                    vector<double> laser_two = laser_list.at(m);
                    double omega_two = laser_two.at(0);
                    // Final state L=0
                    if (final_energy < 3 * max(omega_one, omega_two) + initial_energy) {
                        double E_fi = final_energy - initial_energy;

                        complex<double> term_one =
                                TwoPhotonAmplitude(initial_energy, N_bound_states_l1, two_photon_dipoles_l1,
                                                   bound_energies_l1, two_photon_dipoles_l2.at(end_state), final_energy,
                                                   laser_one, laser_two, laser_one.at(1), laser_two.at(1)) *
                                exp(1i * E_fi * this_delay);
                        complex<double> term_two =
                                TwoPhotonAmplitude(initial_energy, N_bound_states_l1, two_photon_dipoles_l1,
                                                   bound_energies_l1, two_photon_dipoles_l2.at(end_state), final_energy,
                                                   laser_one, laser_two, laser_one.at(1) + this_delay,
                                                   laser_two.at(1)) * exp(1i * omega_one * this_delay);
                        complex<double> term_three =
                                TwoPhotonAmplitude(initial_energy, N_bound_states_l1, two_photon_dipoles_l1,
                                                   bound_energies_l1, two_photon_dipoles_l2.at(end_state), final_energy,
                                                   laser_one, laser_two, laser_one.at(1),
                                                   laser_two.at(1) + this_delay) * exp(1i * omega_two * this_delay);
                        autocorrelation_signal += norm(term_one + term_two + term_three);
                    }
                }
            }
        }

        // L=1 Final State
        for (int end_state = 0; end_state < N_free_states_l1; end_state++) {
            double final_energy = free_energies_l1.at(end_state);
            double E_fi = final_energy - initial_energy;
            for (int n = 0; n < N_gaussians; n++) {
                vector<double> laser = laser_list.at(n);
                autocorrelation_signal += norm(
                        OnePhotonAmplitude(initial_energy, one_photon_dipoles_l1.at(end_state), final_energy, laser,
                                           laser.at(1)) * (1.0 + exp(1i * E_fi * this_delay)));
            }
        }
        signal.emplace_back(autocorrelation_signal);
    }
    signal.shrink_to_fit();
    return signal;
}