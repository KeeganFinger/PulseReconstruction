#include <algorithm>
#include "autocorrelation.h"
#include <H5Cpp.h>
#include <string>
using namespace std;
using namespace H5;

double SI2au_wavelength(double wavelength) {
    return 1239.8 / wavelength / 27.211;
}

int main() {

    // Runtime Parameters - to be put in input file

    double max_intensity = 3e-3;
    vector<int> known_harmonics{}; vector<int> unknown_harmonics{};
    double correlation_delay_min = -1000; double correlation_delay_max = 1000;
    int N_delays = 20001;
    double dt = (correlation_delay_max - correlation_delay_min) / N_delays;
    vector<double> correlation_delay; correlation_delay.reserve(N_delays);
    for (int t = 0; t < N_delays; t++) {
        correlation_delay.emplace_back(t * dt);
    }

    vector<int> reconstruction_gaussian_list{1};
    for (auto &g : reconstruction_gaussian_list) {
        g *= unknown_harmonics.size();
    }
    bool chirp = false; bool fit_color = false; bool cross_correlation = false;

    // Read data from h5 files
//    const H5std_string FILENAME ("CppData.h5");
//    const H5std_string LASER_GROUP ("/LASER");
//    const H5std_string ATOM_GROUP ("/ATOM");
//
//    H5File *datafile = new H5File(FILENAME, H5F_ACC_RDONLY);
//    Group *laser = new Group(datafile->openGroup(LASER_GROUP));
//    Attribute *harm7 = new Attribute(laser->openAttribute("HARM7"));
//    DataType *harm7_type = new DataType(harm7->getDataType());
//
//    vector<vector<double>> gaussian_train_7;
//    harm7->read(*harm7_type,&gaussian_train_7);
//    for (auto i : gaussian_train_7) {
//        for (auto j : i) {
//            cout << j << ' ';
//        }
//        cout << '\n';
//    }
//    return 0;

    int N_free_states_l0; vector<vector<double>> two_photon_dipoles_l0;
    int N_free_states_l2; vector<vector<double>> two_photon_dipoles_l2;
    int N_bound_states_l1; vector<double> two_photon_dipoles_l1;
    int N_free_states_l1; vector<double> one_photon_dipoles_l1;

    double initial_energy;
    vector<double> free_energies_l0;
    vector<double> free_energies_l1;
    vector<double> free_energies_l2;
    vector<double> bound_energies_l1;

    vector<vector<double>> gaussian_train_7;
    vector<vector<double>> gaussian_train_9;
    vector<vector<double>> gaussian_train_11;
    vector<vector<double>> gaussian_train_13;

    vector<vector<double>> experiment;
    for (auto gaussian : unknown_harmonics) {
        if (gaussian == 7) {
            for (auto g : gaussian_train_7) {
                experiment.push_back(g);
            }
        }
        if (gaussian == 9) {
            for (auto g : gaussian_train_9) {
                experiment.push_back(g);
            }
        }
        if (gaussian == 11) {
            for (auto g : gaussian_train_11) {
                experiment.push_back(g);
            }
        }
        if (gaussian == 13) {
            for (auto g : gaussian_train_13) {
                experiment.push_back(g);
            }
        }
    }
    if (experiment.empty()) return 1;
    double rescale_factor = 1;
    for (auto gaussian : experiment) {
        gaussian[2] *= rescale_factor;
        gaussian[3] *= rescale_factor;
    }

    // Reconstruction parameters
    vector<double> omega; omega.reserve(unknown_harmonics.size());
    for (auto color : unknown_harmonics) {
        omega.emplace_back(double(color) * SI2au_wavelength(800));
    }
    int N_windows = 10; int percentages[] = {2, 5, 10, 15, 25, 50, 75, 100};
    double window_width = *max_element(correlation_delay.begin(), correlation_delay.end());
    double sample_step = 0.5 * 2 * M_PI/ (*max_element(omega.begin(),omega.end()));

    vector<double> known_signal = Autocorrelation(initial_energy,
                                    N_free_states_l0,two_photon_dipoles_l0,free_energies_l0,
                                    N_free_states_l2,two_photon_dipoles_l2,free_energies_l2,
                                    N_free_states_l1,one_photon_dipoles_l1,free_energies_l1,
                                    N_bound_states_l1,two_photon_dipoles_l1,bound_energies_l1,
                                    experiment.size(),experiment,N_delays,correlation_delay);
    ofstream outFile;
    outFile.open("autocorrelation.dat");
    for (int i = 0; i < correlation_delay.size(); i++){
        outFile << correlation_delay[i] << "\t" << known_signal[i] << "\n";
    }
    outFile.close();
}
