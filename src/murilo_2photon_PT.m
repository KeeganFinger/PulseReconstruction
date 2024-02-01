%===== Clear Previous Data ==================
close all; clear all;
%%
%===== Read Data ============================
N_states = 500;
data_dir = './data/';
for l=0:2
    %===== Read Energies ====================
    data = h5read([data_dir 'He.h5'],char("/Energy_l"+l));
    if l==0
        nmax = length(data(:,1));
        energies = zeros(nmax,3);
    end
    energies(1:nmax-l,l+1) = data(:,1) + 1i*data(:,2);
    %===== Read States ======================
    data = h5read([data_dir 'He.h5'],char("/Psi_r_l"+l));
    if l==0
        rmax = length(data(:,1,1));
        states = zeros(rmax,nmax,3);
    end
    states(:,1:nmax-l,l+1) = data(:,:,1) + 1i*data(:,:,2);
end % Loop of final L
%========== Read Grid Parameters ============
grid = h5read([data_dir 'parameters.h5'],char("/EPS/r"));
dr = h5read([data_dir 'parameters.h5'],char("/EPS/delta_x"));
%%
%===== Organize Atomic States ===============
%========== L=0 Block =======================
angular_momentum = 0;
l0_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l0_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l0 = length(l0_free_energies);
%========== L=1 Block =======================
angular_momentum = 1;
l1_bound_energies = real(energies(real(energies(:,angular_momentum+1)) < 0,angular_momentum+1));
l1_bound_states = states(:,real(energies(:,angular_momentum+1)) < 0,angular_momentum+1);
N_bound_states_l1 = length(l1_bound_energies);
l1_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l1_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l1 = length(l1_free_energies);
%========== L=2 Block =======================
angular_momentum = 2;
l2_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l2_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l2 = length(l2_free_energies);
%========== Initial State ===================
initial_energy = real(energies(1,1));
initial_state = states(:,1,1);

%%
%===== Precompute Dipoles ===================
one_photon_dipoles_l1 = zeros(length(l1_free_energies),1);
for state = 1:length(l1_free_energies)
    one_photon_dipoles_l1(state) = sum(l1_free_states(:,state).*grid.*initial_state).*sqrt(3).*Wigner3j(1,1,0,0,0,0).^2;
end % Loop over free L=1 states
two_photon_dipoles_l1 = zeros(length(l1_bound_energies),1);
two_photon_dipoles_l0 = zeros(length(l1_bound_energies),length(l0_free_energies));
two_photon_dipoles_l2 = zeros(length(l1_bound_energies),length(l2_free_energies));
for state_one = 1:length(l1_bound_energies)
    for state_two = 1:length(l0_free_energies)
        two_photon_dipoles_l0(state_one,state_two) = sum(l1_bound_states(:,state_one).*grid.*l0_free_states(:,state_two)).*sqrt(3).*Wigner3j(0,1,1,0,0,0).^2;
    end % Loop over free L=0 states
    for state_two = 1:length(l2_free_energies)
        two_photon_dipoles_l2(state_one,state_two) = sum(l1_bound_states(:,state_one).*grid.*l2_free_states(:,state_two)).*sqrt(15).*Wigner3j(2,1,1,0,0,0).^2;
    end % Loop over free L=2 states
    two_photon_dipoles_l1(state_one) = sum(l1_bound_states(:,state_one).*grid.*initial_state).*sqrt(3).*Wigner3j(1,1,0,0,0,0).^2;
end % Loop over bound L=1 states

%%
%===== Pulse Setup ==========================
laser_amplitude_SI  = 5e11;                                                 % [W / cm^2]
laser_amplitude_au  = Laser.SI2au_intensity(laser_amplitude_SI);            % [a.u.]

laser_wavelength    = 800;                                                  % [nm]
photon_energy       = Laser.SI2au_wavelength(laser_wavelength);             % [eV]

laser_FWHM_SI       = 1e4;                                                  % Full-Width at Half-Max of Gaussian [fs]
laser_FWHM_au       = Laser.SI2au_duration(laser_FWHM_SI);                  % Full-Width at Half-Max of Gaussian [a.u.]

start_wavelength    = 800;
end_wavelength      = 800;
laser_chirp         = (Laser.SI2au_wavelength(start_wavelength)...
                        - Laser.SI2au_wavelength(end_wavelength))... 
                        / (2 * sqrt(2 * log(2)) *laser_FWHM_au);            % Change in frequency [1]

laser_position      = 0;                                                    % Location of center of pulse [a.u.]

laser_phase         = 0;                                                    % Phase offset of laser [rad.]

% This is the collection of Gaussian's making your pulse
laser_pulse = [
    Laser(laser_amplitude_au, photon_energy, laser_FWHM_au, laser_chirp, laser_position, laser_phase);
];

%%
%===== Autocorrelation Calculation ==========
calc = @(basis,delay) squeeze(sum(sum(abs( ...
        matrixElementsCalculation(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,delay,[],true,size(basis,1),[])).^2,2),3));

% TODO:
% ~ Plot the input laser pulse and see what it looks like
%     ~ What do typical lasers for this sort of experiment look like? Does yours agree?
%     ~ Hint: look at the Laser class. What function(s) might be useful?
% ~ Calculate an autocorrelation spectrum using your pulse
%     ~ Try using different input lasers. How does the spectrum depend on the laser parameters?
%     ~ Hint: the delay parameter is the x-axis on the plots we've talked about
%     ~ Hint: look at the Laser class. What function can you use to extract the parameters
%     ~ Hint: is your spectrum all zero? What laser parameter can you change to ionize more electrons?