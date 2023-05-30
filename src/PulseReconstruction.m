close all; clear all;
%===== Runtime Parameters ===================
gaussian_expand = false; gaussian_basis_size = 7;

max_intensity = 3e-3;
known_harmonics = [11]; unknown_harmonics = [7];
correlation_delay = linspace(-1000,1000,2000);

reconstruction_gaussian_list = [1 2 4 6];
chirp = false; cross_correlation = true;


%%
%===== Read Data ============================
N_states = 500;
for l=0:2
    data = h5read("./data/He.h5",char("/Energy_l"+l));
    if l==0
        nmax = length(data(:,1));
        energies = zeros(nmax,3);
    end
    energies(1:nmax-l,l+1) = data(:,1) + 1i*data(:,2);
    data = h5read("./data/He.h5",char("/Psi_r_l"+l));
    if l==0
        rmax = length(data(:,1,1));
        states = zeros(rmax,nmax,3);
    end
    states(:,1:nmax-l,l+1) = data(:,:,1) + 1i*data(:,:,2);
end

grid = h5read("./data/parameters.h5",char("/EPS/r"));
dr = h5read("./data/parameters.h5",char("/EPS/delta_x"));
%%

%===== L=0 Block ============================
angular_momentum = 0;
l0_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l0_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l0 = length(l0_free_energies);
%===== L=1 Block ============================
angular_momentum = 1;
l1_bound_energies = real(energies(real(energies(:,angular_momentum+1)) < 0,angular_momentum+1));
l1_bound_states = states(:,real(energies(:,angular_momentum+1)) < 0,angular_momentum+1);
N_bound_states_l1 = length(l1_bound_energies);
l1_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l1_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l1 = length(l1_free_energies);
%===== L=2 Block ============================
angular_momentum = 2;
l2_free_energies = real(energies(real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1));
l2_free_states = states(:,real(energies(:,angular_momentum+1)) >= 0,angular_momentum+1);
N_free_states_l2 = length(l2_free_energies);
%===== Initial State ========================
initial_energy = real(energies(1,1));
initial_state = states(:,1,1);

%%
%===== Precompute Dipoles ===================
one_photon_dipoles_l1 = zeros(length(l1_free_energies),1);
for state = 1:length(l1_free_energies)
    one_photon_dipoles_l1(state) = sum(l1_free_states(:,state).*grid.*initial_state).*sqrt(3).*Wigner3j(1,1,0,0,0,0).^2;
end
two_photon_dipoles_l1 = zeros(length(l1_bound_energies),1);
two_photon_dipoles_l0 = zeros(length(l1_bound_energies),length(l0_free_energies));
two_photon_dipoles_l2 = zeros(length(l1_bound_energies),length(l2_free_energies));
for state_one = 1:length(l1_bound_energies)
    for state_two = 1:length(l0_free_energies)
        two_photon_dipoles_l0(state_one,state_two) = sum(l1_bound_states(:,state_one).*grid.*l0_free_states(:,state_two)).*sqrt(3).*Wigner3j(0,1,1,0,0,0).^2;
    end
    for state_two = 1:length(l2_free_energies)
        two_photon_dipoles_l2(state_one,state_two) = sum(l1_bound_states(:,state_one).*grid.*l2_free_states(:,state_two)).*sqrt(15).*Wigner3j(2,1,1,0,0,0).^2;
    end
    two_photon_dipoles_l1(state_one) = sum(l1_bound_states(:,state_one).*grid.*initial_state).*sqrt(3).*Wigner3j(1,1,0,0,0,0).^2;
end

%%
%===== Prepare Known Waveforms ==============
if gaussian_expand
    gaussian_trains = {}; k = 1;
    for harm = unknown_harmonics
        [FFT,t] = filterHarmonic(harm);
        mask = -1000 < t & t < 1000;
        time = t(mask);
        FFT = FFT(mask)*(max_intensity / max(FFT(mask)));
        gaussian_trains{k} = gaussianExpansion(time,FFT,gaussian_basis_size,harm);
        k = k + 1;
    end
else
    load("./data/experimental_pulses_7g.mat");
end

%%
%===== Prepare 'Experiment' Laser ===========
experiment_components = [];
for harm = unknown_harmonics
    experiment_components = [experiment_components, 'gaussian_train_', num2str(harm), ';'];
end
eval(['experiment = [',experiment_components(1:end-1),'];']);
if cross_correlation
    experiment_components = [];
    for harm = known_harmonics
        experiment_components = [experiment_components, 'gaussian_train_', num2str(harm), ';'];
    end
    eval(['known_laser = [',experiment_components(1:end-1),'];']);
end
rescale_factor = max_intensity / max(abs(experiment.calculate(time)));
temp_params = experiment.params();
temp_params(:,3:4) = temp_params(:,3:4) .* rescale_factor;
experiment = Laser.generate(temp_params,true,false);

%%
%===== Reconstruction Parameters ============
omega = Laser.SI2au_wavelength(800) * unknown_harmonics;
N_windows = 10; percentages = [2 5 10 15 25 50 75 100];
window_width = max(correlation_delay)/N_windows;
max_omega = max(omega) * 4;
max_frequency = max_omega / (2*pi);
sample_step = 0.5 / max_frequency;
options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-10,...
    'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'MaxFunctionEvaluations',1e4,'MaxIterations',5000,'FiniteDifferenceType', ...
    'forward','UseParallel',true,'Display','iter');
%%
%===== Reconstruction Functions =============
if cross_correlation
    calc = @(basis,delay) matrixElementsCalculation_xcorr(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,known_laser.params(),delay,omega,chirp);
else
    calc = @(basis,delay) matrixElementsCalculation(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,delay,omega,chirp);
end

%===== Generate Data to Reconstruct =========
known = calc(experiment.params(),correlation_delay);

%%
%===== Reconstruct Data =====================
filter = @(time) interp1(correlation_delay,known,time);
guesses = cell(size(reconstruction_gaussian_list)); ind = 1;
for N_gaussians = reconstruction_gaussian_list
    if N_gaussians > 1
        single_pulse = false;
    else
        single_pulse = true;
    end
    if length(unknown_harmonics) > 1
        single_color = false;
    else
        single_color = true;
    end
    initial_guess = []; gaussian_spacing = abs(max(correlation_delay)-min(correlation_delay))/N_gaussians/2;
    for i = 0:N_gaussians-1
        initial_guess = [initial_guess; Laser(1.8e-3*(1i),0.5,300,1e-4,i*gaussian_spacing,0);];
    end
    guess = initial_guess.params(single_pulse,single_color,chirp);
    lower_bound = ones(N_gaussians,1) * Laser(1e-4 - 100i,0.2,1,-1,-max(correlation_delay)).params(single_pulse,single_color,chirp);
    upper_bound = ones(N_gaussians,1) * Laser(10 + 100i,1.5,1000,1,max(correlation_delay)).params(single_pulse,single_color,chirp);
    
   for percent = percentages
        disp(percent)
        sample_tau = 0:sample_step:percent*window_width/100;
        tau = sample_tau;
        for window = 1:N_windows-1
            tau = [tau sample_tau + window*window_width];
        end
        sample_data = filter(tau)';
    
        err = @(basis) (calc(basis,tau) - sample_data)./max(abs(sample_data));
    
        new_guess = lsqnonlin(err,guess,lower_bound,upper_bound,options);
        guess = new_guess;
    end
    guesses{ind} = Laser.generate(guess,chirp,single_pulse,omega); ind = ind + 1;
end

filename = ['fit_harm' strjoin(cellstr(num2str(unknown_harmonics','%02d')),'') '_Npulses' strjoin(cellstr(num2str(reconstruction_gaussian_list')),'') '_chirp' num2str(chirp) '.mat'];
save(filename,'guesses','experiment');