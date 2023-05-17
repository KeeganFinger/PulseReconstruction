close all; clear all;

%%
%===== Read Data ============================
N_states = 500;
for l=0:2
    data = h5read("He.h5",char("/Energy_l"+l));
    if l==0
        nmax = length(data(:,1));
        energies = zeros(nmax,3);
    end
    energies(1:nmax-l,l+1) = data(:,1) + 1i*data(:,2);
    data = h5read("He.h5",char("/Psi_r_l"+l));
    if l==0
        rmax = length(data(:,1,1));
        states = zeros(rmax,nmax,3);
    end
    states(:,1:nmax-l,l+1) = data(:,:,1) + 1i*data(:,:,2);
end

grid = h5read("parameters.h5",char("/EPS/r"));
dr = h5read("parameters.h5",char("/EPS/delta_x"));

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
load('C:\Users\finge\MATLAB\Projects\PulseReconstruction\Data\macroscopic_7g_9th.mat','t','FFT');
mask = -1000 < t & t < 1000;
time = t(mask); 
pulse_9 = FFT(mask);
load('C:\Users\finge\MATLAB\Projects\PulseReconstruction\Data\macroscopic_7g_11th.mat','FFT','t');
pulse_11 = FFT(mask);
max_intensity = 3e-3; %Laser.SI2au_intensity(1e14);
pulse_9 = (max_intensity / max(pulse_9)) * pulse_9;
pulse_11 = (max_intensity / max(pulse_11)) * pulse_11;
gaussian_train_9 = gaussianExpansion(time,conj(pulse_9),5,9);
gaussian_train_11 = gaussianExpansion(time,conj(pulse_11),5,11);
% save('experimental_pulses_5g.mat','time','pulse_9','pulse_11','gaussian_train_9','gaussian_train_11');

%%
load('experimental_pulses_5g.mat');
experiment = [gaussian_train_11];%; gaussian_train_11];
gaussian_train = experiment.params();

%%
%===== Reconstruction Parameters ============
single_color = true; chirp = false; xcorr = true;
harmonics = [11]; omega = Laser.SI2au_wavelength(800) * harmonics;
percentages = [2 5 10 15 25 50 75 100];
N_windows = 10;
correlation_delay = linspace(-1000,1000,2000);
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
if xcorr
    known_laser = gaussian_train_9.params();
    calc = @(basis,delay) matrixElementsCalculation_xcorr(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,known_laser,delay,omega,chirp);
else
    calc = @(basis,delay) matrixElementsCalculation(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,delay,omega,chirp);
end

%===== Generate Data to Reconstruct =========
known = calc(gaussian_train,correlation_delay);

%%
% options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-10,...
%     'StepTolerance',1e-10,'OptimalityTolerance',1e-10,...
%     'MaxFunctionEvaluations',1e4,'MaxIterations',50,'FiniteDifferenceType', ...
%     'forward','UseParallel',true,'Display','iter');
% initial_guess = [Laser(1.6e-3,0.5,300,1e-4,-200,0); Laser(1.6e-3,0.5,300,1e-4,400,0)];
% guess = initial_guess.params(single_pulse,single_color,chirp);
% lower_bound = ones(N_gaussians,1) * Laser(1e-4 - 100i,0.2,1,-1,-max(correlation_delay)).params(single_pulse,single_color,chirp);
% upper_bound = ones(N_gaussians,1) * Laser(10 + 100i,1.5,1000,1,max(correlation_delay)).params(single_pulse,single_color,chirp);
% err = @(basis) (calc(basis,correlation_delay) - known)./max(abs(known));
% new_guess = lsqnonlin(err,guess,lower_bound,upper_bound,options);

%%
%===== Reconstruct Data =====================
filter = @(time) interp1(correlation_delay,known,time);
gaussian_list = [2 4]; ind = 1;
guesses = cell(size(gaussian_list));
for N_gaussians = gaussian_list
    if N_gaussians > 1
        single_pulse = false;
    else
        single_pulse = true;
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

filename = ['fit_harm' strjoin(cellstr(num2str(harmonics','%02d')),'') '_Npulses' strjoin(cellstr(num2str(gaussian_list')),'') '_chirp' num2str(chirp) '.mat'];
save(filename,'guesses','experiment');