%TODO
% -One gaussian Helium 9th harmonic
% -Two gaussian Helium 9th harmonic
% -One gaussian Helium 11th harmonic
% -Two gaussian Helium 11th harmonic
% -One gaussian Hydrogen 9th+11th harmonic
% -Two gaussian Hydrogen 9th+11th harmonic
%%
%===== Clear Previous Data ==================
close all; clear all;
%%
%===== Runtime Parameters ===================
%========== Comparison Basis Parameters =====
gaussian_expand = false; gaussian_basis_size = 16;
%========== Laser Parameters ================
max_intensity = 3e-3;
known_harmonics = []; unknown_harmonics = [9 11];
tmin = -1000; tmax = 1000; time = linspace(tmin,tmax,floor(2*tmax)+1);
%========== Reconstruction Parameters =======
reconstruction_gaussian_list = [1] * size(unknown_harmonics,2);
chirp = false; fit_color = false; cross_correlation = false;
checkpoint = true;
N_windows = 10; percentages = [5 10 15 25 50 75 100];
N_iterations = 500; N_epochs = 2;
options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-16,...
    'StepTolerance',1e-16,'OptimalityTolerance',1e-16,...
    'MaxFunctionEvaluations',N_iterations*20,'MaxIterations',N_iterations,...
    'FiniteDifferenceType','forward','UseParallel',true,'Display','iter');
%========== File and Directory Names ========
data_dir = './data/';
filename = ['fit_harm' strjoin(cellstr(num2str(unknown_harmonics','%02d')),'')...
            '_Npulses' strjoin(cellstr(num2str(reconstruction_gaussian_list')),'')...
            '_chirp' num2str(chirp) '.mat'];
%========== Overload Parameters =============
overload_guess = true;
overload_file = 'fit_harm0911_Npulses2_chirp0.mat';
overload_var = 'guesses';
%========== Setup Delay Sweep ===============
tau_max = 1000; dtau = 0.2;
if cross_correlation
	correlation_delay = linspace(-tau_max,tau_max,2*tau_max/dtau+1);
else % Autocorrelations are symmetric
	correlation_delay = linspace(0,tau_max,tau_max/dtau+1);
end

%%
%===== Read Data ============================
N_states = 500;
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
%===== Prepare Known Waveforms ==============
if gaussian_expand
    gaussian_trains = {}; k = 1;
    for harm = unknown_harmonics
        %===== Select Individual Harmonic ===
        [FFT,t] = filterHarmonic([data_dir 'filtered_hhg.txt'],harm);
        mask = tmin < t & t < tmax;
        time = t(mask);
        %===== Rescale Selected Harmonics ===
        FFT = FFT(mask)*(max_intensity / max(FFT(mask)));
        %===== Expanded in Gaussian Basis ===
        gaussian_trains{k} = gaussianExpansion(time,FFT,gaussian_basis_size,harm);
        eval(['gaussian_train_' num2str(harm) '= gaussian_trains{k};']);
        k = k + 1;
    end % Loop over `unknown_harmonics`
    save('expanded_harmonics.mat','gaussian_trains');
else % Use precomputed data if available
    load([data_dir 'experimental_pulses_16g.mat']);
end

%%
%===== Prepare 'Experiment' Laser ===========
experiment_components = [];
for harm = unknown_harmonics
    experiment_components = [experiment_components, 'gaussian_train_', num2str(harm), ';'];
end % Loop over `unknown_harmonics`
eval(['experiment = [',experiment_components(1:end-1),'];']);
if cross_correlation
    experiment_components = [];
    for harm = known_harmonics
        experiment_components = [experiment_components, 'gaussian_train_', num2str(harm), ';'];
    end % Loop over `known_harmonics`
    eval(['known_laser = [',experiment_components(1:end-1),'];']);
end
%===== Rescale Laser Intensity ==============
rescale_factor = max_intensity / max(abs(experiment.calculate(time)));
temp_params = experiment.params();
temp_params(:,3:4) = temp_params(:,3:4) .* rescale_factor;
%===== Generate 'Experiment Laser' ==========
experiment = Laser.generate(temp_params,chirp,false);

%%
%===== Reconstruction Parameters ============
omega = Laser.SI2au_wavelength(800) * unknown_harmonics;
window_width = max(correlation_delay)/N_windows;
max_omega = max(omega) * 4;
max_frequency = max_omega / (2*pi);
sample_step = 0.5 / max_frequency;
fit_color = ~fit_color;

%%
%===== Correlation Functions ================
if cross_correlation
    calc = @(basis,delay,N_gaussians,del_pos) squeeze(sum(sum(abs( ...
        matrixElementsCalculation_xcorr(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,known_laser.params(),delay,omega,chirp,N_gaussians,del_pos)).^2,2),3));
else
    calc = @(basis,delay,N_gaussians,del_pos) squeeze(sum(sum(abs( ...
        matrixElementsCalculation(initial_energy, ...
        N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
        N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
        N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        basis,delay,omega,chirp,N_gaussians,del_pos)).^2,2),3));
end
%%
%===== Generate Data to Reconstruct =========
known = calc(experiment.params(),correlation_delay,size(experiment.params(),1),[]);

%%
%===== Overload Initial Guess ===============
if overload_guess
    load(overload_file,overload_var);
    laser = guesses{1};
end

%%
%===== Reconstruct Data =====================
filter = @(time) interp1(correlation_delay,known,time);
guesses = cell(size(reconstruction_gaussian_list));
ion_error = cell(size(reconstruction_gaussian_list));
field_error = cell(size(reconstruction_gaussian_list)); ind = 1;
try % Try-catch to dump results on error
    for epoch = 1:N_epochs
        for N_gaussians = reconstruction_gaussian_list
            %===== Generate Initial Guess ===========
            if epoch == 1 % Use "random" for first round
                if ~overload_guess % Use random by default
                    initial_guess = []; gaussian_spacing = abs(tmax - tmin) / N_gaussians / 4;
                    for i = 1:N_gaussians
                        initial_guess = [initial_guess; Laser(max_intensity,0.5,500,1e-4,(-1).^i*floor(i/2)*gaussian_spacing,0);];
                    end % Loop over N_gaussians
                else % Use a different initial guess if available
                    initial_guess = laser;
                end
            else % Use previous epoch's best for next round
                initial_guess = Laser.generate(guess,chirp,omega_list);
            end
            guess = initial_guess.params(fit_color,chirp);
            %===== Set Parameter Bounds =============
            lower_bound = ones(N_gaussians,1) * Laser(-10 - 100i,0.2,1,-1,tmin).params(fit_color,chirp);
            upper_bound = ones(N_gaussians,1) * Laser(10 + 100i,1.5,1000,1,tmax).params(fit_color,chirp);
            %===== Determine Unnecessary Parameters =
            if fit_color
                del_pos = [4 3] * size(guess,1) - size(guess,1) + 1;
            else
                del_pos = [5 4] * size(guess,1) - size(guess,1) + 1;
            end
            %===== Reshape Parameter Matrix =========
            guess = reshape(guess,1,[]);
            lower_bound = reshape(lower_bound,1,[]);
            upper_bound = reshape(upper_bound,1,[]);
            %===== Remove Unnecessary Parameters ====
            for pos = del_pos
                guess(pos) = [];
                lower_bound(pos) = [];
                upper_bound(pos) = [];
            end
            
            %===== Begin Fitting Procedure ==========
            for percent = percentages
                disp(['Percent of Data: ' num2str(percent) '%']);
                %===== Windowing Setup ==============
                sample_tau = 0:sample_step:percent*window_width/100;
                tau = sample_tau;
                for window = 1:N_windows-1
                    tau = [tau sample_tau + window*window_width];
                end
                %===== Window Data ==================
                sample_data = filter(tau)';
                
                %===== Objective Function Setup =====
                err = @(basis) abs(calc(basis,tau,N_gaussians,del_pos) - sample_data)./max(abs(sample_data));
                %===== Fit Sampling of Data =========
                new_guess = lsqnonlin(err,guess,lower_bound,upper_bound,options);
                guess = new_guess;
                if checkpoint % Save on each iteration
                    temp_guess = guess;
                    for pos = flip(del_pos)
                        temp_guess = [temp_guess(1:pos-1) 0 temp_guess(pos:end)];
                    end
                    temp_guess = reshape(temp_guess,N_gaussians,[]);
                    omega_list = [];
                    for i=1:length(omega)
                        omega_list = [omega_list omega(i)*ones(N_gaussians/size(unknown_harmonics,2),1)];
                    end
                    current_guess = Laser.generate(temp_guess,chirp,omega_list);
                    save(['checkpoint_' num2str(percent) 'percent.mat'],'current_guess');
                end % Checkpoint
            end % Loop over `percentages
            %===== Undo Removal of Parameters =======
            for pos = flip(del_pos)
                guess = [guess(1:pos-1) 0 guess(pos:end)];
            end
            %===== Reshape Parameter Matrix =========
            guess = reshape(guess,N_gaussians,[]);
        
            %===== Calculate Ionization Error =======
            err = @(basis) abs(calc(basis,correlation_delay,N_gaussians,[]) - known)./max(abs(known));
            %===== Save Error and Guess Parameters ==
            ion_error{ind} = sqrt(dtau * sum(err(guess).^2));
            omega_list = [];
            for i=1:length(omega)
                omega_list = [omega_list omega(i)*ones(N_gaussians/size(unknown_harmonics,2),1)];
            end
            guesses{ind} = Laser.generate(guess,chirp,omega_list);
            ind = ind + 1;
        end % Loop over `reconstruction_gaussian_list`
        save(filename,'guesses','ion_error');
    end
catch % Dump results on error
    save('calculation_dump.mat');
end % Try-catch block