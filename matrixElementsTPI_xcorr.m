function A = matrixElementsTPI_xcorr(initial_energy, ...
    N_states_middle,dipole_middle,middle_energies, ...
    N_states_final,dipole_final,final_energies, ...
    laser_parameters_one,laser_parameters_two, ...
    laser_parameters_three,laser_parameters_four,correlation_delay)

%===== Prepare Variables ====================
A = zeros(2,N_states_final,length(correlation_delay));
position_one = laser_parameters_one(5); % undelayed unknown
position_two = laser_parameters_two(5); % delayed unknown
position_three = laser_parameters_three(5); %undelayed known
position_four = laser_parameters_four(5); %delayed known
omega_one = laser_parameters_one(1);
omega_two = laser_parameters_two(1);
omega_three = laser_parameters_three(1);
omega_four = laser_parameters_four(1);

%===== Final State Summation ================
for end_state = 1:N_states_final
    final_energy = final_energies(end_state);
    if (final_energy < 3 * max([omega_one,omega_two,omega_three,omega_four]) + initial_energy)
        E_fi = final_energy - initial_energy;
        A(1,end_state,1) = final_energy;

        % Undelayed laser from each pulse
        term_one = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_three,laser_parameters_one, ...
            position_three,position_one);
           
        % Delayed laser from known pulse and undelayed from unknown pulse
        term_two = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_four,laser_parameters_one, ...
            position_four+correlation_delay,position_one) ...
            .* exp(1i.*omega_four.*correlation_delay);

        % Undelayed laser from known pulse and delayed from unknown pulse
        term_three = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_three,laser_parameters_two, ...
            position_three,position_two+correlation_delay) ...
            .* exp(1i.*omega_two.*correlation_delay);

        % Delayed laser from known pulse and delayed from unknown pulse
        term_four = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_four,laser_parameters_two, ...
            position_four,position_two) ...
            .* exp(1i.*E_fi.*correlation_delay);

        A(2,end_state,:) = term_one + term_two + term_three + term_four;
    end
end
end

function S = subMatrixElementsTPI(initial_energy, ...
    N_states,dipole_middle,middle_energies, ...
    final_state,dipole_final,final_energy, ...
    laser_parameters_one,laser_parameters_two, ...
    position_one,position_two)

%===== Laser One Parameters =================
frequency_one = laser_parameters_one(1);
period_one = 1/sqrt(1/laser_parameters_one(2)^2 + 1i*laser_parameters_one(6));
amplitude_one = laser_parameters_one(3) + 1i*laser_parameters_one(4);

%===== Laser Two Parameters =================
frequency_two = laser_parameters_two(1);
period_two = 1/sqrt(1/laser_parameters_two(2)^2 + 1i*laser_parameters_two(6));
amplitude_two = laser_parameters_two(3) + 1i*laser_parameters_two(4);

%===== Matrix Element =======================
S=0;
for mid_state = 1:N_states
    middle_energy = middle_energies(mid_state);
    delta_fm = final_energy - middle_energy - frequency_two;
    delta_mi = middle_energy - initial_energy - frequency_one;

    prefactor = -0.25 * pi * amplitude_two * amplitude_one * period_two * period_one;
% when the amplitudes are abs then the spectrum is good
%     prefactor = 1;
    
    laser_factor_real = -0.5.*(period_one.^2 .* delta_mi.^2 + period_two.^2 .* delta_fm.^2);
    laser_factor_imag = 1i.*(position_two.*delta_fm + position_one.*delta_mi);
    erf_factor = 1./sqrt(2) .* 1./sqrt(period_one.^2 + period_two^2) .* (position_two - position_one + 1i.*(delta_fm.*period_two^2 - delta_mi.*period_one^2));
    laser1 = exp(laser_factor_imag);
    laser2 = exp(laser_factor_real);
    laser3 = error_function(erf_factor,laser_factor_real);
    laser = laser1 .* (laser2 + laser3);

    S = S + (prefactor .* dipole_final(mid_state,final_state) .* dipole_middle(mid_state) .* laser);
end
end
