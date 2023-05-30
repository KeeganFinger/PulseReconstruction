function A = matrixElementsTPI(initial_energy, ...
    N_states_middle,dipole_middle,middle_energies, ...
    N_states_final,dipole_final,final_energies, ...
    laser_parameters_one,laser_parameters_two,correlation_delay)

%===== Prepare Variables ====================
A = zeros(2,N_states_final,length(correlation_delay));
position_one = laser_parameters_one(5);
position_two = laser_parameters_two(5);
omega_one = laser_parameters_one(1);
omega_two = laser_parameters_two(1);

%===== Final State Summation ================
for end_state = 1:N_states_final
    final_energy = final_energies(end_state);
    if (final_energy < 3 * max(omega_one,omega_two) + initial_energy)
        E_fi = final_energy - initial_energy;
        A(1,end_state,1) = final_energy;

        term_one = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_one,laser_parameters_two, ...
            position_one,position_two) ...
            .* (1 + exp(1i.*E_fi.*correlation_delay));

        term_two = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_one,laser_parameters_two, ...
            position_one+correlation_delay,position_two) ...
            .* exp(1i.*omega_one.*correlation_delay);

        term_three = subMatrixElementsTPI(initial_energy, ...
            N_states_middle,dipole_middle,middle_energies, ...
            end_state,dipole_final,final_energy, ...
            laser_parameters_one,laser_parameters_two, ...
            position_one,position_two+correlation_delay) ...
            .* exp(1i.*omega_two.*correlation_delay);

        A(2,end_state,:) = term_one + term_two + term_three;
    end
end
end