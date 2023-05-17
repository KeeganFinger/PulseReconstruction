function A = matrixElementsSPI(initial_energy, ...
    N_states,dipole,final_energies, ...
    laser_parameters,correlation_delay)

%===== Final State Summation ================
A = zeros(2,N_states,length(correlation_delay));
position = laser_parameters(5);
for end_state = 1:N_states
    final_energy = final_energies(end_state);
    E_fi = final_energy - initial_energy;
    A(1,end_state,1) = final_energy;

    A(2,end_state,:) = subMatrixElementsSPI(dipole,initial_energy, ...
        final_energy,end_state,laser_parameters,position) ...
        .* (1 + exp(1i.*E_fi.*correlation_delay));
end

end

function S = subMatrixElementsSPI(dipole,initial_energy,final_energy, ...
    final_state,laser_parameters,position)

%===== Laser Parameters =====================
frequency = laser_parameters(1);
period = 1/sqrt(1/laser_parameters(2)^2 + 1i*laser_parameters(6));
amplitude = laser_parameters(3) + 1i*laser_parameters(4);

%===== Matrix Element =======================
% prefactor = -1i .* sqrt(pi/2) .* amplitude .* period; %does this period get affected by the chirp?
prefactor = -1i .* sqrt(pi/2) .* amplitude .* laser_parameters(2);

delta_fi = final_energy - initial_energy - frequency;
laser_factor = -1i .* delta_fi .* position - 0.5 .* (delta_fi.^2 .* period.^2);
laser = exp(laser_factor);

S = prefactor .* dipole(final_state,1,1) .* laser;

end