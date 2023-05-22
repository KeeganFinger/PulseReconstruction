function A = matrixElementsSPI_xcorr(initial_energy, ...
    N_states,dipole,final_energies, ...
    laser_parameters,correlation_delay)

%===== Final State Summation ================
A = zeros(2,N_states,length(correlation_delay));
position = laser_parameters(5);
for end_state = 1:N_states
    final_energy = final_energies(end_state);
    A(1,end_state,1) = final_energy;

    A(2,end_state,:) = subMatrixElementsSPI(dipole,initial_energy, ...
        final_energy,end_state,laser_parameters,position+correlation_delay);
end
end