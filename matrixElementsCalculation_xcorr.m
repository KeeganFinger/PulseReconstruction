function autocorrelation = matrixElementsCalculation_xcorr(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    unknown_laser_list,known_laser_list,correlation_delay,omega_list,chirp)

color_override = [];
for color = omega_list
    for pulse = 1:size(unknown_laser_list,1)/size(omega_list,2)
        color_override = [color_override; color];
    end
end

% Fix data structure to include any parameters excluded from fitting space
if(chirp)
    if size(unknown_laser_list,2) == 4 % Single pulse of known color
        unknown_laser_list = [color_override unknown_laser_list(:,1:end-1) zeros(size(unknown_laser_list,1),1) unknown_laser_list(:,end)];
    elseif size(unknown_laser_list,2) == 5 % Single pulse of unknown color
        if(size(unknown_laser_list,1) > 1)
            unknown_laser_list = [color_override unknown_laser_list(:,:)];
        else
            unknown_laser_list = [unknown_laser_list(:,1:end-1) zeros(size(unknown_laser_list,1),1) unknown_laser_list(:,end)];
        end
    end 
else
    if size(unknown_laser_list,2) == 3 % Single pulse of known color
        unknown_laser_list = [color_override unknown_laser_list zeros(size(unknown_laser_list,1),1)];
    elseif size(unknown_laser_list,2) == 4 % Single pulse of unknown color
        if(size(unknown_laser_list,1) > 1)
            unknown_laser_list = [color_override unknown_laser_list(:,:)];
        else
            unknown_laser_list = [unknown_laser_list(:,:) zeros(size(unknown_laser_list,1),1)];
        end
    end % Default is multuple pulses of multiple unknown colors
    unknown_laser_list = [unknown_laser_list zeros(size(unknown_laser_list,1),1)];
end


max_states = max([N_free_states_l0,N_free_states_l1,N_free_states_l2]);
A = zeros(size(correlation_delay,2), ...
    max_states,3);

l0_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l0);
l1_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l1);
l2_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l2);

% Iterate over each laser in basis for first gaussian
for i = 1:size(unknown_laser_list,1)
for m = 1:size(known_laser_list,1)
    unknown_laser_one = unknown_laser_list(i,:);
    known_laser_one = known_laser_list(m,:);
    % Iterate over each laser in basis for second gaussian
    for j = 1:size(unknown_laser_list,1)
    for n = 1:size(known_laser_list,1)
        unknown_laser_two = unknown_laser_list(j,:);
        known_laser_two = known_laser_list(n,:);
        % Calculate the two-photon transition that ends in an L=0 state
        T0 = matrixElementsTPI_xcorr(initial_energy, ...
                N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
                N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
                unknown_laser_one,unknown_laser_two, ...
                known_laser_one,known_laser_two,correlation_delay);
        A(:,:,1) = A(:,:,1) + [squeeze(T0(2,:,:))' l0_padding];
        % Calculate the two-photon transition that ends in an L=2 state
        T2 = matrixElementsTPI_xcorr(initial_energy, ...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
            unknown_laser_one,unknown_laser_two, ...
            known_laser_one,known_laser_two,correlation_delay);
        A(:,:,3) = A(:,:,3) + [squeeze(T2(2,:,:))' l2_padding];
    end
    end
    % Calculate the one-photon transition
    O = matrixElementsSPI_xcorr(initial_energy, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        unknown_laser_one,known_laser_one,correlation_delay);
    A(:,:,2) = A(:,:,2) + [squeeze(O(2,:,:))' l1_padding];
end
end
autocorrelation = squeeze(sum(sum(abs(A).^2,2),3));
end