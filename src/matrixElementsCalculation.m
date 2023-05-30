function autocorrelation = matrixElementsCalculation(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    laser_list,correlation_delay,omega_list,chirp)

color_override = [];
for color = omega_list
    for pulse = 1:size(laser_list,1)/size(omega_list,2)
        color_override = [color_override; color];
    end
end

% Fix data structure to include any parameters excluded from fitting space
if(chirp)
    if size(laser_list,2) == 4 % Single pulse of known color
        laser_list = [color_override laser_list(:,1:end-1) zeros(size(laser_list,1),1) laser_list(:,end)];
    elseif size(laser_list,2) == 5 % Single pulse of unknown color
        if(size(laser_list,1) > 1)
            laser_list = [color_override laser_list(:,:)];
        else
            laser_list = [laser_list(:,1:end-1) zeros(size(laser_list,1),1) laser_list(:,end)];
        end
    end 
else
    if size(laser_list,2) == 3 % Single pulse of known color
        laser_list = [color_override laser_list zeros(size(laser_list,1),1)];
    elseif size(laser_list,2) == 4 % Single pulse of unknown color
        if(size(laser_list,1) > 1)
            laser_list = [color_override laser_list(:,:)];
        else
            laser_list = [laser_list(:,:) zeros(size(laser_list,1),1)];
        end
    end % Default is multuple pulses of multiple unknown colors
    laser_list = [laser_list zeros(size(laser_list,1),1)];
end


max_states = max([N_free_states_l0,N_free_states_l1,N_free_states_l2]);
A = zeros(size(correlation_delay,2), ...
    max_states,3);

l0_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l0);
l1_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l1);
l2_padding = zeros(size(correlation_delay,2),max_states-N_free_states_l2);

% Iterate over each laser in basis for first gaussian
for i = 1:size(laser_list,1)
    laser_one = laser_list(i,:);
    % Iterate over each laser in basis for second gaussian
    for j = 1:size(laser_list,1)
        laser_two = laser_list(j,:);
        % Calculate the two-photon transition that ends in an L=0 state
        T0 = matrixElementsTPI(initial_energy, ...
                N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
                N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
                laser_one,laser_two,correlation_delay);
        A(:,:,1) = A(:,:,1) + [squeeze(T0(2,:,:))' l0_padding];
        % Calculate the two-photon transition that ends in an L=2 state
        T2 = matrixElementsTPI(initial_energy, ...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
            laser_one,laser_two,correlation_delay);
        A(:,:,3) = A(:,:,3) + [squeeze(T2(2,:,:))' l2_padding];
    end
    % Calculate the one-photon transition
    O = matrixElementsSPI(initial_energy, ...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
        laser_one,correlation_delay);
    A(:,:,2) = A(:,:,2) + [squeeze(O(2,:,:))' l1_padding];
end
autocorrelation = squeeze(sum(sum(abs(A).^2,2),3));
end