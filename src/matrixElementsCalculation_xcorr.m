function crosscorrelation = matrixElementsCalculation_xcorr(initial_energy, ...
    N_free_states_l0,two_photon_dipoles_l0,l0_free_energies, ...
    N_free_states_l2,two_photon_dipoles_l2,l2_free_energies, ...
    N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies, ...
    N_free_states_l1,one_photon_dipoles_l1,l1_free_energies, ...
    unknown_laser_list,known_laser_list,correlation_delay,omega_list,chirp,...
    N_gaussians,del_pos)

for pos = flip(del_pos)
    laser_list = [laser_list(1:pos-1) 0 laser_list(pos:end)];
end
laser_list = reshape(laser_list,N_gaussians,[]);

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

% Single photon
for i = 1:size(unknown_laser_list,1)
    unknown_laser = unknown_laser_list(i,:);
    O_unknown = matrixElementsSPI_xcorr(initial_energy,...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies,...
        unknown_laser,zeros(size(correlation_delay)));
    A(:,:,2) = A(:,:,2) + [squeeze(O_unknown(2,:,:))' l1_padding];
end
for i = 1:size(known_laser_list,1)
    known_laser = known_laser_list(i,:);
    O_known = matrixElementsSPI_xcorr(initial_energy,...
        N_free_states_l1,one_photon_dipoles_l1,l1_free_energies,...
        known_laser,correlation_delay);
    O_known(2,:,:) = squeeze(O_known(2,:,:)) .* exp(1i.*ones(size(O_known,2),1)*known_laser(1)*correlation_delay);
    A(:,:,2) = A(:,:,2) + [squeeze(O_known(2,:,:))' l1_padding];
end

% Double photon
for i = 1:size(unknown_laser_list,1)
    laser_one = unknown_laser_list(i,:);
    for j = 1:size(unknown_laser_list,1)
        laser_two = unknown_laser_list(j,:);
        T0_ff = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),zeros(size(correlation_delay)));
        A(:,:,1) = A(:,:,1) + [squeeze(T0_ff(2,:,:))' l0_padding];
        T2_ff = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),zeros(size(correlation_delay)));
        A(:,:,3) = A(:,:,3) + [squeeze(T2_ff(2,:,:))' l2_padding];
    end
end
for i = 1:size(unknown_laser_list,1)
    laser_one = unknown_laser_list(i,:);
    for j = 1:size(known_laser_list,1)
        laser_two = known_laser_list(j,:);
        T0_fg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),correlation_delay);
        T0_fg(2,:,:) = squeeze(T0_fg(2,:,:)) .* exp(1i.*ones(size(T0_fg,2),1)*laser_two(1)*correlation_delay);
        A(:,:,1) = A(:,:,1) + [squeeze(T0_fg(2,:,:))' l0_padding];
        T2_fg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            zeros(size(correlation_delay)),correlation_delay);
        T2_fg(2,:,:) = squeeze(T2_fg(2,:,:)) .* exp(1i.*ones(size(T2_fg,2),1)*laser_two(1)*correlation_delay);
        A(:,:,3) = A(:,:,3) + [squeeze(T2_fg(2,:,:))' l2_padding];
    end
end
for i = 1:size(known_laser_list,1)
    laser_one = known_laser_list(i,:);
    for j = 1:size(unknown_laser_list,1)
        laser_two = unknown_laser_list(j,:);
        T0_gf = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            correlation_delay,zeros(size(correlation_delay)));
        T0_gf(2,:,:) = squeeze(T0_gf(2,:,:)) .* exp(1i.*ones(size(T0_gf,2),1)*laser_one(1)*correlation_delay);
        A(:,:,1) = A(:,:,1) + [squeeze(T0_gf(2,:,:))' l0_padding];
        T2_gf = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            correlation_delay,zeros(size(correlation_delay)));
        T2_gf(2,:,:) = squeeze(T2_gf(2,:,:)) .* exp(1i.*ones(size(T2_gf,2),1)*laser_one(1)*correlation_delay);
        A(:,:,3) = A(:,:,3) + [squeeze(T2_gf(2,:,:))' l2_padding];
    end
end
for i = 1:size(known_laser_list,1)
    laser_one = known_laser_list(i,:);
    for j = 1:size(known_laser_list,1)
        laser_two = known_laser_list(j,:);
        T0_gg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l0,two_photon_dipoles_l0,l0_free_energies,...
            laser_one,laser_two,...
            correlation_delay,correlation_delay);
        T0_gg(2,:,:) = squeeze(T0_gg(2,:,:)) .* exp(1i.*ones(size(T0_gg,2),1)*(laser_one(1)+laser_two(1))*correlation_delay);
        A(:,:,1) = A(:,:,1) + [squeeze(T0_gg(2,:,:))' l0_padding];
        T2_gg = matrixElementsTPI_xcorr(initial_energy,...
            N_bound_states_l1,two_photon_dipoles_l1,l1_bound_energies,...
            N_free_states_l2,two_photon_dipoles_l2,l2_free_energies,...
            laser_one,laser_two,...
            correlation_delay,correlation_delay);
        T2_gg(2,:,:) = squeeze(T2_gg(2,:,:)) .* exp(1i.*ones(size(T2_gg,2),1)*(laser_one(1)+laser_two(1))*correlation_delay);
        A(:,:,3) = A(:,:,3) + [squeeze(T2_gg(2,:,:))' l2_padding];
    end
end

crosscorrelation = squeeze(sum(sum(abs(A).^2,2),3));
% crosscorrelation = squeeze(sum(abs(sum(A,2)).^2,3));
end