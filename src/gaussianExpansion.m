function basis = gaussianExpansion(time,data,N_gaussians,harmonic)

if ~(size(time) == size(data))
    return;
end

color = Laser.au2SI_wavelength(800)*harmonic;
known_color = true;
chirp = true;

N = length(data);

spacing = floor(N/N_gaussians/2); midpoint = floor(N/2); offset = 0;
basis = [Laser(real(data(midpoint)),color,500,0,0)];
for gaussian = 2:N_gaussians
    position = midpoint + (-1)^gaussian * spacing * floor(gaussian/2);
    basis = [basis; Laser(data(position),color,500,0,time(position)+offset)];
end

options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-16,...
    'StepTolerance',1e-16,'OptimalityTolerance',1e-10,...
    'MaxFunctionEvaluations',1e7,'MaxIterations',10000,'FiniteDifferenceType', ...
    'forward','UseParallel',true ,'Display','iter-detailed');
lower_bound = ones(N_gaussians,1) * Laser(-1,0.3,50,-1,min(time)).params(known_color,chirp);
upper_bound = ones(N_gaussians,1) * Laser(1,1,1000,1,max(time)).params(known_color,chirp);
basis = basis.params(known_color,chirp);
if(known_color)
    lower_bound(:,3) = -Inf; upper_bound(:,3) = Inf;
    del_pos = [4 3] * size(basis,1) - size(basis,1) + 1;
else
    lower_bound(:,4) = -Inf; upper_bound(:,4) = Inf;
    del_pos = [5 4] * size(basis,1) - size(basis,1) + 1;
end

basis = reshape(basis,1,[]);
lower_bound = reshape(lower_bound,1,[]);
upper_bound = reshape(upper_bound,1,[]);
for pos = del_pos
    basis(pos) = [];
    lower_bound(pos) = [];
    upper_bound(pos) = [];
end

err = @(gaussians) (abs(Laser.generate(fix_params(gaussians,del_pos,N_gaussians),chirp,color).calculate(time) - data) ...
    )./ max(abs(data));

basis = lsqnonlin(err, basis, lower_bound, upper_bound, options);

basis = Laser.generate(basis,chirp,color);

end

function L = fix_params(laser_list,del_pos,N_gaussians)
    for pos = flip(del_pos)
        laser_list = [laser_list(1:pos-1) 0 laser_list(pos:end)];
    end
    L = reshape(laser_list,N_gaussians,[]);
end