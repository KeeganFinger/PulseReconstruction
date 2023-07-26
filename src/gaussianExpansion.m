function basis = gaussianExpansion(time,data,N_gaussians,harmonic)

if ~(size(time) == size(data))
    return;
end

color = Laser.au2SI_wavelength(800)*harmonic;
known_color = true;
chirp = true;

N = length(data);

spacing = floor(length(time)/N_gaussians);
basis = [];
for gaussian = 1:N_gaussians
    position = spacing * gaussian;
    basis = [basis; Laser(data(position),color,300,0,time(position))];
end

options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-16,...
    'StepTolerance',1e-16,'OptimalityTolerance',1e-10,...
    'MaxFunctionEvaluations',1e7,'MaxIterations',500,'FiniteDifferenceType', ...
    'forward','UseParallel',true ,'Display','iter-detailed');
lower_bound = ones(N_gaussians,1) * Laser(-1,0.3,50,-1,min(time)).params(known_color,chirp);
upper_bound = ones(N_gaussians,1) * Laser(1,1,1000,1,max(time)).params(known_color,chirp);
basis = basis.params(known_color,chirp);
if(known_color)
    lower_bound(:,3) = -Inf; upper_bound(:,3) = Inf;
else
    lower_bound(:,4) = -Inf; upper_bound(:,4) = Inf;
end

err_imag = @(gaussians) (imag(Laser.generate(gaussians,chirp,color*ones(size(gaussians,1),1)).calculate(time) - data) ...
    )./ max(abs(data));

basis = lsqnonlin(err_imag, basis, lower_bound, upper_bound, options);

err_abs = @(gaussians) (abs(Laser.generate(gaussians,chirp,color*ones(size(gaussians,1),1)).calculate(time) - data) ...
    )./ max(abs(data));

basis = lsqnonlin(err_abs, basis, lower_bound, upper_bound, options);

basis = Laser.generate(basis,chirp,color*ones(size(basis,1),1));

end
