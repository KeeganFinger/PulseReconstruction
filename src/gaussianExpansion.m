function basis = gaussianExpansion(time,data,N_gaussians,harmonic)

if ~(size(time) == size(data))
    return;
end

color = Laser.au2SI_wavelength(800)*harmonic;
known_color = true;
chirp = true;

N = length(data);
% basis = [];
% if N_gaussians == 1
%     position = floor(N/2);
%     basis = [Laser(data(position),color,500,0,time(position))];
% else
%     for position = floor(N/N_gaussians/N_gaussians):floor(N/N_gaussians):N
%         basis = [basis; Laser(abs(data(position)),color,500,0,time(position))];
%     end
% end
spacing = floor(N/N_gaussians/2); midpoint = floor(N/2); offset = 0;
basis = [Laser(data(midpoint),color,500,0,time(midpoint)+offset)];
for gaussian = 2:N_gaussians
    position = midpoint + (-1)^gaussian * spacing * floor(gaussian/2);
    basis = [basis; Laser(data(position),color,500,0,time(position)+offset)];
end

err = @(gaussians) (imag(Laser.generate(gaussians,chirp,false,color).calculate(time) - data) ...
    )./ max(abs(data));

options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-16,...
    'StepTolerance',1e-16,'OptimalityTolerance',1e-10,...
    'MaxFunctionEvaluations',1e6,'MaxIterations',50000,'FiniteDifferenceType', ...
    'forward','UseParallel',true ,'Display','iter-detailed');
lower_bound = ones(N_gaussians,1) * Laser(-1,0.3,50,-1,min(time)).params(false,known_color,chirp);
upper_bound = ones(N_gaussians,1) * Laser(1,1,1000,1,max(time)).params(false,known_color,chirp);
if(known_color)
    lower_bound(:,3) = -Inf; upper_bound(:,3) = Inf;
else
    lower_bound(:,4) = -Inf; upper_bound(:,4) = Inf;
end
basis = lsqnonlin(err, basis.params(false,known_color,chirp), lower_bound, upper_bound, options);

basis = Laser.generate(basis,chirp,false,color);

end

function Z = out2(gaussians,time,chirp,harmonic)
    color = Laser.SI2au_wavelength(800) * harmonic;
    [~,Z] = Laser.generate(gaussians,chirp,false,color).calculate(time);
end