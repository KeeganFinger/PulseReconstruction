function [harmonic,time] = filterHarmonic(harm)

load("Multiphoton\500000c_0um.mat","angle_saves");
harmonic_order = angle_saves.H;
spectrum = angle_saves.E/max(abs(angle_saves.E));

%===== Filter Selected Harmonic =============
beta = 0.35;
filter = exp(-0.5*(harmonic_order-harm).^2/(2*beta.^2));
spectral_data = filter .* spectrum;

%===== Inverse Fourier Transfrom ============
fundamental_frequency = Laser.SI2au_wavelength(800);
padding = 2^15;
spectral_data = [blackman(length(spectral_data)) .* spectral_data; zeros(padding,1)];
harmonic = fft(spectral_data);
harmonic = flip(harmonic(1:end)');
harmonic = harmonic/max(real(harmonic))*3e-3;

max_time = 2*pi/(harmonic_order(2)-harmonic_order(1))/fundamental_frequency;
time = 0:max_time/length(spectral_data):max_time;
time = time(1:end-1) - max(time)/2;
time = flip(-time);
end