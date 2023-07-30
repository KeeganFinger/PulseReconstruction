clear all; close all;
% load("./results/results_0612.mat");

time = linspace(-1000,1000,20001);
options = optimoptions(@lsqnonlin,'FunctionTolerance',1e-18,...
    'StepTolerance',1e-18,'OptimalityTolerance',1e-15,...
    'MaxFunctionEvaluations',1e7,'MaxIterations',50000,'FiniteDifferenceType', ...
    'forward','UseParallel',true ,'Display','iter-detailed');

% Unchirped 9th Harmonic
load("./results/new_results/fit_harm09_Npulses1_chirp0.mat")
harm9_exact = experiment.calculate(time);
fit_guess = @(params) abs(exp(1i*params(1)) * guesses{1}.calculate(time - params(2)) - harm9_exact) / (max(abs(harm9_exact)));
params = lsqnonlin(fit_guess, [0 0], [], [], options);
harm9_chirp0_1g_vals = exp(1i*params(1)) * guesses{1}.calculate(time - params(2));

load("./results/new_results/fit_harm09_Npulses24_chirp0.mat");
fit_guess = @(params) abs(exp(1i*params(1)) * guesses{2}.calculate(time - params(2)) - harm9_exact) / (max(abs(harm9_exact)));
params = lsqnonlin(fit_guess, [0 0], [], [], options);
harm9_chirp0_4g_vals = exp(1i*params(1)) * guesses{2}.calculate(time - params(2));

load("fit_harm09_Npulses2_chirp0.mat");
fit_guess = @(params) abs(exp(1i*params(1)) * flip(conj(guesses{1}.calculate(time - params(2)))) - harm9_exact) / (max(abs(harm9_exact)));
params = lsqnonlin(fit_guess, [0 0], [], [], options);
harm9_chirp0_2g_vals = exp(1i*params(1)) * flip(conj(guesses{1}.calculate(time - params(2))));

load("./results/fit_harm09_Npulses6_chirp0.mat");
harm9_chirp0_6g_vals = guesses{1}.calculate(time);

load("./results/short_results/fit_harm11_Npulses1246_chirp0.mat");
harm11_exact = experiment.calculate(time);
harm11_chirp0_1g_vals = guesses{1}.calculate(time);
harm11_chirp0_2g_vals = guesses{2}.calculate(time);
harm11_chirp0_4g_vals = guesses{3}.calculate(time);
harm11_chirp0_6g_vals = guesses{4}.calculate(time);

% load("./results/fit_harm11_Npulses1_chirp0.mat");
% harm11_exact = experiment.calculate(time);
% harm11_chirp0_1g_vals = guesses{1}.calculate(time);
% 
% load("./results/fit_harm11_Npulses2_chirp0.mat");
% % params = guesses{1}.params();
% % params(2,3:4) = params(2,3:4) * 1.07;
% % guesses{1} = Laser.generate(params,false,0);
% harm11_chirp0_2g_vals = guesses{1}.calculate(time);
% 
% load("./results/fit_harm11_Npulses24_chirp0.mat");
% harm11_chirp0_4g_vals = guesses{2}.calculate(time);
% 
% load("./results/fit_harm11_Npulses6_chirp0.mat");
% harm11_chirp0_6g_vals = guesses{1}.calculate(time);
% 
load("./results/fit_harm1113_Npulses2_chirp0.mat");
harm1113_exact = experiment.calculate(time);
harm1113_chirp0_1g_vals = guesses{1}.calculate(time);

% Chirped data
load("./results/new_results/fit_harm09_Npulses124_chirp1.mat");
fit_guess = @(params) imag(exp(1i*params(1)) * guesses{1}.calculate(time - params(2)) - harm9_exact) / (max(abs(harm9_exact)));
params = lsqnonlin(fit_guess, [0 -100], [], [], options);
harm9_chirp1_1g_vals = exp(1i*params(1)) * guesses{1}.calculate(time - params(2));
fit_guess = @(params) imag(exp(1i*params(1)) * flip(guesses{2}.calculate(time - params(2))) - harm9_exact) / (max(abs(harm9_exact)));
params = lsqnonlin(fit_guess, [0 -300], [], [], options);
harm9_chirp1_2g_vals = exp(1i*params(1)) * flip(guesses{2}.calculate(time - params(2)));
fit_guess = @(params) abs(exp(1i*params(1)) * guesses{3}.calculate(time - params(2)) - harm9_exact) / (max(abs(harm9_exact)));
params = lsqnonlin(fit_guess, [0 0], [], [], options);
harm9_chirp1_4g_vals = exp(1i*params(1)) * guesses{3}.calculate(time - params(2));

load("./results/fit_harm09_Npulses46_chirp1.mat");
% harm9_chirp1_4g_vals = guesses{1}.calculate(time);
harm9_chirp1_6g_vals = guesses{2}.calculate(time);

load("./results/fit_harm11_Npulses12_chirp1.mat");
harm11_chirp1_1g_vals = guesses{1}.calculate(time);
harm11_chirp1_2g_vals = guesses{2}.calculate(time);


%%
left_ylabel = "$$|\tilde{f}(t)|$$ (a.u.)";
right_ylabel = "$$\omega(t)$$ (a.u.)";

left_ycolor = "#0072BD";
right_ycolor = "#A2142F";

left_ylim_9 = [0 3.2e-3];
right_ylim_9 = [0.47 0.55];
x_lim_9 = [-700 900];

left_ylim_11 = [0 3.2e-3];
right_ylim_11 = [0.59 0.67];
x_lim_11 = [-700 900];

left_ylim_1113 = [-6e-3 6e-3];
right_ylim_1113 = [0 1];
x_lim_1113 = [-700 900];

dt = abs(time(2)-time(1));

%%
% 9th Harmonic Autocorrelation
harm9_exp_phase = unwrap(angle(conj(harm9_exact)));
harm9_exp_omega = zeros(size(harm9_exp_phase));

harm9_chirp0_1g_phase = unwrap(angle(conj(harm9_chirp0_1g_vals)));
harm9_chirp0_1g_omega = zeros(size(harm9_chirp0_1g_phase));
harm9_chirp0_2g_phase = unwrap(angle(conj(harm9_chirp0_2g_vals)));
harm9_chirp0_2g_omega = zeros(size(harm9_chirp0_2g_phase));
harm9_chirp0_4g_phase = unwrap(angle(conj(harm9_chirp0_4g_vals)));
harm9_chirp0_4g_omega = zeros(size(harm9_chirp0_4g_phase));
harm9_chirp0_6g_phase = unwrap(angle(harm9_chirp0_6g_vals));
harm9_chirp0_6g_omega = zeros(size(harm9_chirp0_6g_phase));

harm9_chirp1_1g_phase = unwrap(angle(conj(harm9_chirp1_1g_vals)));
harm9_chirp1_1g_omega = zeros(size(harm9_chirp1_1g_phase));
harm9_chirp1_2g_phase = unwrap(angle(conj(harm9_chirp1_2g_vals)));
harm9_chirp1_2g_omega = zeros(size(harm9_chirp1_2g_phase));
harm9_chirp1_4g_phase = unwrap(angle(conj(harm9_chirp1_4g_vals)));
harm9_chirp1_4g_omega = zeros(size(harm9_chirp1_4g_phase));
harm9_chirp1_6g_phase = unwrap(angle(harm9_chirp1_6g_vals));
harm9_chirp1_6g_omega = zeros(size(harm9_chirp1_6g_phase));

%%
% 11th Harmonic Autocorrelation
harm11_exp_phase = unwrap(angle(conj(harm11_exact)));
harm11_exp_omega = zeros(size(harm11_exp_phase));

harm11_chirp0_1g_phase = unwrap(angle(conj(harm11_chirp0_1g_vals)));
harm11_chirp0_1g_omega = zeros(size(harm11_chirp0_1g_phase));
harm11_chirp0_2g_phase = unwrap(angle(conj(harm11_chirp0_2g_vals)));
harm11_chirp0_2g_omega = zeros(size(harm11_chirp0_2g_phase));
harm11_chirp0_4g_phase = unwrap(angle(conj(harm11_chirp0_4g_vals)));
harm11_chirp0_4g_omega = zeros(size(harm11_chirp0_4g_phase));
harm11_chirp0_6g_phase = unwrap(angle(harm11_chirp0_6g_vals));
harm11_chirp0_6g_omega = zeros(size(harm11_chirp0_6g_phase));

harm11_chirp1_1g_phase = unwrap(angle(conj(harm11_chirp1_1g_vals)));
harm11_chirp1_1g_omega = zeros(size(harm11_chirp1_1g_phase));
harm11_chirp1_2g_phase = unwrap(angle(conj(harm11_chirp1_2g_vals)));
harm11_chirp1_2g_omega = zeros(size(harm11_chirp1_2g_phase));
% harm11_chirp1_4g_phase = unwrap(angle(conj(harm11_chirp1_4g_vals)));
% harm11_chirp1_4g_omega = zeros(size(harm11_chirp1_4g_phase));
% harm11_chirp1_6g_phase = unwrap(angle(harm11_chirp1_6g_vals));
% harm11_chirp1_6g_omega = zeros(size(harm11_chirp1_6g_phase));

% 11th + 13th Autocorrelation
harm1113_exp_phase = unwrap(angle(conj(harm1113_exact)));
harm1113_exp_omega = zeros(size(harm1113_exp_phase));

harm1113_chirp0_1g_phase = unwrap(angle(conj(harm1113_chirp0_1g_vals)));
harm1113_chirp0_1g_omega = zeros(size(harm1113_chirp0_1g_phase));

for t = 2:length(time)-1
    harm9_exp_omega(t) = (harm9_exp_phase(t+1) - harm9_exp_phase(t-1))/(2 * dt);
    harm11_exp_omega(t) = (harm11_exp_phase(t+1) - harm11_exp_phase(t-1))/(2 * dt);
    harm1113_exp_omega(t) = (harm1113_exp_phase(t+1) - harm1113_exp_phase(t-1))/(2 * dt);
end

for t = 2:length(time)-1
    harm9_chirp0_1g_omega(t) = (harm9_chirp0_1g_phase(t+1) - harm9_chirp0_1g_phase(t-1))/(2 * dt);
    harm9_chirp0_2g_omega(t) = (harm9_chirp0_2g_phase(t+1) - harm9_chirp0_2g_phase(t-1))/(2 * dt);
    harm9_chirp0_4g_omega(t) = (harm9_chirp0_4g_phase(t+1) - harm9_chirp0_4g_phase(t-1))/(2 * dt);
    harm9_chirp0_6g_omega(t) = -(harm9_chirp0_6g_phase(t+1) - harm9_chirp0_6g_phase(t-1))/(2 * dt);

    harm9_chirp1_1g_omega(t) = (harm9_chirp1_1g_phase(t+1) - harm9_chirp1_1g_phase(t-1))/(2 * dt);
    harm9_chirp1_2g_omega(t) = -(harm9_chirp1_2g_phase(t+1) - harm9_chirp1_2g_phase(t-1))/(2 * dt);
    harm9_chirp1_4g_omega(t) = (harm9_chirp1_4g_phase(t+1) - harm9_chirp1_4g_phase(t-1))/(2 * dt);
    harm9_chirp1_6g_omega(t) = -(harm9_chirp1_6g_phase(t+1) - harm9_chirp1_6g_phase(t-1))/(2 * dt);


    harm11_chirp0_1g_omega(t) = (harm11_chirp0_1g_phase(t+1) - harm11_chirp0_1g_phase(t-1))/(2 * dt);
    harm11_chirp0_2g_omega(t) = (harm11_chirp0_2g_phase(t+1) - harm11_chirp0_2g_phase(t-1))/(2 * dt);
    harm11_chirp0_4g_omega(t) = (harm11_chirp0_4g_phase(t+1) - harm11_chirp0_4g_phase(t-1))/(2 * dt);
    harm11_chirp0_6g_omega(t) = -(harm11_chirp0_6g_phase(t+1) - harm11_chirp0_6g_phase(t-1))/(2 * dt);

    harm11_chirp1_1g_omega(t) = (harm11_chirp1_1g_phase(t+1) - harm11_chirp1_1g_phase(t-1))/(2 * dt);
    harm11_chirp1_2g_omega(t) = (harm11_chirp1_2g_phase(t+1) - harm11_chirp1_2g_phase(t-1))/(2 * dt);
%     harm11_chirp1_4g_omega(t) = (harm11_chirp1_4g_phase(t+1) - harm11_chirp1_4g_phase(t-1))/(2 * dt);
%     harm11_chirp1_6g_omega(t) = -(harm11_chirp1_6g_phase(t+1) - harm11_chirp1_6g_phase(t-1))/(2 * dt);

    harm1113_chirp0_1g_omega(t) = (harm1113_chirp0_1g_phase(t+1) - harm1113_chirp0_1g_phase(t-1))/(2 * dt);
end
%%
figure(); tiles_9_chirp0 = tiledlayout(2,2);
title(tiles_9_chirp0,"9th Harmonic Reconstruction - No Chirp",'Interpreter','latex');
tiles_9_chirp0.Padding = 'compact'; tiles_9_chirp0.TileSpacing = 'compact';
%%
ax1_9_chirp0 = nexttile; hold on;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1_9_chirp0,time,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax1_9_chirp0,time,abs(harm9_chirp0_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1_9_chirp0,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax1_9_chirp0,time,harm9_chirp0_1g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax2_9_chirp0 = nexttile; hold on;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp0,time,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax2_9_chirp0,time,abs(harm9_chirp0_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp0,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax2_9_chirp0,time,harm9_chirp0_2g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax3_9_chirp0 = nexttile; hold on;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp0,time,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax3_9_chirp0,time,abs(harm9_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp0,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax3_9_chirp0,time,harm9_chirp0_4g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax4_9_chirp0 = nexttile; hold on;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp0,time,abs(harm9_exact),'--','Color',left_ycolor);
% plot(ax4_9_chirp0,time,abs(harm9_chirp0_6g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp0,time,harm9_exp_omega,'--','Color',right_ycolor)
% plot(ax4_9_chirp0,time,harm9_chirp0_6g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
linkaxes([ax1_9_chirp0,ax2_9_chirp0],'y')
linkaxes([ax3_9_chirp0,ax4_9_chirp0],'y')
linkaxes([ax1_9_chirp0,ax3_9_chirp0],'x')
linkaxes([ax2_9_chirp0,ax4_9_chirp0],'x')

%%
figure(); tiles_9_chirp1 = tiledlayout(2,2);
title(tiles_9_chirp1,"9th Harmonic Reconstruction - Chirp",'Interpreter','latex');
tiles_9_chirp1.Padding = 'compact'; tiles_9_chirp1.TileSpacing = 'compact';
%%
ax1_9_chirp1 = nexttile; hold on;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1_9_chirp1,time,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax1_9_chirp1,time,abs(harm9_chirp1_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1_9_chirp1,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax1_9_chirp1,time,harm9_chirp1_1g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax2_9_chirp1 = nexttile; hold on;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp1,time,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax2_9_chirp1,time,abs(harm9_chirp1_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp1,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax2_9_chirp1,time,harm9_chirp1_2g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax3_9_chirp1 = nexttile; hold on;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp1,time,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax3_9_chirp1,time,abs(harm9_chirp1_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp1,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax3_9_chirp1,time,harm9_chirp1_4g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax4_9_chirp1 = nexttile; hold on;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp1,time,abs(harm9_exact),'--','Color',left_ycolor);
% plot(ax4_9_chirp1,time,flip(abs(harm9_chirp1_6g_vals)),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp1,time,harm9_exp_omega,'--','Color',right_ycolor)
% plot(ax4_9_chirp1,time,flip(harm9_chirp1_6g_omega),'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
linkaxes([ax1_9_chirp1,ax2_9_chirp1],'y')
linkaxes([ax3_9_chirp1,ax4_9_chirp1],'y')
linkaxes([ax1_9_chirp1,ax3_9_chirp1],'x')
linkaxes([ax2_9_chirp1,ax4_9_chirp1],'x')
% 
% %%
% figure(); tiles_11_chirp0 = tiledlayout(2,2);
% title(tiles_11_chirp0,"11th Harmonic Reconstruction - No Chirp",'Interpreter','latex');
% tiles_11_chirp0.Padding = 'compact'; tiles_11_chirp0.TileSpacing = 'compact';
% %%
% ax1_11_chirp0 = nexttile; hold on;
% title("1 Gaussian",'Interpreter','latex');
% yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax1_11_chirp0,time,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax1_11_chirp0,time,abs(harm11_chirp0_1g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
% plot(ax1_11_chirp0,time,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax1_11_chirp0,time,harm11_chirp0_1g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% hold off; ax2_11_chirp0 = nexttile; hold on;
% title("2 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax2_11_chirp0,time,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax2_11_chirp0,time,abs(flip(harm11_chirp0_2g_vals)),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax2_11_chirp0,time,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax2_11_chirp0,time,flip(harm11_chirp0_2g_omega),'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% hold off; ax3_11_chirp0 = nexttile; hold on;
% title("4 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax3_11_chirp0,time,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax3_11_chirp0,time,abs(harm11_chirp0_4g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax3_11_chirp0,time,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax3_11_chirp0,time,harm11_chirp0_4g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% hold off; ax4_11_chirp0 = nexttile; hold on;
% title("6 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax4_11_chirp0,time,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax4_11_chirp0,time,abs(harm11_chirp0_6g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax4_11_chirp0,time,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax4_11_chirp0,time,harm11_chirp0_6g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% linkaxes([ax1_11_chirp0,ax2_11_chirp0],'y')
% linkaxes([ax3_11_chirp0,ax4_11_chirp0],'y')
% linkaxes([ax1_11_chirp0,ax3_11_chirp0],'x')
% linkaxes([ax2_11_chirp0,ax4_11_chirp0],'x')
% 
% %%
% figure(); tiles_11_chirp1 = tiledlayout(2,2);
% title(tiles_11_chirp1,"11th Harmonic Reconstruction - Chirp",'Interpreter','latex');
% tiles_11_chirp1.Padding = 'compact'; tiles_11_chirp1.TileSpacing = 'compact';
% %%
% ax1_11_chirp1 = nexttile; hold on;
% title("1 Gaussian",'Interpreter','latex');
% yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax1_11_chirp1,time,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax1_11_chirp1,time,abs(harm11_chirp1_1g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
% plot(ax1_11_chirp1,time,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax1_11_chirp1,time,harm11_chirp1_1g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% hold off; ax2_11_chirp1 = nexttile; hold on;
% title("2 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax2_11_chirp1,time,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax2_11_chirp1,time,abs(harm11_chirp1_2g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax2_11_chirp1,time,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax2_11_chirp1,time,harm11_chirp1_2g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% hold off; ax3_11_chirp1 = nexttile; hold on;
% title("4 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax3_11_chirp1,time,abs(harm11_exact),'--','Color',left_ycolor);
% % plot(ax3_11_chirp1,time,abs(harm11_chirp1_4g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax3_11_chirp1,time,harm11_exp_omega,'--','Color',right_ycolor)
% % plot(ax3_11_chirp1,time,harm11_chirp1_4g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% hold off; ax4_11_chirp1 = nexttile; hold on;
% title("6 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax4_11_chirp1,time,abs(harm11_exact),'--','Color',left_ycolor);
% % plot(ax4_11,time,abs(harm11_chirp1_4g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax4_11_chirp1,time,harm11_exp_omega,'--','Color',right_ycolor)
% % plot(ax4_11,time,harm11_chirp1_6g_omega,'-','Color',right_ycolor)
% xlim(x_lim_11); ylim(right_ylim_11);
% %%
% linkaxes([ax1_11_chirp1,ax2_11_chirp1],'y')
% linkaxes([ax3_11_chirp1,ax4_11_chirp1],'y')
% linkaxes([ax1_11_chirp1,ax3_11_chirp1],'x')
% linkaxes([ax2_11_chirp1,ax4_11_chirp1],'x')
% %%
% figure(); tiles_1113_chirp0 = tiledlayout(2,2);
% title(tiles_1113_chirp0,"11th + 13th Harmonic Reconstruction - Chirp",'Interpreter','latex');
% tiles_1113_chirp0.Padding = 'compact'; tiles_1113_chirp0.TileSpacing = 'compact';
% %%
% ax1_1113_chirp0 = nexttile; hold on;
% title("1 Gaussian",'Interpreter','latex');
% yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
% ylim(left_ylim_1113); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax1_1113_chirp0,time,abs(harm1113_exact),'--','Color',left_ycolor);
% plot(ax1_1113_chirp0,time,abs(harm1113_chirp0_1g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
% plot(ax1_1113_chirp0,time,harm1113_exp_omega,'--','Color',right_ycolor)
% plot(ax1_1113_chirp0,time,harm1113_chirp0_1g_omega,'-','Color',right_ycolor)
% xlim(x_lim_1113); ylim(right_ylim_1113);
% %%
% hold off; ax2_1113_chirp0 = nexttile; hold on;
% title("2 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_1113); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax2_1113_chirp0,time,abs(harm1113_exact),'--','Color',left_ycolor);
% % plot(ax2_1113_chirp0,time,abs(harm1113_chirp0_2g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax2_1113_chirp0,time,harm1113_exp_omega,'--','Color',right_ycolor)
% % plot(ax2_1113_chirp0,time,harm1113_chirp0_2g_omega,'-','Color',right_ycolor)
% xlim(x_lim_1113); ylim(right_ylim_1113);
% %%
% hold off; ax3_1113_chirp0 = nexttile; hold on;
% title("4 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_1113); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax3_1113_chirp0,time,abs(harm1113_exact),'--','Color',left_ycolor);
% % plot(ax3_1113_chirp0,time,abs(harm1113_chirp0_4g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax3_1113_chirp0,time,harm1113_exp_omega,'--','Color',right_ycolor)
% % plot(ax3_1113_chirp0,time,harm1113_chirp0_4g_omega,'-','Color',right_ycolor)
% xlim(x_lim_1113); ylim(right_ylim_1113);
% %%
% hold off; ax4_1113_chirp0 = nexttile; hold on;
% title("6 Gaussians",'Interpreter','latex');
% yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
% ylim(left_ylim_1113); xlabel("Time (a.u.)",'Interpreter','latex');
% plot(ax4_1113_chirp0,time,abs(harm1113_exact),'--','Color',left_ycolor);
% % plot(ax4_1113,time,abs(harm1113_chirp0_4g_vals),'-','Color',left_ycolor);
% yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
% plot(ax4_1113_chirp0,time,harm1113_exp_omega,'--','Color',right_ycolor)
% % plot(ax4_1113,time,harm1113_chirp0_6g_omega,'-','Color',right_ycolor)
% xlim(x_lim_1113); ylim(right_ylim_1113);
% %%
% linkaxes([ax1_1113_chirp0,ax2_1113_chirp0],'y')
% linkaxes([ax3_1113_chirp0,ax4_1113_chirp0],'y')
% linkaxes([ax1_1113_chirp0,ax3_1113_chirp0],'x')
% linkaxes([ax2_1113_chirp0,ax4_1113_chirp0],'x')