close all; clear all;
load("./results/results_0612.mat");

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

data_dir = './data/'; max_intensity = 3e-3;

[FFT,t] = filterHarmonic([data_dir '500000c_0um.mat'],9);
mask = -1000 < t & t < 1000;
time_9_exact = t(mask);
harm9_exact = FFT(mask)*(max_intensity / max(FFT(mask)));
dt_9_exact = abs(time_9_exact(2)-time_9_exact(1));

[FFT,t] = filterHarmonic([data_dir '500000c_0um.mat'],11);
mask = -1000 < t & t < 1000;
time_11_exact = t(mask);
harm11_exact = FFT(mask)*(max_intensity / max(FFT(mask)));
dt_11_exact = abs(time_11_exact(2)-time_11_exact(1));

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
% harm9_chirp0_6g_phase = unwrap(angle(harm9_chirp0_6g_vals));
% harm9_chirp0_6g_omega = zeros(size(harm9_chirp0_6g_phase));

harm9_chirp1_1g_phase = unwrap(angle(conj(harm9_chirp1_1g_vals)));
harm9_chirp1_1g_omega = zeros(size(harm9_chirp1_1g_phase));
harm9_chirp1_2g_phase = unwrap(angle(conj(harm9_chirp1_2g_vals)));
harm9_chirp1_2g_omega = zeros(size(harm9_chirp1_2g_phase));
% harm9_chirp1_4g_phase = unwrap(angle(conj(harm9_chirp1_4g_vals)));
% harm9_chirp1_4g_omega = zeros(size(harm9_chirp1_4g_phase));
% harm9_chirp1_6g_phase = unwrap(angle(harm9_chirp1_6g_vals));
% harm9_chirp1_6g_omega = zeros(size(harm9_chirp1_6g_phase));

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
% harm11_chirp0_6g_phase = unwrap(angle(harm11_chirp0_6g_vals));
% harm11_chirp0_6g_omega = zeros(size(harm11_chirp0_6g_phase));

harm11_chirp1_1g_phase = unwrap(angle(conj(harm11_chirp1_1g_vals)));
harm11_chirp1_1g_omega = zeros(size(harm11_chirp1_1g_phase));
harm11_chirp1_2g_phase = unwrap(angle(conj(harm11_chirp1_2g_vals)));
harm11_chirp1_2g_omega = zeros(size(harm11_chirp1_2g_phase));
% harm11_chirp1_4g_phase = unwrap(angle(conj(harm11_chirp1_4g_vals)));
% harm11_chirp1_4g_omega = zeros(size(harm11_chirp1_4g_phase));
% harm11_chirp1_6g_phase = unwrap(angle(harm11_chirp1_6g_vals));
% harm11_chirp1_6g_omega = zeros(size(harm11_chirp1_6g_phase));

for t = 2:length(time_9_exact)-1
    harm9_exp_omega(t) = (harm9_exp_phase(t+1) - harm9_exp_phase(t-1))/(2 * dt_9_exact);
    harm11_exp_omega(t) = (harm11_exp_phase(t+1) - harm11_exp_phase(t-1))/(2 * dt_11_exact);
end

for t = 2:length(time)-1
    harm9_chirp0_1g_omega(t) = (harm9_chirp0_1g_phase(t+1) - harm9_chirp0_1g_phase(t-1))/(2 * dt);
    harm9_chirp0_2g_omega(t) = (harm9_chirp0_2g_phase(t+1) - harm9_chirp0_2g_phase(t-1))/(2 * dt);
    harm9_chirp0_4g_omega(t) = (harm9_chirp0_4g_phase(t+1) - harm9_chirp0_4g_phase(t-1))/(2 * dt);
%     harm9_chirp0_6g_omega(t) = -(harm9_chirp0_6g_phase(t+1) - harm9_chirp0_6g_phase(t-1))/(2 * dt);

    harm9_chirp1_1g_omega(t) = (harm9_chirp1_1g_phase(t+1) - harm9_chirp1_1g_phase(t-1))/(2 * dt);
    harm9_chirp1_2g_omega(t) = (harm9_chirp1_2g_phase(t+1) - harm9_chirp1_2g_phase(t-1))/(2 * dt);
%     harm9_chirp1_4g_omega(t) = (harm9_chirp1_4g_phase(t+1) - harm9_chirp1_4g_phase(t-1))/(2 * dt);
%     harm9_chirp1_6g_omega(t) = -(harm9_chirp1_6g_phase(t+1) - harm9_chirp1_6g_phase(t-1))/(2 * dt);


    harm11_chirp0_1g_omega(t) = (harm11_chirp0_1g_phase(t+1) - harm11_chirp0_1g_phase(t-1))/(2 * dt);
    harm11_chirp0_2g_omega(t) = (harm11_chirp0_2g_phase(t+1) - harm11_chirp0_2g_phase(t-1))/(2 * dt);
    harm11_chirp0_4g_omega(t) = (harm11_chirp0_4g_phase(t+1) - harm11_chirp0_4g_phase(t-1))/(2 * dt);
%     harm11_chirp0_6g_omega(t) = -(harm11_chirp0_6g_phase(t+1) - harm11_chirp0_6g_phase(t-1))/(2 * dt);

    harm11_chirp1_1g_omega(t) = (harm11_chirp1_1g_phase(t+1) - harm11_chirp1_1g_phase(t-1))/(2 * dt);
    harm11_chirp1_2g_omega(t) = (harm11_chirp1_2g_phase(t+1) - harm11_chirp1_2g_phase(t-1))/(2 * dt);
%     harm11_chirp1_4g_omega(t) = (harm11_chirp1_4g_phase(t+1) - harm11_chirp1_4g_phase(t-1))/(2 * dt);
%     harm11_chirp1_6g_omega(t) = -(harm11_chirp1_6g_phase(t+1) - harm11_chirp1_6g_phase(t-1))/(2 * dt);
end
%%
figure(); tiles_9_chirp0 = tiledlayout(2,2);
title(tiles_9_chirp0,"9th Harmonic Reconstruction - No Chirp",'Interpreter','latex');
tiles_9_chirp0.Padding = 'compact'; tiles_9_chirp0.TileSpacing = 'compact';
%%
ax1_9_chirp0 = nexttile; hold on; shift = -32.5;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1_9_chirp0,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax1_9_chirp0,time - shift,abs(harm9_chirp0_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1_9_chirp0,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax1_9_chirp0,time - shift,harm9_chirp0_1g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax2_9_chirp0 = nexttile; hold on; shift = 67.6;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp0,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax2_9_chirp0,time - shift,abs(harm9_chirp0_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp0,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax2_9_chirp0,time - shift,harm9_chirp0_2g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax3_9_chirp0 = nexttile; hold on; shift = 81;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp0,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax3_9_chirp0,time -shift,abs(harm9_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp0,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax3_9_chirp0,time - shift,harm9_chirp0_4g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax4_9_chirp0 = nexttile; hold on; shift = 0;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp0,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
% plot(ax4_9,time - shift,abs(harm9_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp0,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
% plot(ax4_9,time - shift,harm9_chirp0_6g_omega,'-','Color',right_ycolor)
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
ax1_9_chirp1 = nexttile; hold on; shift = -57.5;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1_9_chirp1,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax1_9_chirp1,time - shift,abs(harm9_chirp1_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1_9_chirp1,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax1_9_chirp1,time - shift,harm9_chirp1_1g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax2_9_chirp1 = nexttile; hold on; shift = 78.7+15.1;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp1,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
plot(ax2_9_chirp1,time - shift,abs(harm9_chirp1_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2_9_chirp1,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax2_9_chirp1,time - shift,harm9_chirp1_2g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax3_9_chirp1 = nexttile; hold on; shift = 0;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp1,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
% plot(ax3_9_chirp1,time -shift,abs(harm9_chirp1_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3_9_chirp1,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
% plot(ax3_9_chirp1,time - shift,harm9_chirp1_4g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
hold off; ax4_9_chirp1 = nexttile; hold on; shift = 0;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_9); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp1,time_9_exact,abs(harm9_exact),'--','Color',left_ycolor);
% plot(ax4_9,time - shift,abs(harm9_chirp1_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4_9_chirp1,time_9_exact,harm9_exp_omega,'--','Color',right_ycolor)
% plot(ax4_9,time - shift,harm9_chirp1_6g_omega,'-','Color',right_ycolor)
xlim(x_lim_9); ylim(right_ylim_9);
%%
linkaxes([ax1_9_chirp1,ax2_9_chirp1],'y')
linkaxes([ax3_9_chirp1,ax4_9_chirp1],'y')
linkaxes([ax1_9_chirp1,ax3_9_chirp1],'x')
linkaxes([ax2_9_chirp1,ax4_9_chirp1],'x')

%%
figure(); tiles_11_chirp0 = tiledlayout(2,2);
title(tiles_11_chirp0,"11th Harmonic Reconstruction - No Chirp",'Interpreter','latex');
tiles_11_chirp0.Padding = 'compact'; tiles_11_chirp0.TileSpacing = 'compact';
%%
ax1_11_chirp0 = nexttile; hold on; shift = -66.8;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1_11_chirp0,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
plot(ax1_11_chirp0,time - shift,abs(harm11_chirp0_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1_11_chirp0,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
plot(ax1_11_chirp0,time - shift,harm11_chirp0_1g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
hold off; ax2_11_chirp0 = nexttile; hold on; shift = 79.6;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2_11_chirp0,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
plot(ax2_11_chirp0,time - shift,abs(harm11_chirp0_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2_11_chirp0,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
plot(ax2_11_chirp0,time - shift,harm11_chirp0_2g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
hold off; ax3_11_chirp0 = nexttile; hold on; shift = -76;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3_11_chirp0,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
plot(ax3_11_chirp0,time -shift,abs(harm11_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3_11_chirp0,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
plot(ax3_11_chirp0,time - shift,harm11_chirp0_4g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
hold off; ax4_11_chirp0 = nexttile; hold on; shift = 0;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4_11_chirp0,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax4_11,time - shift,abs(harm11_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4_11_chirp0,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax4_11,time - shift,harm11_chirp0_6g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
linkaxes([ax1_11_chirp0,ax2_11_chirp0],'y')
linkaxes([ax3_11_chirp0,ax4_11_chirp0],'y')
linkaxes([ax1_11_chirp0,ax3_11_chirp0],'x')
linkaxes([ax2_11_chirp0,ax4_11_chirp0],'x')

%%
figure(); tiles_11_chirp1 = tiledlayout(2,2);
title(tiles_11_chirp1,"11th Harmonic Reconstruction - Chirp",'Interpreter','latex');
tiles_11_chirp1.Padding = 'compact'; tiles_11_chirp1.TileSpacing = 'compact';
%%
ax1_11_chirp1 = nexttile; hold on; shift = -77.3;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1_11_chirp1,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
plot(ax1_11_chirp1,time - shift,abs(harm11_chirp1_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1_11_chirp1,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
plot(ax1_11_chirp1,time - shift,harm11_chirp1_1g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
hold off; ax2_11_chirp1 = nexttile; hold on; shift = 44.8;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2_11_chirp1,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
plot(ax2_11_chirp1,time - shift,abs(harm11_chirp1_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2_11_chirp1,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
plot(ax2_11_chirp1,time - shift,harm11_chirp1_2g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
hold off; ax3_11_chirp1 = nexttile; hold on; shift = 0;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3_11_chirp1,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax3_11_chirp1,time -shift,abs(harm11_chirp1_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3_11_chirp1,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax3_11_chirp1,time - shift,harm11_chirp1_4g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
hold off; ax4_11_chirp1 = nexttile; hold on; shift = 0;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim_11); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4_11_chirp1,time_11_exact,abs(harm11_exact),'--','Color',left_ycolor);
% plot(ax4_11,time - shift,abs(harm11_chirp1_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4_11_chirp1,time_11_exact,harm11_exp_omega,'--','Color',right_ycolor)
% plot(ax4_11,time - shift,harm11_chirp1_6g_omega,'-','Color',right_ycolor)
xlim(x_lim_11); ylim(right_ylim_11);
%%
linkaxes([ax1_11_chirp1,ax2_11_chirp1],'y')
linkaxes([ax3_11_chirp1,ax4_11_chirp1],'y')
linkaxes([ax1_11_chirp1,ax3_11_chirp1],'x')
linkaxes([ax2_11_chirp1,ax4_11_chirp1],'x')
