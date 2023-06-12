close all; clear all;
load("./results/results_0612.mat");

left_ylabel = "$$|\tilde{f}(t)|$$ (a.u.)";
right_ylabel = "$$\omega(t)$$ (a.u.)";

left_ycolor = "#0072BD";
right_ycolor = "#A2142F";

left_ylim = [0 3.2e-3];
right_ylim = [0.47 0.55];
x_lim = [-700 900];

dt = abs(time(2)-time(1));
harm9_exp_phase = unwrap(angle(conj(harm9_exp_vals)));
harm9_exp_omega = zeros(size(harm9_exp_phase));
harm9_chirp0_1g_phase = unwrap(angle(conj(harm9_chirp0_1g_vals)));
harm9_chirp0_1g_omega = zeros(size(harm9_chirp0_1g_phase));
harm9_chirp0_2g_phase = unwrap(angle(conj(harm9_chirp0_2g_vals)));
harm9_chirp0_2g_omega = zeros(size(harm9_chirp0_2g_phase));
harm9_chirp0_4g_phase = unwrap(angle(conj(harm9_chirp0_4g_vals)));
harm9_chirp0_4g_omega = zeros(size(harm9_chirp0_4g_phase));
% harm9_chirp0_6g_phase = unwrap(angle(harm9_chirp0_6g_vals));
% harm9_chirp0_6g_omega = zeros(size(harm9_chirp0_6g_phase));
for t = 2:length(time)-1
    harm9_exp_omega(t) = (harm9_exp_phase(t+1) - harm9_exp_phase(t-1))/(2 * dt);
    harm9_chirp0_1g_omega(t) = (harm9_chirp0_1g_phase(t+1) - harm9_chirp0_1g_phase(t-1))/(2 * dt);
    harm9_chirp0_2g_omega(t) = (harm9_chirp0_2g_phase(t+1) - harm9_chirp0_2g_phase(t-1))/(2 * dt);
    harm9_chirp0_4g_omega(t) = (harm9_chirp0_4g_phase(t+1) - harm9_chirp0_4g_phase(t-1))/(2 * dt);
%     harm9_chirp0_6g_omega(t) = -(harm9_chirp0_6g_phase(t+1) - harm9_chirp0_6g_phase(t-1))/(2 * dt);
end

figure(); tiles = tiledlayout(2,2);
title(tiles,"9th Harmonic Reconstruction",'Interpreter','latex');
tiles.Padding = 'compact'; tiles.TileSpacing = 'compact';
%%
ax1 = nexttile; hold on; shift = 0;
title("1 Gaussian",'Interpreter','latex');
yyaxis left; ylabel(left_ylabel,'Interpreter','latex');
ylim(left_ylim); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax1,time,abs(harm9_exp_vals),'--','Color',left_ycolor);
plot(ax1,time - shift,abs(harm9_chirp0_1g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel(right_ylabel,'Interpreter','latex');
plot(ax1,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax1,time - shift,harm9_chirp0_1g_omega,'-','Color',right_ycolor)
xlim(x_lim); ylim(right_ylim);
%%
hold off; ax2 = nexttile; hold on; shift = 93.3;
title("2 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax2,time,abs(harm9_exp_vals),'--','Color',left_ycolor);
plot(ax2,time - shift,abs(harm9_chirp0_2g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax2,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax2,time - shift,harm9_chirp0_2g_omega,'-','Color',right_ycolor)
xlim(x_lim); ylim(right_ylim);
%%
hold off; ax3 = nexttile; hold on; shift = 0;
title("4 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax3,time,abs(harm9_exp_vals),'--','Color',left_ycolor);
plot(ax3,time -shift,abs(harm9_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax3,time,harm9_exp_omega,'--','Color',right_ycolor)
plot(ax3,time - shift,harm9_chirp0_4g_omega,'-','Color',right_ycolor)
xlim(x_lim); ylim(right_ylim);
%%
hold off; ax4 = nexttile; hold on; shift = 0;
title("6 Gaussians",'Interpreter','latex');
yyaxis left; ylabel("$$|\tilde{f}(t)|$$ (a.u.)",'Interpreter','latex');
ylim(left_ylim); xlabel("Time (a.u.)",'Interpreter','latex');
plot(ax4,time,abs(harm9_exp_vals),'--','Color',left_ycolor);
% plot(ax4,time - shift,abs(harm9_chirp0_4g_vals),'-','Color',left_ycolor);
yyaxis right; ylabel("$$\omega(t)$$ (a.u.)",'Interpreter','latex');
plot(ax4,time,harm9_exp_omega,'--','Color',right_ycolor)
% plot(ax2,time - shift,harm9_chirp0_6g_omega,'-','Color',right_ycolor)
xlim(x_lim); ylim(right_ylim);
%%
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
linkaxes([ax1,ax3],'x')
linkaxes([ax2,ax4],'x')
