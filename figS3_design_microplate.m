%% Created and last modified (5/16/22) by Hyungseok Kim (hskimm@mit.edu)
clc; clear; close all;
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultAxesfontWeight', 'normal')
set(0,'DefaultTextFontName', 'Arial')
set(0,'DefaultTextFontSize', 12)
set(0,'DefaultTextfontWeight', 'normal')

linest = {'k-','k--','k-.','k:'};

%% Algae only! No bacteria

% load data
days = [0 2 4 6 8 11 14 17 20];
pt_hex = readmatrix('data/pt_test.txt');
pt_phy = readmatrix('data/pt_flask.txt');

N_t = 100; % time discretiztion by 1/N_t (day)
dt = 86400/N_t; % time increment [s]
D_gel = 3.47e-10; % diffusivity of glycine [m2/s]
t = 1.5e-3; % [m]
Q = 7.2e-6 / 86400; % DOC [umol/s-cell] (Seymour2017)

C_mem = cell(3,4); % stores for each alpha and N
C_mem_phy = cell(3,4);
% alphas = [0.2, 0.5, 1]; % diameter ratio (surr) - (center)
alphas = [0.344]; % used for experiment
% Ns = [3,4,6,9]; % number of surrounding wells
Ns = [6]; % used for experiment
figs = cell(3,2);

for a = 1 % change alpha
    alpha = alphas(a);
    fprintf("alpha=%g\n", alpha);
    D = 0.045 / (1+4*alpha);  % m  (diam. center)
    fig1 = figure();
    fig2 = figure();
    for b = 1 % change N
        N = Ns(b);
        C = zeros(1, 3); % DOC initialize
        C_phy = 0;
        C_cumul = C;
        C_cumul_phy = C_phy;
        for i = 1:N_t*20
            rho = call_culture(pt_hex, i*dt, days);
            rho_phy = call_culture(pt_phy, i*dt, 0:13);
            rho = rho * 1e+3; % cells/l
            rho_phy = rho_phy * 1e+3; % cells/l
            C(1) = rho*Q*dt ...
                + 4*D_gel*dt/(D*t) * C(2) ...
                + (1 - 4*D_gel*dt/(D*t)) * C(1);  % center
            C(2) = 4*D_gel*dt/(D*alpha^2*N*t) * (C(1) - C(2)) ...
                + 4*D_gel*dt/(pi*D*alpha*t) * (C(3) - C(2)) ...
                + C(2);  % inner well
            C(3) = 4*D_gel*dt/(pi*D*alpha*t) * (C(2) - 2*C(3)) ...
                + C(3);  % figureer well
            C_phy = C_phy + Q * rho_phy * dt;
            if rem(i,N_t)==0
                C_cumul = [C_cumul; C]; % uM
                C_cumul_phy = [C_cumul_phy; C_phy]; % uM
            end
        end
        C_mem{a,b} = C_cumul;
        C_mem_phy{a,b} = C_cumul_phy;
        figure(1); plot(0:20, C_cumul(:,2), linest{b}, 'LineWidth', 0.4); hold on;
        plot(0:13, C_cumul_phy(1:14), 'r-', 'LineWidth', 0.5); hold on;
        figure(2); plot(0:20, C_cumul(:,3), linest{b}, 'LineWidth', 0.4); hold on;
        plot(0:13, C_cumul_phy(1:14), 'r-', 'LineWidth', 0.5); hold on;
        writematrix(C_cumul, sprintf('DesignResults/DOC_pm_alph_%g_N_%g.csv', alpha, N));
    end
%     figs{a,1} = fig1;
%     figs{a,2} = fig2;
%     figure(1);
%     xticks([0 5 10 15 20]);
%     set(gca,'TickDir','out');
%     set(gca, 'Units', 'inches', 'Position', [2, 2, 2, 2]);
%     grid on;
%     figure(2);
%     xticks([0 5 10 15 20]);
%     set(gca,'TickDir','out');
%     set(gca, 'Units', 'inches', 'Position', [2, 2, 2, 2]);
%     grid on;
%     saveas(fig1,sprintf('figure/alpha_%g_ring1_conc.eps',alpha),'epsc');
%     saveas(fig2,sprintf('figure/alpha_%g_ring2_conc.eps',alpha),'epsc');
%     close all;
end

% if cond, writematrix(C_mem,'hex-DOC2.csv','Delimiter','tab');
% else, writematrix(C_mem,'hex-nitrate.csv','Delimiter','tab');
% end

%%%%%%%%%%%%%
function [rho] = call_culture(pt0, time, days) % [cells/ml]
% pt0 : reference cell number per volume
% time : timepoint of interest
    time = time / 86400;
    j = 1;
    if time > days(end)
        rho = pt0(end);
    else
        while time > days(j)
            j = j+1;
        end    
        frac = (time-days(j-1))/(days(j)-days(j-1));
        rho = pt0(j-1) * (pt0(j)/pt0(j-1))^frac;
    end
end