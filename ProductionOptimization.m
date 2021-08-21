% Geothermics Manuscript
% Production Optimization
% Jeremy Patterson
% July 2018

%% Clean Environment
close all; clear; clc

%% Specify Directory

%School Computer
% addpath (genpath('/Users/jpatterson7/Documents/Google Drive/Patterson_Workspace/Supporting.Functions'))

%Mike's computer
% addpath(genpath('/Users/cardiff/Google Drive/WorkSync/SharedSpaces/Patterson_Workspace/Supporting.Functions'))

%Home Computer
addpath (genpath('/Users/jeremypatterson/Documents/Google Drive/Patterson_Workspace/Supporting.Functions'))

%% Physical Constants

g = 9.81;       % Acceleration due to gravity [m/s^2]

% Water Properties
rho_w = 983;   % Water Density [kg/m^3]
mu = 4.7e-4;    % Dynamic viscosity [kg/(m s)]
Cp_w = 4000;    % Water Specific Heat [J/(kg K)]
 
% Hydraulic Properties
T = 1.3e-2;     % Transmissivity [m^2/s] 

% Reservoir Geometry
num_frac = 1 : 10; % Number of Fractures

r_bore = 0.25;  % Borehole Radius [m]
L = 1500;       % Distance Between Wells [m]
Thick = 1000;   % Reservoir Thickness [m]

% Operating Conditions
Q = log(0.2); % Total Volumetric Flow Rate [m^3/s]

% Thermal Properties
lambda = 3;     % Reservoir Rock Thermal Conductivity [W/(m k)]
rho_r = 2500;   % Reservoir Rock Density [kg/m^3]
Cp_r = 1000;    % Reservoir Rock Specific Heat [J/(kg K)]

rhoCp_r = rho_r * Cp_r; % Volumetric Heat Capacity [J/(m^3 K)]

TempInj = 80;   % Injection Water Temperature [deg C]
TempInit = 190; % Initial Production Water Temperature [deg C]

% Time Discretization
num_years = 30;
num_steps = 100;

% Cost Variables

AnnDisc = 0.07;     % Annual Discount Rate [Mines and Nathwani, 2013]
Price = 0.10 * 1e3; % Energy Price [$ / MWh] 

%% Frictional Head Losses

r = 0.33; % Pipe radius [m]
epsilon = 0.045 * 1e-3; % Roughness Coefficient [m]
Mult = 1.5;             % Pipe Tortuosity Multiplier

eff = [0.1; 0.7];       % [Power Plant Efficiency; Pump Efficiency]
InstCap = 26;           % Installed Plant Capacity [MW]

%% Present Value Optimization (Base Case)

tic
PresVal = @(Q) CostOptim(g, mu, epsilon, Thick, L, r, r_bore, T, lambda, rhoCp_r, rho_w, Cp_w, TempInj,...
    TempInit, exp(Q), num_frac, eff(1), eff(2), num_years, InstCap, Price, AnnDisc, Mult);

Q0 = Q;
options = optimset('Display' , 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000);

[Q_opt] = fminsearch(PresVal, Q0, options);
toc

%% Forward Model w/ Optimal Production Rate 

[ENPV, NPV, FracTemp, Energy, ThermPowerRate, ParasiticPower, Revenue, OM, Profit, Time] = CostOptim(g, mu, epsilon, Thick, L, r, r_bore, T, lambda, rhoCp_r, rho_w, Cp_w, TempInj,...
    TempInit, exp(Q_opt), num_frac, eff(1), eff(2), num_years, InstCap, Price, AnnDisc, Mult);

%% T Sensitivity (Figure 8) 

% T = [T*1e-1 T T*1e1];
% for i = 1 : numel(T)
%     [ENPV(i), NPV(:,i), FracTemp(:,i), Energy(:,i), ThermalPowerRate(:,i), ParasiticPower(:,i), Revenue(:,i), OM(:,i), Profit(:,i), Time] = CostOptim(g, mu, epsilon, Thick, L, r, r_bore, T(i), lambda, rhoCp_r, rho_w, Cp_w, TempInj,...
%         TempInit, exp(Q_opt), 2, eff(1), eff(2), num_years, InstCap, Price, AnnDisc, Mult);
% end
%% Optimal Q L vs T Parameter Space (Figure 7)

% Tvec = linspace(log(1.3e-3), log(1.3e-1), 50);
% Lvec = linspace(500, 2500, 50);
% [Tv, Lv] = meshgrid(Tvec, Lvec);
% 
% 
% for i = 1 : numel(Tvec)
%     for j = 1 : numel(Lvec)
%         
%         PresVal = @(Q) CostOptim(g, mu, epsilon, Thick, Lv(i,j), r, r_bore, exp(Tv(i,j)), lambda, rhoCp_r, rho_w, Cp_w, TempInj,...
%             TempInit, exp(Q), num_frac, eff(1), eff(2), num_years, InstCap, Price, AnnDisc, Mult);
%         
%         [Q_opt(i,j)] = fminsearch(PresVal, Q0, options);
%         
%         [ENPV(i,j)] = CostOptim(g, mu, epsilon, Thick, Lv(i,j), r, r_bore, exp(Tv(i,j)), lambda, rhoCp_r, rho_w, Cp_w, TempInj,...
%             TempInit, exp(Q_opt(i,j)), num_frac, eff(1), eff(2), num_years, InstCap, Price, AnnDisc, Mult);
%     end
% end
% ENPV = -ENPV;
% 
% save TLQmesh.mat Tv Lv Q_opt ENPV

%% Figures

% School
print_dir = '/Users/jpatterson7/Documents/Google Drive/Patterson_Workspace/Drafts/Geothermics 2019/Figures/';
% Home
% print_dir = '/Users/jeremypatterson/Documents/Google Drive/Patterson_Workspace/Drafts/Geothermics 2019/Figures/';
 
% figure(2)
% clf
% ax = gca;
% plot(Time, ThermPowerRate(:,2), 'r', 'LineWidth', 2);
% hold on
% plot([Time(1) Time(end)], [ParasiticPower(1) ParasiticPower(1)], 'k', 'LineWidth', 2);
% plot([Time(1) Time(end)], [InstCap InstCap], 'b', 'LineWidth', 2);
% xlabel('Time (years)')
% ylabel('Power Production Rate (MW/hr)')
% ylim([0 60])
% leg = legend('Produced Power', 'Parasitic Power', 'Plant Capacity');
% ax.FontSize = 16;
% leg.FontSize = 14;
% set(2, 'Position', [0 0 938 750])
% file = 'fig2';
% print_file = [print_dir file];
% print(print_file, '-djpeg', '-r1200')

% figure(3)
% clf
% ax = gca;
% plot(Time, Revenue(:,2) * 1e-6, 'k', 'LineWidth', 2);
% hold on
% plot(Time, OM(:,2) * 1e-6, 'r', 'Linewidth', 2)
% xlabel('Time (years)')
% ylabel('Dollars in millions')
% leg = legend('Annual Revenue', 'Annual O&M Costs');
% ax.FontSize = 18;
% leg.FontSize = 18;
% set(3, 'Position', [0 0 938 750])
% file = 'fig3';
% print_file = [print_dir file];
% print(print_file, '-djpeg', '-r1200')

% figure(5)
% clf
% ax = gca;
% h1 = plot(Time, FracTemp(:,1), 'LineWidth', 2);
% hold on
% h2 = plot(Time, FracTemp(:,2), 'LineWidth', 2);
% h3 = plot(Time, FracTemp(:,3), 'LineWidth', 2);
% ylim([80 200])
% xlabel('Time (years)')
% ylabel('Temperature (\circC)')
% leg = legend([h1 h2 h3],...
%     {[num2str(num_frac(1)),' fractures'], ...
%      [num2str(num_frac(2)),' fractures'],...
%      [num2str(num_frac(3)),' fractures']}, 'Location', 'NorthEast');
% ax.FontSize = 18;
% leg.FontSize = 14;
% grid on
% set(5, 'Position', [0 0 938 750])
% file = 'fig5';
% print_file = [print_dir file];
% print(print_file, '-djpeg', '-r1200')

% figure(6)
% clf
% subplot(1,2,1)
% ax = gca;
% h1 = plot(Time, Profit(:,1), 'LineWidth', 2);
% hold on
% h2 = plot(Time, Profit(:,2), 'LineWidth', 2);
% h3 = plot(Time, Profit(:,3), 'LineWidth', 2);
% plot([Time(299) Time(299) Time(299)], NPV(1:3), 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.5)
% ylim([0 200])
% xlabel('Time (years)')
% ylabel('Profit ($ in millions)')
% leg = legend([h1 h2 h3],...
%     {[num2str(num_frac(1)),' fractures'], ...
%      [num2str(num_frac(2)),' fractures'],...
%      [num2str(num_frac(3)),' fractures']}, 'Location', 'NorthEast');
% ax.FontSize = 18;
% leg.FontSize = 14;
% grid on

% subplot(1,2,2)
% ax = gca;
% plot(num_frac, NPV, 'k--o', 'LineWidth', 1.5, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% ylim([100 200])
% xlabel('Number of Fractures')
% ylabel('Net Present Value ($ in millions)')
% ax.FontSize = 18;
% grid on
% set(6, 'Position', [0 0 2244 750])
% file = 'fig6';
% print_file = [print_dir file];
% print(print_file, '-djpeg', '-r1200')

load TLQmesh.mat
figure(7)
clf
% subplot(1,2,1)
ax = gca;
surface(Tv, Lv, (Q_est./2))
hold on
surface(Tv, Lv, ENPV)
colormap(jet)
xlabel('ln (T)')
ylabel('Well Spacing (m)')
zlabel('Q_{opt} (m^3/s)')
grid on
ax.FontSize = 18;
view([-45 35])

% subplot(1,2,2)
% ax = gca;
% surface(Tv, Lv, ENPV)
% colormap(jet)
% xlabel('ln (T)')
% ylabel('Well Spacing (m)')
% zlabel('ENPV ($ in millions)')
% zlim([60 200])
% grid on
% ax.FontSize = 18;
% view([-45 35])
% set(7, 'Position', [0 0 2244 750])
% file = 'fig7';
% print_file = [print_dir file];
% print(print_file, '-djpeg', '-r1200')

% figure(8)
% clf
% subplot(1,2,1)
% ax = gca;
% h1 = plot(Time, FracTemp(:,1), 'or', 'MarkerFaceColor', 'r', 'LineWidth', 2);
% hold on
% h2 = plot(Time, FracTemp(:,2), '^b', 'MarkerFaceColor', 'b', 'LineWidth', 2);
% h3 = plot(Time, FracTemp(:,3), 'dk', 'MarkerFaceColor', 'k', 'LineWidth', 2);
% ylim([100 200])
% xlabel('Time (years)')
% ylabel('Water Temperature (\circC)')
% leg = legend([h1 h2 h3],...
%     {['T = ', num2str(T(1)),' m^2/s'], ...
%      ['T = ', num2str(T(2)),' m^2/s'],...
%      ['T = ', num2str(T(3)),' m^2/s']}, 'Location', 'NorthEast');
% ax.FontSize = 18;
% leg.FontSize = 14;
% 
% subplot(1,2,2)
% ax = gca;
% h1 = plot(Time, Energy(:,1), 'LineWidth', 2);
% hold on
% h2 = plot(Time, Energy(:,2), 'LineWidth', 2);
% h3 = plot(Time, Energy(:,3), 'LineWidth', 2);
% ylim([5 30])
% xlabel('Time (years)')
% ylabel('Power Production Rate (MW/hr)')
% leg = legend([h1 h2 h3],...
%     {['T = ', num2str(T(1)),' m^2/s'], ...
%      ['T = ', num2str(T(2)),' m^2/s'],...
%      ['T = ', num2str(T(3)),' m^2/s']}, 'Location', 'NorthEast');
% ax.FontSize = 18;
% leg.FontSize = 14;
% set(8, 'Position', [0 0 2244 750])
% file = 'fig8';
% print_file = [print_dir file];
% print(print_file, '-djpeg', '-r1200')

% figure
% clf
% ax = gca;
% plot(Time, Profit, 'LineWidth', 1.5)
% xlabel('Time (years)')
% ylabel('Profits ($ in Millions)')
% ax.FontSize = 12;