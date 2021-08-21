% Bodvarson and Tsang (1982) 
% Thermal Front Movement Analytical Solution
% Code author: Jeremy Patterson
% Last updated: Aug 2019

% This function calculates thermal breakthrough at a distance from an
% injection well using the analytical solution developed by Bodvarsson and
% Tsang (1982).
% 
% OUTPUTS 
% FracTemp [numt x numf] - Matrix of water temperature at a specified point
% at each time step. Each column represents a realization with a specified
% number of fractures. Each row represents an individual time step
%
% INPUTS 
% L - Distance from injection well [m]
% D - Total reservoir thickness [m]
% q - Volumetric flow rate per fracture [m^3/s]
% b - Fracture aperture [m]
% rho_w - Water density [kg/m^3]
% Cp_w - Water specific heat capacity [J/(kg K)]
% lambda - Reservoir thermal conductivity [W/(m k)]
% rhoCp_r - Reservoir volumetric heat capacity []
% temp_init - Initial water temperature [deg C] 
% temp_inj - Injection water temperature [deg C]
% t [numt x 1] - Time vector [s]

function [FracTemp] = BodvTemp(L, D, q, b, rho_w, Cp_w, lambda, rhoCp_r, temp_init, temp_inj, t)
%% Dimensionless Parameters
numt = numel(t);
% Dimensionless Parameters
theta = (rho_w .* Cp_w .* b) ./ (rhoCp_r .* D);
eta = 0; % Dimensionless Vertical Distance [-]
xi = (lambda * pi * L^2 .* (2+theta)) ./ (rho_w * Cp_w .* q .* D); 

tau = zeros(numt, numel(theta));
for i = 1 : numel(theta)
tau(:,i) = (lambda .* t) ./ (rhoCp_r * D(i)^2); % Dimensionless Time [-]
end
clear i

numf = numel(theta);
num_eta = numel(eta);

%% Simulated Water Temperatures

% Inverse Laplace Parameters
err = 1e-6;
Nk = 40;
tmax_acc = 2 * max(tau);


for k = 1 : num_eta
    for j = 1 : numf
        if k == 1
            % Fracture Temperature
            Temp = @(p) (1/p) * exp(-(((theta(j)*p) + 2*sqrt(p) * tanh(sqrt(p))) .* xi(j))/ 2 + theta(j));
            for i = 1 : numt
                FracTempNonD(i,j,k) = invlaplace_dehoog2(Temp, tau(i,j), tmax_acc(j), err, Nk);
                FracTemp(i,j) = (FracTempNonD(i,j,k) .* (temp_inj - temp_init)) + temp_init; 
            end
        else
%             % Matrix Temperature
%             Temp = @(p) 1/p*exp(-((theta(j)*p + 2*p^.5*tanh(p^.5))*xi(j)/(2+theta(j))))...
%                 *(cosh(p^.5*eta(k)) - sinh(p^.5*eta(k))*tanh(p^.5));
%             
%             for i = 1 : numt
%                 MatrixTempNonD(i,j,k-1) = invlaplace_dehoog2(Temp, tau(i,j), tmax_acc(j), err, Nk);
%                 MatrixTemp(i,j,k-1) = (MatrixTempNonD(i,j,k-1) .* (temp_inj - temp_init)) + temp_init; 
%             end
        end
    end
end


