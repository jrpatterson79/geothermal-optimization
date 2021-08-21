% Parasitic Power Loss
% Code Author: Jeremy Patterson
% Last Updated: August 2018
%
% This code determines the parasitic power requirements to pump water from
% a production well to an injection well based on the head gradient between
% the two wells. The head gradient is determined using an analytical
% solution presented in Haitjema (1995). 
%
% OUTPUTS:
% ParasiticPower - Power consumption rate [MW]
%
% INPUTS:
%
% g - Acceleration due to gravity [m/s^2]
% mu - Dynamic viscosity [kg/(m s)]
% epsilon - Pipe roughness coefficient [m]
% L - Distance between wells [m]
% r - Pipe radius [m] 
% r_bore - Borehole radius [m]
% T - Total fracture transmissivity [m^2/s]
% b - Fracture aperture [m]
% rho_w - Water density [kg / m^3]
% Q - Total volumetric flow rate [m^3/s]
% eff_pump - Pump efficiency [-]
% Mult - Tortuosity multiplier to account for extra pipe length greater
%        than the well pair distance [-]

function ParasiticPower = PowerLoss(g, mu, epsilon, L, r, r_bore, T, b, rho_w, Q, eff_pump, Mult)

%% Frictional Head Losses
PipeLength = Mult * L;
Area = pi * r^2; % Pipe cross-sectional area [m^2]
Vel = Q / Area; % Fluid velocity [m/s^2]

Re = (rho_w * Vel * r * 2) / mu; % Reynolds Number [-]

f0 = 2e-2; % Friction Factor Initial Guess

f = Colebrook(Re, epsilon, r*2, f0); % Friction Factor [-]
FricLoss = f * (PipeLength / (r * 2)) * (Vel^2 / (2 * g)); % Friction Head Loss [m]

%% Head Difference Between Wells
Wells = [-(L/2) (L/2)]; %[Injection Well (m) Production Well (m)]
Loc = [Wells(1)+r_bore Wells(2)-r_bore]; % [Injection Borehole Wall (m) Production Borehole Wall (m)]

Phi = @(x, Well1, Well2) (Q / (4 * pi)) * log((x - Well1)^2/(x - Well2)^2);

% Discharge Potential at Injection Borehole Wall
InjPhi = Phi(Loc(2), Wells(1), Wells(2));

% Discharge Potential at Production Borehole Wall
ProdPhi = Phi(Loc(1), Wells(1), Wells(2));

% Hydraulic Head at Injection and Production wells [m]
InjHead = (InjPhi + ((1/2) .* T .* b)) ./ T;
ProdHead = (ProdPhi + ((1/2) .* T .* b)) ./ T;

% Head Difference Between Wells [m]
dh = (InjHead - ProdHead) + FricLoss;

% Power Required to Pump Water [MW]
ParasiticPower = ((rho_w * dh * Q * g) ./ eff_pump) * 1e-6;
