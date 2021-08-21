% Net Present Value Optimization
% Code Author: Jeremy Patterson
% Last Update: August 2019
%
% This code takes given reservoir parameters and determines the production
% rate which optimizes the expected net present value for a given number
% realizations. The base case scenario in the paper is based on 10
% realizations with the number of fractures in the reservoir incresasing
% incrementally from 1 to 10. 
% 
% 
% All matrices [numt x numf] in size follow the convention that each row
% represents an individual monthly time step, and each column represents an
% individual realization with a given number of fractures.
%
% OUTPUTS:
% ENPV [numf x 1] - Expected net present value for each realization [$ in
%                   millions]
% NPV - [numt x numf] - Monthly net present value [$ in millions]
% FracTemp [numt x numf] - Water temperature in fracture center using
%                          Bodvarsson and Tsang (1982) solution. [deg C] 
% PowerRate - [numt x numf] - Electrical power production rate [MW].  
% ThermalPowerRate - [numt x numf] - Thermal power production rate [MW].
% Parasitic Power - [numt x numf] - Parasitic power consumption rate [MW].
% Revenue - [numt x numf] - Monthly revenue [$ in millions]
% OM - [numt x numf] - Montlhly operational costs [$ in millions]
% Profit - [numt x numf] - Monthly profits
% Time - [numt x 1] - Individual time steps [s]
% 
% INPUTS:
% g - Acceleration due to gravity [m/s^2]
% mu - Dynamic viscosity [kg/(m s)]
% epsilon - Roughness coefficient [m]
% Thick - Reservoir thickness [m]
% L - Distance between wells [m]
% r - Pipe radius [m]
% r_bore - Borehole radius [m]
% T - Summative fracture transmissivity [m^2/s]
% lambda - Reservoir thermal conductivity [W/(m k)]
% rhoCp_r - Reservoir volumetric heat capacity [J/(m^3 K)]
% rho_w - Water density [kg/m^3]
% Cp_w - Water specific heat capacity [J/(kg K)]
% temp_inj - Injection water temperature [deg C]
% temp_init - Initial water temperature [deg C]
% Q - Total volumetric flow rate [m^3/s]
% num_frac [numf x 1] - Vector of number of fractures for each realization
% eff_plant - Power plant efficiency [-]
% eff_pump - Pump efficiency [-]
% num_years - Number of years to simulate
% InstCap - Installed power plant production capacity [MW]
% Price - Electricity selling price [$ / MW]
% AnnDisc - Annual discount rate used for NPV analysis [-]
% Mult - Pipe tortuosity multiplier

function varargout = CostOptim(g, mu, epsilon, Thick, L, r, r_bore, T, lambda, rhoCp_r, rho_w, Cp_w, temp_inj,...
   temp_init, Q, num_frac, eff_plant, eff_pump, num_years, InstCap, Price, AnnDisc, Mult)
%% Thermal Breakthrough Parameters
% Time Discretization
dt = 86400 * 30;
max_t = 86400 * 365 * num_years;
t = (0 : dt : max_t)';
numt = numel(t);

Time = (t ./ (86400 * 365));

% Fracture Properties
D = Thick ./ num_frac;                                   % Fracture Half Spacing [m]
b = ((12 .* T .* mu) ./ (rho_w .* g.* num_frac)).^(1/3); % Fracture Aperture [m]
q = Q ./ num_frac;                                       % Volumetric flow rate per fracture [m^3/s]

numf = numel(num_frac);

%% Physical Reservoir Model

% Fluid Temperature Model
[FracTemp] = BodvTemp(L, D, q, b, rho_w, Cp_w, lambda, rhoCp_r, temp_init, temp_inj, t);
Delta_T = FracTemp - temp_inj;  % Temperature Difference [deg C]

% Power Production Model
ThermalPower = Q * rho_w * Cp_w .* Delta_T * 1e-6;  % Thermal Power [MW]
ParasiticPower = PowerLoss(g, mu, epsilon, L, r, r_bore, T, b, rho_w, Q, eff_pump, Mult); % Parasitic Power Loss [MW]

NetPower = ThermalPower - ParasiticPower;

% Quantifies Revenue Loss to Overproduction
% ProducedEnergy = NetPower * eff_plant * (dt / 3600);

ProducedEnergy = zeros(size(NetPower));
for i = 1 : numel(NetPower)
    if NetPower(i) * eff_plant > InstCap
        ProducedEnergy(i) = InstCap * (dt / 3600);
    else
        ProducedEnergy(i) = NetPower(i) * eff_plant * (dt / 3600);
    end
end

ThermalPowerRate = ThermalPower * eff_plant;
%% Economic Cost Model
% Monthly Operation and Maintenance Costs
OM = 20 * exp(-0.0025 *(InstCap-5)) * ProducedEnergy;

Revenue = zeros(numt,numf);
for j = 1 : numf
    for i = 1 : numt
        Revenue(i,j) = (ProducedEnergy(i,j) * Price) / (1+AnnDisc)^Time(i);
    end
end

Profit = cumsum(Revenue - OM) * 1e-6;
NPV = max(Profit);
ENPV = -mean(NPV);
PowerRate = ProducedEnergy / (dt/3600);

varargout = {ENPV, NPV, FracTemp, PowerRate, ThermalPowerRate, ParasiticPower, Revenue, OM, Profit, Time};
