% Colebrook Equation
% Code Author: Jeremy Patterson
% Last Update: Aug 2018

% This code determines the empirical f parameter in the Colebrook equation.
% 
% OUTPUTS:
% f_est = Empirical f parameter in Colebrook equation
%
% INPUTS: 
% Re - Reynolds Number [-]
% epsilon - Pipe roughness factor [m]
% D - Pipe diameter [m]
% f0 - Initial guess at f parameter

function [f_est] = Colebrook(Re, epsilon, D, f0)

left = @(f) 1 / sqrt(f);
right = @(f) -2 * log10(((epsilon/D)/3.7) + (2.51/(Re*sqrt(f))));

obj_fxn = @(f) norm(left(f) - right(f))^2;

f_est = fminsearch(obj_fxn, f0);
