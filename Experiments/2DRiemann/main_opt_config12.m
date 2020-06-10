clear
close all
clc
addpath(genpath('matlab_cfd_shock_tube_TVDRK2'))
my_case = 12;
k_d = 2;
k_w = 4;
kM = [k_d, k_w];
Gamma_0_x = 0.01/k_d;
Gamma_0_t = 0.01/k_d;
Gamma_0_M = [Gamma_0_x, Gamma_0_t];
%% Code Description
% ******************************************************************************
% solving 1D Euler equations:
%
%
%   Q_{t} - {F(Q)}_{x} = 0
%
%   x \in [0, xmax]
%   t \in [0, tmax]
%
% -----------------------------------------------------------------
%%
multiplier  = 1;
xmax = 1;
ymax = 1;
tmax = 0.25;
Nx   = 150*xmax*multiplier;
Ny 	 = 150*ymax*multiplier;
Nt 	 = 2000*tmax*multiplier;
%%
main_opt