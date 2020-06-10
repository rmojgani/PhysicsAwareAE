clear
clc
close all
% ------------ Config 03
load('saved/config03/opt_FOM_150x150down_50x50x50_t0d80_G10d00_6')
% =========================================================================
lambda_rpca = 10;
q_grid = q;
clear q
kM = 2:2:20;
tmax = 0.8;
Nt 	 = 2000*tmax;
% t_fine = linspace(0, tmax, Nt)';
% % % dx = x_fine(2)-x_fine(1);
% dt = t_fine(2)-t_fine(1);
interptype = 'linear'; % note that spline causes oscilation 
coarse_to_fine_grid_interp = 'linear';
if_sub_IC = 0; fprintf('--> Not subtracting inital condition \n')
% if_sub_IC = 1; fprintf('--> With subtracting inital condition \n')
t_save_period = 1;
if_sym_grid = 1;
myinterp = myinterp_fun('interp2'); % interp2 or griddata
%%
main_solve.m