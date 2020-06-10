% clear 
% clc
% close all

% ------------ Config 03
load('saved/config03/opt_FOM_150x150down_50x50x50_t0d80_G10d00_6')
% ------------ Config 12
% load('saved/config12/t0d25/config12_opt_FOM_150x150down_50x50x50_t0d80_G10d00_6')
% load('saved/config12/t0d30/config12_opt_FOM_150x150down_50x50x50_t0d80_G10d00_11')
% ------------ Config 15
% load('saved/config15/opt_FOM_150x150down_50x50x50_t0d80_G05d00_7')% seems the better one
% ------------ Triple Point Shock
% load('saved/trp/opt_FOM_400x200down_40x60x50_t0d25_G10d00_2')
%%
addpathhere
addpath_host([])
q_grid = q;
interptype = 'linear'; % note that spline causes oscilation 
coarse_to_fine_grid_interp = 'linear';
if or( my_case == 12, my_case == 3)
    if_sym_grid = 1;
else
    if_sym_grid = 1;
end
myinterp = myinterp_fun('interp2'); % interp2 or griddata
%% ---------------------- Up-sampling the grid
fprintf('Up-sampling the grid \n')
size_x = size(M_downsampled.rho, 1);
size_y = size(M_downsampled.rho, 2);
size_t = size(M_downsampled.rho, 3);
[x_morphing_fine , y_morphing_fine, ...
 x_f             , y_f            ] = q_coarse_2_grid_fine_XY(size_x, size_y, size_t, q_grid, ...
                                                                 x_coarse, y_coarse, t_coarse, x_fine, y_fine, t_fine, ...
                                                                 k_d, Flag, if_sym_grid, coarse_to_fine_grid_interp);
%% Grid
fig_grid = figure('units','normalized','outerposition',[0 0 1 1]);
size_t_plot = 12; prediction_ratio=1;
Ax = zeros(1,size_t_plot);
tcounter = 1;
for tcount = floor(linspace(1,prediction_ratio*Nt,size_t_plot))
    Ax(tcounter) = subplot(3,4,tcounter);hold all
    nn=20;
    surface('EdgeColor','red','FaceColor','none',...
        'CData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) ),...
        'ZData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) ),...
        'YData',y_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount),...
        'XData',x_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount)); axis equal
%     nn=60;
%     grey_factor=0.75;
%     surface('EdgeColor',grey_factor*[1,1,1],'FaceColor','none',...'LineStyle','--',...
%         'CData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) )-1,...
%         'ZData', ones(size(floor(linspace(1,Ny,nn)),2), size(floor(linspace(1,Nx,nn)),2) )-1,...
%         'YData',y_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount),...
%         'XData',x_morphing_fine(floor(linspace(1,Ny,nn)),floor(linspace(1,Nx,nn)),tcount)); axis equal
    axis off
    title( ['t: ',num2str(tcount),'/',num2str(Nt)] )
    xlabel('x')
    ylabel('y')
    tcounter=tcounter+1;
end