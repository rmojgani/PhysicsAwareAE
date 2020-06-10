%% Finer Grid
multiplier  = 1;
t_save_period = t_save_period * multiplier;
% Nx   = 150 * multiplier;
% Ny 	 = 150 * multiplier;
Nt 	 = 2000*tmax * multiplier;
% x_fine = linspace(0, xmax, Nx)';
% y_fine = linspace(0, ymax, Ny)';
% t_fine = linspace(0, tmax, Nt)';
% dx = x_fine(2)-x_fine(1);
dt = t_fine(2)-t_fine(1);
% clear m
%%
t_save_range = [1:t_save_period:Nt-1, Nt];
t_save_size = length(t_save_range);
%% ---------------------- Eulerian - solve matlab_cfd_shock_tube_TVDRK2\main.m
fprintf('Eulerian - Implicit \n')
clear t
main
save('savenow','-v7.3')
%%
fields = fieldnames(M_XY);
for varcount = 1:numel(fields)
    M.(fields{varcount}) = reshape( M_XY.(fields{varcount}), [Nx*Ny, t_save_size]);
end
% %% ---------------------- Up-sampling the grid
fprintf('Up-sampling the grid \n')
size_x = size(M_downsampled.rho, 1);
size_y = size(M_downsampled.rho, 2);
size_t = size(M_downsampled.rho, 3);
[x_morphing_fine , y_morphing_fine, ...
 x_f             , y_f            ] = q_coarse_2_grid_fine_XY(size_x, size_y, size_t, q_grid, ...
                                                                 x_coarse, y_coarse, t_coarse, x_fine, y_fine, t_fine, ...
                                                                 k_d, Flag, if_sym_grid, coarse_to_fine_grid_interp);
% %% ---------------------- ALE FOM [on predefined grid movement] - Interpolation [To be changed]
fprintf('ALE FOM - With interpolation (Calc. ALE ) \n')
fields = fieldnames(M_XY);
for varcount = 1:numel(fields)
    t_save = 0;
    for tcount = t_save_range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t_save = t_save + 1;
         M_ALE_XY_T = ...
            myinterp( x_f, y_f, ...
            M_XY.(fields{varcount})(:,:,t_save)', ...
            x_morphing_fine(:,:,t_save_range(t_save)), y_morphing_fine(:,:,t_save_range(t_save)) );
        M_ALE_XY.(fields{varcount})(:,:,t_save) = M_ALE_XY_T';
    end
end
clear M_ALE_XY_T M_XY
% ---------------------- Interpolate ALE FOM on Eulerian grid
% fprintf('ALE FOM - With interpolation (P(ALE) on Eulerian ) \n')
% for varcount = 1:numel(fields)
%     t_save = 0;
%     for tcount = t_save_range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         t_save = t_save + 1;
%         M_ALE_on_stationary_grid_XY.(fields{varcount})(:,:,t_save) = ...
%             myinterp( x_morphing_fine(:,:,t_save_range(t_save)), y_morphing_fine(:,:,t_save_range(t_save)),...
%             M_ALE_XY.(fields{varcount})(:,:,t_save),...
%             x_f, y_f);%, 'linear', 0);
%     end
% end
% clear M_ALE_on_stationary_grid_XY
% %% ----------------------ALE Basis
fprintf('Basis on morphing grid \n')
for varcount = 1:numel(fields)
       M_ALE = reshape( M_ALE_XY.(fields{varcount}), [Nx*Ny, t_save_size] );%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [UU.(fields{varcount}), SS.(fields{varcount}) ,VV.(fields{varcount}) ] = svd(M_ALE, 'econ');
       
end
clear M_ALE_XY
%%
kMS=[1,kM];
for k = kM
    kcount = find(kM == k);
    fprintf('Reducing size of basis functions \n')
    for varcount = 1:numel(fields)
        [U.(fields{varcount}), V.(fields{varcount})] = ...
            USV_reduce(UU.(fields{varcount}), SS.(fields{varcount}), VV.(fields{varcount}), k);
        U_XY.(fields{varcount}) = reshape( U.(fields{varcount}), [Nx, Ny, k]);
    end
    %% ---------------------- ALE Proj + ROM [with predefined grid movement] - Basis on Eulerian grid
    tic
    fprintf(['ALE ', num2str(k),'-dim Proj and ROM \n'])
    for varcount = 1:numel(fields)
         a_0.(fields{varcount}) =   V.(fields{varcount})(:,1);
    end
    [M_a_ALE, M_ALE_ROM_on_stationary_grid, M_ALE_Proj_on_stationary_grid] = ...
        Euler_2d_Rie_ROM_basis_on_Eulerian_U(a_0, Nx, Ny, Nt, dt, gamma, U_XY, V, t_save_period, ...
                                                x_morphing_fine, y_morphing_fine, x_f, y_f, myinterp);
    toc
    %% ---------------------- Eulerian Projection
    fprintf(['Eulerian ', num2str(k),'-dim Projection \n'])
    for varcount = 1:numel(fields)      
        M_Eulerian.(fields{varcount}) = reshape( M.(fields{varcount}), [Nx*Ny, t_save_size] );%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [UU_Eul, SS_Eul, VV_Eul] = svd(bsxfun(@minus, M_Eulerian.(fields{varcount}), M_Eulerian.(fields{varcount})(:,1)*if_sub_IC), 'econ');
        
        [U_Eul.(fields{varcount}), V_Eul.(fields{varcount}) ] = USV_reduce(UU_Eul, SS_Eul, VV_Eul, k);
        M_Eulerian_Proj.(fields{varcount}) = bsxfun(@plus, U_Eul.(fields{varcount}) * V_Eul.(fields{varcount}), M_Eulerian.(fields{varcount})(:,1)*if_sub_IC );
    end
    %% ---------------------- Eulerian ROM
    fprintf(['Eulerian ', num2str(k),'-dim ROM \n'])
    for varcount = 1:numel(fields)
         a_0.(fields{varcount}) =   V_Eul.(fields{varcount})(:,1);
    end
    for varcount = 1:numel(fields)
         M_0.(fields{varcount}) = M_Eulerian.(fields{varcount})(:,1)*if_sub_IC;
    end
    M_a = Euler_2d_Rie_ROM_IC(a_0, Nx, Ny, Nt, dt, gamma, U_Eul, M_0, t_save_period);
    for varcount = 1:numel(fields)
         M_Eulerian_ROM.(fields{varcount}) = bsxfun(@plus, U_Eul.(fields{varcount}) * M_a.(fields{varcount}), M_0.(fields{varcount})  );
    end
    %%
    fields = fieldnames(M_Eulerian);
    for filecount = 1:numel(fields)
        
        error.Eul_Proj.(fields{filecount})(kcount) = ...
            norm(  M_Eulerian_Proj.(fields{filecount})' - M_Eulerian.(fields{filecount})'  , 'fro');
        
        error.Eul_ROM.(fields{filecount})(kcount) = ...
            norm(  M_Eulerian_ROM.(fields{filecount})' - M_Eulerian.(fields{filecount})'  , 'fro');
        
        error.ALE_Proj_PALE_2_Eul.(fields{filecount})(kcount) = ...
            norm(  M_ALE_Proj_on_stationary_grid.(fields{filecount})' - M_Eulerian.(fields{filecount})'  ,'fro');
        
        error.ALE_ROM_PALE_2_Eul.(fields{filecount})(kcount) = ...
            norm(  M_ALE_ROM_on_stationary_grid.(fields{filecount})' - M_Eulerian.(fields{filecount})'  , 'fro');
        
    end
end
%%
display('Saving the results')

clear M_Q
clear UU SS VV

display('Saving the results -- All')
save('main_solve.mat','-v7.3')
display('The End')