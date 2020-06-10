x_fine = linspace(0, xmax, Nx)';
y_fine = linspace(0, ymax, Ny)';
t_fine = linspace(0, tmax, Nt)';
dt = t_fine(2)-t_fine(1);
%%
t_save_period = 100;
t_save_range = [1:t_save_period:Nt-1, Nt];
t_save_size = length(t_save_range);
%%
main
fields = fieldnames(M_XY);
for varcount = 1:numel(fields)
    M.(fields{varcount}) = reshape( M_XY.(fields{varcount}), [Nx*Ny,Nt]);
end
%%
k = 4;
% [pM, uM, M.rho] = Q_2_primitive(M_Q, Nx, gamma);
for varcount = 1:numel(fields)
    % [U,S,V]=svd(M,'econ')
    [U.(fields{varcount}), S.(fields{varcount}), V.(fields{varcount}) ] = ...
        svd(M.rho, 'econ');
    % [U,S,V]=USV_reduce(U,S,V,k)
    [U.(fields{varcount}), V.(fields{varcount}) ] = ...
        USV_reduce(U.(fields{varcount}), S.(fields{varcount}), V.(fields{varcount}), k);
    % norm(M - U*V,'fro')
    norm(M.(fields{varcount}) - U.(fields{varcount})*V.(fields{varcount}), 'fro')
end
clear M S
%% Downsample
xbox_size = 50;
ybox_size = 50;
tbox_size = 50;
for varcount = 1:numel(fields)
%     [M_Q_downsampled.rho, x_coarse, y_coarse, t_coarse, size_x, size_y, size_t] = ...
%     downsampleXY(M_Q, x_fine, y_fine, t_fine, xbox_size, ybox_size, tbox_size);
    [M_downsampled.(fields{varcount}), x_coarse, y_coarse, t_coarse, size_x, size_y, size_t, X_c, Y_c] = ...
        downsampleXY(M_XY.(fields{varcount}), x_fine, y_fine, t_fine, xbox_size, ybox_size, tbox_size);
end
%%
figure()
count = 1;
for varcount = 1:numel(fields)
    subplot(2, numel(fields), varcount);
    surf(X, Y, M_XY.(fields{varcount})(:,:,end) )
    axis equal, shading flat
    view([0 90])
    title((fields{varcount}))
    
    subplot(2, numel(fields), numel(fields)+varcount);
    surf(X_c, Y_c, M_downsampled.(fields{varcount})(:,:,end) )
    axis equal, shading flat
    view([0 90])
    title((fields{varcount}))

end
drawnow
%%
file_name = ['case',num2str(my_case),'_',num2str(Nx),'x',num2str(Ny),'x',num2str(Nt),...
    '_down_',num2str(size_x),'x',num2str(size_y),'x',num2str(size_t),'_t0d',num2str(tmax*100)];
save([file_name,'c'], '-v7.3')
% clear M_XY % Clear
%% Recontruct the initial condition and problem parameters
% Optimization parameters
MaxIter = 200;
SaveIter = 10;
%% Inital condition for basis functions
Dx  = sparse(fd_normal(size_x,3,x_coarse,1));
U_d_x(:,1:k_d) = Dx*reshape(M_downsampled.rho(:,1,1:k_d),[size_x,k_d]);
[U_d_x,~,~]    = svd(U_d_x,'econ'); % orthonormalize U_x
U_d_x = U_d_x(2:1:end-1,:); % first & last row zero

V_d_x = zeros(k_d,size_t);

q_x = [U_d_x(:);V_d_x(:)];
Flag.x = [length(U_d_x(:)), length(V_d_x(:))];
%% Inital condition for basis functions
Dy  = sparse(fd_normal(size_y,3,y_coarse,1));
U_d_y(:,1:k_d) = Dy*reshape(M_downsampled.rho(1,:,1:k_d),[size_y,k_d]);
[U_d_y,~,~]    = svd(U_d_y,'econ'); % orthonormalize U_x
U_d_y = U_d_y(2:1:end-1,:); % first & last row zero

V_d_y = zeros(k_d,size_t);

q_y = [U_d_y(:);V_d_y(:)];
Flag.y = [length(U_d_y(:)), length(V_d_y(:))];
%%
q = [q_x; q_y];
%%
% Limit the U_d and V_d elements
lb = [];
ub =  [];
%% Regularization terms (Smoothness)
Gamma_x = sparse(fd_normal(size_x, 3, x_coarse, 2));
Gamma_y = sparse(fd_normal(size_y, 3, y_coarse, 2));
Gamma_t = sparse(fd_normal(size_t, 3, t_coarse, 2));
%% Minimun grid spacing
min_grid_spacing_fine_factor = 10;
min_grid_spacing_fine = 1/Nx/min_grid_spacing_fine_factor;
%% ---------------------- LSQ-nonlinear
delete(gcp('nocreate'))
options = optimoptions(...
    @lsqnonlin,...%     'lsqnonlin',
    'Algorithm',    'levenberg-marquardt',...
    'Display',      'iter-detailed',...
    'MaxFunEvals',  1e6,...
    'MaxIter',      1,...
    'ScaleProblem', 'Jacobian',...
    'FinDiffType',  'forward',...
    'PlotFcn',      @optimplotresnorm,... % plotting let us stop/pause lsqonlin at any point
    'StepTolerance',1e-12,...
    'FunctionTolerance',1e-12,...
    'OptimalityTolerance',1e-12,...
    'TolFun',       1e-12,...
    'TolX',         1e-12);%,...
    %     'InitDamping',  1);% options.SubproblemAlgorithm='cg'
if and(par ~= 0,~isempty(par)); options.UseParallel=true; end
%%
options.MaxIter = 20;
F = @(q) objective_ALE_multiv_XY(q,...
    kM, Flag, ...
    x_coarse, y_coarse, M_downsampled,...
    Gamma_0_M, Gamma_x, Gamma_y, Gamma_t);
[q, resnorm_old] = lsqnonlin(F,q,lb,ub,options); % added to drop the first iterations from the Residual plot
save(['opt_',file_name], '-v7.3')
%% Optimizer and auto-save
options.MaxIter = SaveIter;
exitFlag_autosave = 0;
for loop_count = 1:1:floor(MaxIter/SaveIter)
    if exitFlag_autosave == 0
        
        display(['loop #',num2str(loop_count)])
        
        penalty_in_min_spacing = 10^(loop_count-1);
        
        F = @(q) objective_ALE_multiv_XY_plus_min_spacing(q,...
            kM, Flag, ...
            x_coarse, y_coarse, M_downsampled,...
            Gamma_0_M, Gamma_x, Gamma_y, Gamma_t, ...
            min_grid_spacing_fine, penalty_in_min_spacing, Nx, Ny);%, penalty_in_Dxx_fine, x_0, Nx, Dxx_func);
        
        [q, resnorm] = lsqnonlin(F,q,lb,ub,options);
        
        save(['opt_',file_name,'_G10d00','_',num2str(loop_count)], '-v7.3') % Auto-save
        if resnorm < options.FunctionTolerance; exitFlag_autosave=1 ;end
        if abs(resnorm_old-resnorm) < options.FunctionTolerance; exitFlag_autosave=1 ;end
        resnorm_old = resnorm;
        options.InitDamping= options.InitDamping/10;
    end
end
delete(gcp('nocreate'))