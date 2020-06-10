%% ************************************************************************
M_lag_j = zeros(Nx*Ny , 1); % Pre-allocation
M_lag_vec = zeros(Nx*Ny , Nt); % Pre-allocation
for j = 1:Nt
    x_m = reshape(x_m_v(:,j),[Nx,Ny]);
    y_m = reshape(y_m_v(:,j),[Nx,Ny]);
    M_lag_j = griddata_fill(x_e, y_e, M(:,:,j), ...
                              x_m, y_m, griddata_method);
    M_lag_vec(:,j) = M_lag_j(:); 
end
%% [SVD or simply interpolation] of snapshot on moving grid
% M_tilde_lag_vec = bsxfun(@times, M_lag_vec(:,1), ones(Nx*Ny, Nt) );
[U_lag, S_lag, V_lag] = svd(M_lag_vec,'econ');
[U_lag, V_lag ] = USV_reduce(U_lag, S_lag ,V_lag, n_mag);
M_tilde_lag_vec = U_lag * V_lag;
%% Reconstructed (snapshot on moving grid) interpolated on Eulerian grid
M_tilde_moving = zeros(Nx, Ny , Nt); % Pre-allocation
for j = 1:Nt
    x_m = reshape(x_m_v(:,j),[Nx,Ny]);
    y_m = reshape(y_m_v(:,j),[Nx,Ny]);
    w_lag_j = reshape(M_tilde_lag_vec(:,j),[Nx,Ny]);
    M_tilde_moving(:,:,j)  = griddata_fill(x_m, y_m, w_lag_j, ...
                                             x_e, y_e, griddata_method);
end