function [f] = objective_rotA_proj2(q,...
        k, Flag, ...%x_coarse, y_coarse, ...
        M, ...
        x_e     , y_e     , ...
        x_e_c   , y_e_c   , ...
        Nx    , Ny    , Nt    ,...
        size_x, size_y, size_t,...
        Gamma_x, Gamma_t, ...
        griddata_method)
n_pos = k(1);
n_mag = k(2);

length_q.x = (size_x*size_y + size_t) * n_pos;

[ U_x_pos , V_x_pos , ~] = q2UV( q(           1:1:length_q.x), size_x*size_y, size_t, n_pos , Flag.x);
[ U_y_pos , V_y_pos , ~] = q2UV( q(length_q.x+1:1:end       ), size_x*size_y, size_t, n_pos , Flag.y);
%V_x_pos(:,1)=0;
%V_y_pos(:,1)=0;
%%
[x,y] = Ucoarse_2_xyfine(U_x_pos, V_x_pos, U_y_pos, V_y_pos,...
    size_x,size_y,size_t,n_pos,...
    Nx, Ny, ...
    x_e, y_e, x_e_c, y_e_c);

% M_tilde = zeros(Nx, Ny, Nt); % Pre-allocation
M_tilde_lag = zeros(Nx*Ny , Nt); % Pre-allocation
for j = 1:Nt
    xx = reshape(x(:,j),[Nx,Ny]);
    yy = reshape(y(:,j),[Nx,Ny]);
    M_tilde_j = griddata_fill(x_e, y_e, M(:,:,j), ...
                                         xx, yy, griddata_method);
    M_tilde_lag(:,j) = M_tilde_j(:);
end

[U_lag,S_lag,V_lag]=svd(M_tilde_lag,'econ');

M_tilde_lag = U_lag(:,1:n_mag)*S_lag(1:n_mag,1:n_mag)*V_lag(:,1:n_mag)';

M_tilde_lag_back = zeros(Nx, Ny , Nt); % Pre-allocation
for j = 1:Nt
    xx = reshape(x(:,j),[Nx,Ny]);
    yy = reshape(y(:,j),[Nx,Ny]);
    ww_lag = reshape(M_tilde_lag(:,j),[Nx,Ny]);
    M_tilde_lag_back(:,:,j)  = griddata_fill(xx, yy, ww_lag, ...
                                             x_e, y_e, griddata_method);
end

f_1 = M - M_tilde_lag_back;
% f_1 =  rho_GM(M - M_tilde_lag_back);
f_2 = Gamma_x*U_x_pos;
f_3 = Gamma_x*U_y_pos;
f_4 = Gamma_t*V_x_pos';
f_5 = Gamma_t*V_y_pos';

f = [f_1(:);f_2(:);f_3(:);f_4(:);f_5(:)];
%f = [f_1(:);f_4(:);f_5(:)];%f_4(:);f_5(:)];
end
