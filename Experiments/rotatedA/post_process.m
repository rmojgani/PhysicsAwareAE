%% Load file
file_name = 'A50CMAX90|3_n_mag_1_n_pos_2_cp7x7|G100_100|00XX.mat';
load(file_name)
%%
n_mag = 1;
length_q.x = (size_x*size_y + size_t) * n_pos;

[ U_x_pos , V_x_pos , ~] = q2UV( q(           1:1:length_q.x), size_x*size_y, size_t, n_pos , Flag.x);
[ U_y_pos , V_y_pos , ~] = q2UV( q(length_q.x+1:1:end       ), size_x*size_y, size_t, n_pos , Flag.y);
% V_y_pos = V_x_pos;
% V_x_pos(:,1)=0;
% V_y_pos(:,1)=0;

[x,y] = Ucoarse_2_xyfine(U_x_pos, V_x_pos, U_y_pos, V_y_pos,...
    size_x,size_y,size_t,n_pos,...
    Nx, Ny,...
    x_e, y_e, x_e_c, y_e_c);

M_tilde_lag = zeros(Nx*Ny , Nt); % Pre-allocation
for j = 1:Nt
    xx = reshape(x(:,j),[Nx,Ny]);
    yy = reshape(y(:,j),[Nx,Ny]);
    M_tilde_j = griddata(x_e,y_e,M(:,:,j),xx,yy,'linear');
    M_tilde_j(isnan(M_tilde_j))=1;
    M_tilde_lag(:,j) = M_tilde_j(:); % --- 
end
% %%
[U_lag,S_lag,V_lag]=svd(M_tilde_lag,'econ');

M_tilde_lag = U_lag(:,1:n_mag)*S_lag(1:n_mag,1:n_mag)*V_lag(:,1:n_mag)';
% %%
M_tilde_lag_back = zeros(Nx, Ny , Nt); % Pre-allocation
for j = 1:Nt
    xx = reshape(x(:,j),[Nx,Ny]);
    yy = reshape(y(:,j),[Nx,Ny]);
    ww_lag = reshape(M_tilde_lag(:,j),[Nx,Ny]);
    M_tilde_j_back = griddata(xx,yy,ww_lag,x_e,y_e,'linear');
    M_tilde_j_back(isnan(M_tilde_j_back))=1;
    M_tilde_lag_back(:,:,j) =  M_tilde_j_back;
end
% %%
M_tilde_moving=reshape(M_tilde_lag_back, [Nx,Ny, Nt]);
% %% ROM on Eulerian/Stationary Grid
[U,S,V] = svd(MC,'econ');
U = U(:,1:n_mag);
V = S(1:n_mag,1:n_mag)*V(:,1:n_mag)';

UV = U*V;
M_tilde = reshape(UV,[Nx,Ny,Nt]);
%%
matrix_subplotter