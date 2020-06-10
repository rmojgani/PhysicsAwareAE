function [f, M_x_ALE, M_y_ALE, size_x, size_y] = objective_ALE_multiv_XY(q,...
    kM, Flag, ...
    x_0_downsampled, y_0_downsampled, X_E_FOM_multi, ...
    Gamma_0_M, Gamma_x, Gamma_y, Gamma_t)
fields = fieldnames(X_E_FOM_multi);

k_d = kM(1);
k_w = kM(2);
size_x = size(X_E_FOM_multi.(fields{1}),1);
size_x = size_x-2; % first row zero

size_y = size(X_E_FOM_multi.(fields{1}),2);
size_y = size_y-2; % first row zero

size_t = size(X_E_FOM_multi.(fields{1}),3);

length_q_x = (size_x + size_t)*k_d;
%%
[ U_d_x , V_d_x , ~] = q2UV( q(1:1:length_q_x), size_x, size_t, k_d , Flag.x);
U_d_x = padarray(U_d_x,1,0,'both'); % first & last row zero
size_x = size_x+2; % first & last row zero - added row
d_x = U_d_x * V_d_x;
%%
[ U_d_y , V_d_y , ~] = q2UV( q(length_q_x+1:1:end), size_y, size_t, k_d , Flag.y);
U_d_y = padarray(U_d_y,1,0,'both'); % first & last row zero
size_y = size_y+2; % first & last row zero - added row
d_y = U_d_y * V_d_y; %d_x;% 
%%
% Preallocations
% M_w_ALE_Proj_on_Eulerian =  zeros(numel(fields)*size_x, size_y, size_t );
% X_E_FOM = zeros(numel(fields)*size_x, size_y, size_t );
% M_w_ALE_XY = zeros(size_x, size_y, size_t);
% M_w_ALE_Proj_on_Eulerian_multi_XY = zeros(size_x, size_y, size_t);
% x_m = zeros(size_x, size_y, size_t);
% y_m = zeros(size_x, size_y, size_t);
% --------------
for varcount = 1:numel(fields)
    %% ---------------------- ALE morphing grid, solved with interpolation
    [x_f, y_f] = meshgrid(x_0_downsampled, y_0_downsampled);
    % find w on reduced morphing grid
    M_x_ALE = bsxfun(@plus, x_0_downsampled, d_x);
    M_y_ALE = bsxfun(@plus, y_0_downsampled, d_y);

    X_E_FOM_multi.(fields{varcount}) = permute(X_E_FOM_multi.(fields{varcount}) ,[2 1 3]);

    for tcount = 1:1:size_t
        [x_m(:,:,tcount), y_m(:,:,tcount)] = meshgrid(M_x_ALE(:,tcount), M_y_ALE(:,tcount));
        M_w_ALE_XY(:,:,tcount) = griddata( x_f, y_f, X_E_FOM_multi.(fields{varcount})(:,:,tcount),...
            x_m(:,:,tcount), y_m(:,:,tcount));%, 'linear', 0);
    end
    M_w_ALE_XY(isnan(M_w_ALE_XY))=10; % Remove NaN valaues caused by griddata and replace with a overshoot
    %% ---------------------- Projection of ALE
    M_w_ALE = reshape( M_w_ALE_XY, [size_x*size_y,size_t]);
    [U_w, S_w ,V_w ]  = svd(M_w_ALE, 'econ');
    [U_w, V_w] = USV_reduce(U_w, S_w ,V_w, k_w);
    M_w_ALE_Proj = U_w*V_w;
    %% ---------------------- ALE_Proj on Eulerian Grid
    M_w_ALE_Proj_XY = reshape( M_w_ALE_Proj, [size_y, size_x, size_t]);
    % --------------
    for tcount = 1:1:size_t
        M_w_ALE_Proj_on_Eulerian_multi_XY(:,:,tcount) = ...
        griddata( x_m(:,:,tcount), y_m(:,:,tcount), M_w_ALE_Proj_XY(:,:,tcount), x_f, y_f);%, 'linear', 0);
    end
    M_w_ALE_Proj_on_Eulerian((varcount-1)*size_y+1:1:varcount*size_y,:,:) = M_w_ALE_Proj_on_Eulerian_multi_XY;
    X_E_FOM((varcount-1)*size_y+1:1:varcount*size_y,:,:) = X_E_FOM_multi.(fields{varcount});
end
%%x
% Distance between FOM and Proj
f1 = X_E_FOM - M_w_ALE_Proj_on_Eulerian; % w in equation
% Smoothness of basis
f2 = Gamma_0_M(1).*[Gamma_x*U_d_x ; Gamma_y*U_d_y ];
f3 = Gamma_0_M(2).*[Gamma_t*V_d_x'; Gamma_t*V_d_y'];

f = [ f1(:); f2(:); f3(:)];
end