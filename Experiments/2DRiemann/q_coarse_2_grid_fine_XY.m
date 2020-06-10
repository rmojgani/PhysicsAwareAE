function [x_morphing_fine , y_morphing_fine, ...
        x_f             , y_f        ,...    
        U_d_x_fine      , V_d_x_fine , ...
        U_d_y_fine      , V_d_y_fine ] = q_coarse_2_grid_fine_XY(size_x, size_y, size_t, q_grid, ...
                                                                 x_coarse, y_coarse, t_coarse, x_fine, y_fine, t_fine, ...
                                                                 k_d, Flag, if_sym, coarse_to_fine_grid_interp)
    size_x = size_x-2; % first row zero
    size_y = size_y-2; % first row zero
    length_q_x = (size_x + size_t)*k_d;
    %% ---------------------- Up-sampling the grid
if if_sym
    %%  Symmetric grid
    [ U_d_x , V_d_x , ~] = q2UV( q_grid(1:1:length_q_x), size_x, size_t, k_d , Flag.x);
    U_d_x = padarray(U_d_x,1,0,'both'); % first & last row zero
%     size_x = size_x+2; % first & last row zero - added row
%     d_x = U_d_x * V_d_x;
    U_d_y = U_d_x;
    V_d_y = V_d_x;
   
elseif ~if_sym
    %%  Assymmetric grid
    [ U_d_x , V_d_x , ~] = q2UV( q_grid(1:1:length_q_x), size_x, size_t, k_d , Flag.x);
    U_d_x = padarray(U_d_x,1,0,'both'); % first & last row zero
%     size_x = size_x+2; % first & last row zero - added row
%     d_x = U_d_x * V_d_x;

    [ U_d_y , V_d_y , ~] = q2UV( q_grid(length_q_x+1:1:end), size_y, size_t, k_d , Flag.y);
    U_d_y = padarray(U_d_y,1,0,'both'); % first & last row zero
%     size_y = size_y+2; % first & last row zero - added row
%     d_y = U_d_y * V_d_y;

end
[U_d_x_fine, V_d_x_fine, U_d_y_fine, V_d_y_fine, x_morphing_fine, y_morphing_fine] = ...
    grid_coarse_2_fine_XY(U_d_x, V_d_x, U_d_y, V_d_y, ...
    x_coarse, y_coarse, t_coarse, ...
    x_fine  , y_fine  , t_fine  , coarse_to_fine_grid_interp);
[x_f, y_f] = meshgrid(x_fine, y_fine);