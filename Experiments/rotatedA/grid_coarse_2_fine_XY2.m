function [U_d_x_out, V_d_x_out, U_d_y_out, V_d_y_out, x_morphing_out, y_morphing_out] = ...
    grid_coarse_2_fine_XY2(U_d_x_in, V_d_x_in, U_d_y_in, V_d_y_in, x_in, y_in, t_in, x_out, y_out, t_out, interp)

U_d_x_out = interp1(x_in, U_d_x_in , x_out, interp) ;
V_d_x_out = interp1(t_in, V_d_x_in', t_out, interp)';

U_d_y_out = interp1(y_in, U_d_y_in , y_out, interp) ;
V_d_y_out = interp1(t_in, V_d_y_in', t_out, interp)';


d_x_out = U_d_x_out * V_d_x_out;
d_y_out = U_d_y_out * V_d_y_out;

M_x_ALE = bsxfun(@plus, x_out, d_x_out);
M_y_ALE = bsxfun(@plus, y_out, d_y_out);

%     x_m = zeros(size_y, size_x, size_t);
%     y_m = zeros(size_y, size_x, size_t);
for tcount = 1:1:length(t_out)
    [x_morphing_out(:,:,tcount), y_morphing_out(:,:,tcount)] = meshgrid(M_x_ALE(:,tcount), M_y_ALE(:,tcount));
end

% x_morphing_out = bsxfun(@plus, x_out, d_x_out);
% y_morphing_out = bsxfun(@plus, y_out, d_y_out);

end