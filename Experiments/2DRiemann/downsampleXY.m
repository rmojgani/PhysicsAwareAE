function [X_downsampled, x_coarse, y_coarse, t_coarse, size_x, size_y, size_t, X_c, Y_c] = ...
    downsampleXY(X, x_fine, y_fine, t_fine, xbox_size, ybox_size, tbox_size)

x_coarse = linspace(0, x_fine(end), xbox_size+1)' ; % Divides x in xbox_size boxes, so we get xbox_size+1 dividers
y_coarse = linspace(0, y_fine(end), ybox_size+1)' ; % Divides y in ybox_size boxes, so we get ybox_size+1 dividers
t_coarse = linspace(0, t_fine(end), tbox_size+1)' ;

[X_f, Y_f, T_f] = meshgrid(  x_fine,   y_fine,   t_fine);
[X_c, Y_c, T_c] = meshgrid(x_coarse, y_coarse, t_coarse);

X = permute(X,[2 1 3]);

X_downsampled = interp3(X_f, Y_f, T_f, X, X_c, Y_c, T_c, 'spline');

X_downsampled = permute(X_downsampled,[2 1 3]);
[size_x , size_y, size_t] = size(X_downsampled);
X_c = X_c(:,:,1)';
Y_c = Y_c(:,:,1)';

end


% 
% x_coarse = linspace(0, x_fine(end), xbox_size+1)' ; % Divides x in xbox_size boxes, so we get xbox_size+1 dividers
% y_coarse = linspace(0, y_fine(end), ybox_size+1)' ; % Divides y in ybox_size boxes, so we get ybox_size+1 dividers
% t_coarse = linspace(0, t_fine(end), tbox_size+1)' ;
% 
% [x_f, y_f, t_f] = meshgrid(  x_fine,   y_fine,   t_fine);
% [x_c, y_c, t_c] = meshgrid(x_coarse, y_coarse, t_coarse);
% 
% for n=1:1:size(t_f,3)
%     X2(:,:,n) = X(:,:,n)';
% end
% 
% X_downsampled = interp3(x_f, y_f, t_f, X2, x_c, y_c, t_c, 'spline');
% 
% [size_y , size_x, size_t] = size(X_downsampled);

% %%
% function [X_downsampled, x_coarse, y_coarse, t_coarse, size_x, size_y, size_t] = downsampleXY(X, x_fine, y_fine, t_fine, xbox_size, ybox_size, tbox_size)
% 
% x_coarse = linspace(0, x_fine(end), xbox_size+1)' ; % Divides x in xbox_size boxes, so we get xbox_size+1 dividers
% y_coarse = linspace(0, y_fine(end), ybox_size+1)' ; % Divides y in ybox_size boxes, so we get ybox_size+1 dividers
% t_coarse = linspace(0, t_fine(end), tbox_size+1)' ;
% 
% [x_f, y_f, t_f] = meshgrid(  x_fine,   y_fine,   t_fine);
% [x_c, y_c, t_c] = meshgrid(x_coarse, y_coarse, t_coarse);
% 
% X_downsampled = interp3(x_f, y_f, t_f, X, x_c, y_c, t_c, 'spline');
% 
% [size_y , size_x, size_t] = size(X_downsampled);
% end