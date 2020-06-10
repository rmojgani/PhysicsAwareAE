function [M, M_lag, x_m_vec, y_m_vec]  = imrotate_via_griddata(Img, angleM, griddata_method)
% Input :             Img: Image on Eulerian grid
%         griddata_method:
% -------------------------------------------------------------------------
% Output:               M: Snapshot of rotated image on Eulerian grid
%                 x_m_vec: x, moving grid in vector form
%                 y_m_vec: y, moving grid in vector form
% =========================================================================
[Nx, Ny] = size(Img);    % snapshot size
Nt = length(angleM);

[x_e, y_e] = meshgrid(1:1:Nx, 1:1:Ny);
%% Build the rotated grid
x_m_vec = zeros( Nx*Ny, Nt ); % Pre-allocation
y_m_vec = zeros( Nx*Ny, Nt ); % Pre-allocation
for icount = 1:1:length(angleM)
    [Xr, Yr] = rotate2Dmesh( x_e, y_e, -angleM(icount) );
    x_m_vec(:,icount) = Xr(:);
    y_m_vec(:,icount) = Yr(:);
end
%%
M_lag = zeros( Nx, Ny , Nt); % Pre-allocation
for j = 1:1:Nt
    x_m = reshape( x_m_vec(:,j), [ Nx, Ny] );
    y_m = reshape( y_m_vec(:,j), [ Nx, Ny] );
    M_lag(:,:,j) = griddata_fill(x_m, y_m, Img, ...
                                      x_e, y_e, griddata_method);
end
%%
M = reshape( M_lag, [ Nx, Ny, Nt] );
end