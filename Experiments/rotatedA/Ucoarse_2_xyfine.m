function [x,y] = Ucoarse_2_xyfine(U_x_pos, V_x_pos, U_y_pos, V_y_pos,...
    size_x,size_y,size_t,n_pos,...
    Nx, Ny, ...
    x_e, y_e, x_e_c, y_e_c)
U_x = reshape(U_x_pos, size_x*size_y, n_pos);
V_x = reshape(V_x_pos, n_pos        , size_t);
U_y = reshape(U_y_pos, size_x*size_y, n_pos);
V_y = reshape(V_y_pos, n_pos        , size_t);%V_x;% V_y_pos = reshape(q( sum(Flag(1:3))+1:sum(Flag(1:4))    ,1), n_pos, Nt);

d_x = U_x * V_x;
d_y = U_y * V_y;

%  approach I : interpolating the 2d grid
for icount=1:1:size_t
    d_x_fine_m = interp2(x_e_c,y_e_c,reshape(d_x(:,icount),[size_x,size_y]),x_e,y_e);
    d_y_fine_m = interp2(x_e_c,y_e_c,reshape(d_y(:,icount),[size_x,size_y]),x_e,y_e);
    d_x_fine(:,icount) = d_x_fine_m(:);
    d_y_fine(:,icount) = d_y_fine_m(:);
end
%

x = bsxfun(@plus, d_x_fine, reshape(x_e,[Nx*Ny,1]) );
y = bsxfun(@plus, d_y_fine, reshape(y_e,[Nx*Ny,1]) );

% approach II (not implemented):
% U_x_fine = interp1(linspace(0,1,size_x*size_y),U_x,linspace(0,1,Nx*Ny),'linear');
% V_x_fine = interp1(linspace(0,1,size_t)       ,V_x',linspace(0,1,Nt)  ,'linear')';
% U_y_fine = interp1(linspace(0,1,size_x*size_y),U_y,linspace(0,1,Nx*Ny),'linear');
% V_y_fine = interp1(linspace(0,1,size_t)       ,V_y',linspace(0,1,Nt)  ,'linear')';
% 
% 
% U_x_pos = U_x_fine;
% V_x_pos = V_x_fine;
% U_y_pos = U_y_fine;
% V_y_pos = V_y_fine;
% 
% d_x_fine = U_x_fine * V_x_fine;
% d_y_fine = U_y_fine * V_y_fine;
end