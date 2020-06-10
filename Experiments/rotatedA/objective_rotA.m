function [f] = objective_rotA(q,k,Nx,Ny,Nt,x_e,y_e,Gamma_x,Gamma_t,Flag,M)
  
U_w = reshape(q(            1:Flag(1)                   ,1),Nx*Ny,k);
V_w = reshape(q(           Flag(1)+1:sum(Flag(1:2))     ,1),k,Nt);

U_x = reshape(q(    sum(Flag(1:2))+1:sum(Flag(1:3))     ,1),Nx*Ny,k);
V_xy = reshape(q(    sum(Flag(1:3))+1:sum(Flag(1:4))     ,1),k,Nt);

U_y = reshape(q(    sum(Flag(1:4))+1:sum(Flag(1:5))     ,1),Nx*Ny,k);


x = bsxfun(@plus, U_x*V_xy, reshape(x_e,[Nx*Ny,1]) );
y = bsxfun(@plus, U_y*V_xy, reshape(y_e,[Nx*Ny,1]));
w = bsxfun(@plus, U_w*V_w, 0*reshape(M(:,:,1),[Nx*Ny,1]));

M_tilde = zeros(Nx, Ny, Nt); % Pre-allocation
for j = 1:Nt
    xx = reshape(x(:,j),[Nx,Ny]);
    yy = reshape(y(:,j),[Nx,Ny]);
    ww = reshape(w(:,j),[Nx,Ny]);
    M_tilde_j = griddata(xx,yy,ww,x_e,y_e,'linear');
    M_tilde_j(isnan(M_tilde_j))=1;
    M_tilde(:,:,j) = M_tilde_j;
end

f_1 = M - M_tilde;
f_2 = Gamma_x*U_x;
f_3 = Gamma_x*U_y;
f_4 = Gamma_t*V_xy';
f = [f_1(:);f_2(:);f_3(:);f_4(:)];
% f =  norm(f_1(:));
end

