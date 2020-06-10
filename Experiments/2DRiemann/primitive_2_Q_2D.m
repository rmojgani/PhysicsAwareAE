function q = primitive_2_Q_2D(p, u, v, rho, gamma)
q(:,:,1) = rho;
q(:,:,2) = rho.*u;
q(:,:,3) = rho.*v;
q(:,:,4) = p./(gamma-1) + 0.5*rho.*(u.^2 + v.^2);
end