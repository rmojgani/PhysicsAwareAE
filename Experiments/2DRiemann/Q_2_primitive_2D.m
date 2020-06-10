function [p, u, v, rho] = Q_2_primitive_2D(q, gamma)
rho = q(:,:,1);
u = q(:,:,2)./q(:,:,1);
v = q(:,:,3)./q(:,:,1);
p = ( gamma - 1 ) * ( q(:,:,4) - 0.5 * rho .* ( u.^2 + v.^2 ) );
end