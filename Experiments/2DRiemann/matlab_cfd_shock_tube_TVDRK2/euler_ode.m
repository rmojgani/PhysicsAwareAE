function r = euler_ode(q, Nx, Ny)

global gamma dx dy
% global Nx Ny Nz
global Dx Dy Dxx Dyy Dxxxx Dyyyy

N = Nx*Ny;

q(:,2    ,:)  = q(:,1  ,:);
q(:,end-1,:)  = q(:,end,:);

q(2    ,:,:)  = q(1    ,:,:);
q(end-1,:,:)  = q(end,:,:);

p = (gamma - 1)*(q(:,:,4) - 0.5*(q(:,:,2).^2 + q(:,:,3).^2)./q(:,:,1));

% nu = 1.3*(dx/2);%1e-4*(dx/2).*(abs(Dxx*p) + abs(p*Dyy'));% ;%0.005 + 0.5*0.5*1e5*(dx.^5).*(abs(Dxx*p) + abs(p*Dyy'));% 
% nu = 0.001;%(dx.^5).*(abs(Dxxxx*p) + abs(p*Dyyyy'))+0.0001;
nu = (dx.^5).*(abs(Dxxxx*p) + abs(p*Dyyyy'));%+0.0001;

r = zeros(Nx,Ny,4);

r(:,:,1) =     - Dx*(  q(:,:,2)                             );
r(:,:,2) =     - Dx*( (q(:,:,2).*q(:,:,2))./q(:,:,1)  + p   );
r(:,:,3) =     - Dx*( (q(:,:,2).*q(:,:,3))./q(:,:,1)        );
r(:,:,4) =     - Dx*( (q(:,:,2)./q(:,:,1)).*(q(:,:,4) + p)  );


r(:,:,1) = r(:,:,1) - (  q(:,:,3)                             )*Dy';
r(:,:,2) = r(:,:,2) - ( (q(:,:,3).*q(:,:,2))./q(:,:,1)        )*Dy';
r(:,:,3) = r(:,:,3) - ( (q(:,:,3).*q(:,:,3))./q(:,:,1)  + p   )*Dy';
r(:,:,4) = r(:,:,4) - ( (q(:,:,3)./q(:,:,1)).*(q(:,:,4) + p)  )*Dy';

for n = 1:4
  r(:,:,n) = r(:,:,n) + nu.*(Dxx*q(:,:,n) + q(:,:,n)*Dyy');
end
