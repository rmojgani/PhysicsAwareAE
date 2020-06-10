p = (gamma-1)*(q(:,:,4) - 0.5*(q(:,:,2).^2 + q(:,:,3).^2)./q(:,:,1));
c = sqrt(gamma*p./q(:,:,1));

nu =  dx/2;%0.005 + 0.5*0.5*1e5*(dx.^5).*(abs(Dxx*p) + abs(p*Dyy'));%((dx.^5)./2).*(abs(Dxxxx*p) + abs(p*Dyyyy'));

spx = c + abs(q(:,:,2)./q(:,:,1));
spy = c + abs(q(:,:,3)./q(:,:,1));
CFL = 0.25;
dt1 = min(min(min([dx./spx dy./spy])));
dt2 = min(min(min([dx dy])))^2/(2*max(max(nu))); 
dt = CFL*min(dt1,dt2);
