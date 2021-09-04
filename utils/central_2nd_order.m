function D = central_2nd_order(x,PeriodicFlag)

	% Returns the second-order, central finite difference 
	% approximation of the first derivative

Nx = max(size(x));
dx = x(2) - x(1);

V1 = ones([1,Nx-1]);
V2 = ones([1,Nx-2]);
V3 = ones([1,Nx-3]);

D = (-diag(V1,-1) + diag(V1,1))/(2*dx);

if(PeriodicFlag == 0)
	D(1,:)	= [-1 1 0*V2]/dx;
	D(Nx,:)	= [0*V2 -1 1]/dx;
elseif(PeriodicFlag == 1) 
	D(1,:)	= [0 1 0*V3 -1]/(2*dx);
	D(Nx,:)	= [1 0*V3 -1 0]/(2*dx);
else
	error('Incorrect PeriodicFlag');
end
