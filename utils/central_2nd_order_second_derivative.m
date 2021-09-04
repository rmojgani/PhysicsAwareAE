function D = central_2nd_order_second_derivative(x,PeriodicFlag)

	% Returns the second-order, central finite difference 
	% approximation of the secon derivative

Nx = max(size(x));
dx = x(2) - x(1);

V0 = ones([1,Nx-0]);
V1 = ones([1,Nx-1]);
V2 = ones([1,Nx-2]);
V3 = ones([1,Nx-3]);

D = (diag(V1,-1) -2*diag(V0,0) + diag(V1,1))/(dx^2);

if(PeriodicFlag == 0)
	D(1,:)	= [1 -2 1 0*V3]/(dx^2);
	D(Nx,:)	= [0*V3 1 -2 1]/(dx^2);
elseif(PeriodicFlag == 1) 
	D(1,:)	= [-2 1 0*V3 1]/(dx^2);
	D(Nx,:)	= [1 0*V3 1 -2]/(dx^2);
else
	error('Incorrect PeriodicFlag');
end
