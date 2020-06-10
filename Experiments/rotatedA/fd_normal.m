function Dm = fd_normal(N,order,x,derivative)
% Calculates FD derivative matrix. The parameters are:
%  N           mesh size
%  order+1     ODD only!!! (3,5,7,...)
%  x           mesh points
%  derivative  ^th derivative (1,2,...)

NN=N;
for n=1:(order-1)/2
    w=weights(x(n),x(1:order),derivative);
    Dm(n,1:order)=w(derivative+1,:);
end

for n=((order-1)/2)+1:NN-(order-1)/2
    w=weights(x(n),x(n-((order-1)/2):n+((order-1)/2)),derivative);
    Dm(n,n-((order-1)/2):n+((order-1)/2))=w(derivative+1,:);
end
for n=(NN-((order-1)/2))+1:NN
    w=weights(x(n),x(NN-order+1:NN),derivative);
    Dm(n,NN-order+1:NN)=w(derivative+1,:);
end

