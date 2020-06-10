function [U, V] = USV_reduce(U, S, V, k)
V = S * V';
U = U(  :  , 1:k);
V = V(1:k, :    );
end