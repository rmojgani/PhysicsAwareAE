function [Xr, Yr] = rotate2Dmesh(X,Y,theta)
% Rotate 2D Grid around center point

center = [mean(mean(X)); mean(mean(Y))];

% Pre-allocate
Xr = 0*X;
Yr = 0*Y;
% 
for icount = 1:1:size(X,1)
    for jcount = 1:1:size(X,2)
         pr = rotate2D([X(icount,jcount);Y(icount,jcount)], theta, center);
         Xr(icount,jcount) = pr(1);
         Yr(icount,jcount) = pr(2);
    end
end
end