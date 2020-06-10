function pr = rotate2D(p, theta, center)

if nargin <= 2
    pr = my_rotz(theta)*[p;0];
else
    p = p - center;
    pr = rotate2D(p, theta);
    pr = pr + center;
end
pr = pr(1:2);
end