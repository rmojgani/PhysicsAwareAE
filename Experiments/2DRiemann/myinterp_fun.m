function myinterp = myinterp_fun(interp_type)
% Choose.....
% 1. make it easier to uniformly use one type of interpolaion through the
% code
% 2. make it possible to choose the different interpolation functions
% Note: Now, test for 2D
%
switch interp_type
    case 'griddata'
        disp('griddata')
        myinterp = @(xin, yin, vin, xout, yout)  griddata(xin, yin, vin, xout, yout);
    case 'interp2'
        disp('interp2')
        myinterp = @(xin, yin, vin, xout, yout)  interp2(xin, yin, vin, xout, yout);
    otherwise
        disp('Are you sure you have selected the right interpolation? I guess no, so we changed it to default:')
        myinterp = myinterp_fun('griddata');
end
end