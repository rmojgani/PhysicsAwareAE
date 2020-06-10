function M_d = griddata_fill(x_s, y_s, v_on_s, x_d, y_d, griddata_method)
% x_s and y_s are grid of source
% v_s is the data on source
% x_d and y_d are grid of destination

M_d = griddata(x_s, y_s, v_on_s, ...
                x_d, y_d, griddata_method);
                        
M_d( isnan(M_d) ) = 1;

end