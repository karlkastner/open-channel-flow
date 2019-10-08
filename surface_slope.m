% Wed  2 May 14:20:51 CEST 2018
%% surface slope for uniform stationary flow
function dzs_dx = surface_slope(Q,H,W,C)
	dzs_dx = Q.^2./(C.^2.*W.^2.*H.^3);
end

