% 2020-08-25 22:21:28.352073422 +0800
% inverse cumulative distribution (quantiles)
% of a sinusoidal hydrograph
function Q = inv_hydrograph(F,Qmin,Qmax)
	Q = 0.5*(Qmin+Qmax) + 0.5*(Qmax-Qmin)*quantile_sin(F);
end
