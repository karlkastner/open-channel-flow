% 2020-08-25 16:49:40.529102696 +0800
function F = cdf_hydrograph(Q,Qmin,Qmax)
	F = 1/(pi)*acos(1-2*(Q-Qmin)/(Qmax-Qmin));
end

