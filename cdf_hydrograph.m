function F = cdf_hydrograph(Q,Qmin,Qmax)
	F = 1/(pi)*acos(1-2*(Q-Qmin)/(Qmax-Qmin));
end

