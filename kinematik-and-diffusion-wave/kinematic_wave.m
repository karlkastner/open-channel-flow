% Thu 26 Jan 17:42:56 CET 2023
function eh = kinematic_wave(t,L,eh,h0,C0,S0,dflag,bcfun)
	if (nargin()<7)
		dlag = true;
	end
	nt = length(t);
	nx = length(eh);
	u0 = C0*sign(S0).*sqrt(h0*abs(S0));
	eh = [eh,zeros(nx,nt-1)];
	D1 = derivative_matrix_1_1d(nx,L,+sign(S0),'dirichlet');
	D2 = derivative_matrix_2_1d(nx,L,2,'dirichlet');
	dx = L/nx;
	if (abs(u0)*(t(2)-t(1)) > dx)
		warning("time step limit exeeded");
	end

	for idx=2:length(t)
		dt = t(idx)-t(idx-1);
		if (dflag)
			d = h0*u0/(2*S0)*(D2*eh(:,idx-1));
		else
			d = 0;
		end
		eh(:,idx) = eh(:,idx-1) - dt*(d + 3/2*u0*(D1*eh(:,idx-1)));
		if (nargin()>7)
			bcval = bcfun(t(idx));
			eh(1,idx) = bcval;
		end
	end
end

