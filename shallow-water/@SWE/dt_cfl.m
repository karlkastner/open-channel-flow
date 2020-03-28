% 2017-11-05 15:02:25.431518701 +0100
%% determine time step required by cfl
function dt = dt_cfl(obj,q,dx)
	g = obj.g;

	n = length(q)/2;


	if (obj.QAflag)	
		A = q(1:n);
		Q = q(n+1:2*n);
		h = A./obj.w;
		u = Q./A;
	else
		h = q(1:n);
		u = q(n+1:2*n)./h;
	end

	%h_max = max(abs(h));
	%u_max = max(u);

	% maximum celerity
	% safety
	h    = abs(h);
	cmax = max(abs(u) + sqrt(g*h));

	% maximum time step statisfying CFL stability condition
	dt = dx/cmax;
end

