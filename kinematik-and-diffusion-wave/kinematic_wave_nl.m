% Wed 25 Jan 09:41:16 CET 2023
% TODO for 2D do it for every column
% TODO rename, as this is actually a nl-diffusion wave
function h = kinematic_wave_nl(t,x,h,zb, C, a, p);
	explicit = true;
	nx = length(x);
	nt = length(t);
	h = [cvec(h),zeros(nx,nt-1)];
	for idx=2:nt
		dt = t(idx)-t(idx-1);
		D1 = advection_nl_matrix(x,h(:,idx-1),zb,C);
		if (explicit)
			h(:,idx) = h(:,idx-1) + dt*(p - a.*h(:,idx-1) + D1*h(:,idx-1));
		else % implicit
			I = speye(length(x));
			h = dt*p + (I - dt*(diag(sparse(a)) + D))\h;
		end
	end % for idx
end % kineamtic_wave_nl


