% Sun 10 Nov 12:12:53 +08 2019
function [zs] = solve_matrix(obj,x,zs,Q0,Qt,Qmid,Qhr,chezy,width,zb,x0,z0)
	if (isempty(zs))
		h = zs-zb;
		if (x0 == x(1))
			h = z0-zb(1);
		else
			h = z0-zb(end);
		end
	end
	
	% dh_dx  = (S_f - S_b)./(1 - F2)
	% dzs_dx = S_b + (S_f - S_b)./(1 - F2);

	D1     = derivative_matrix_1_1d(x,[],-1);
	dzb_dx = D1*zb;
	% at segment centre
	dzb_dx = dzb_dx(2:end);

	h_    = mid(h);
	Qmid_ = mid(Qmid); 
	Qhr_  = mid(Qhr); 
	%Q0 = mid(Qmid); 
	Qt_    = mid(Qt); 
	C_     = mid(C); 
	W_     = mid(W);

	zs_old = zs;

	for idx=1:10

	dh_dx = obj.dh_dx_(h_,Q0,Qt_,Qmid_,Qhr_,C_,W_,dzb_dx);
	rhs = [z0; dhdx-dzb_dx];

	% boundary conditions
	switch (x0)
	case {x(1)}
		A(1,1:2) = [1,0];
	case {x(end)}
		% note : permute first and last row for faster computation
		A(1,1:2) = 0;
		A(1,end) = 1;
	otherwise
		error('here');
	end
	% solve
	zs = A \ zs;
	norm(zs-zs_old)
	zs_old = zs;
	end
end

