% 2018-05-20 12:25:00.909349325 +0200
%% solve SWE for statinary flow (dU/dt = dQ/dt = 0)
function solve_statinary(obj)
	% TODO inifun
	%x0           = [obj.h; obj.qx; obj.qy];
	n  = prod(obj.mesh.n);
	%x0 = [-obj.zb;zeros(2*n,1)];
	x0 = [-obj.zb;1e-7*randn(2*n,1)];
	sopt         = struct();
	sopt.abstol  = 1e-3;
	sopt.maxiter = 20;
	
	% TODO use gauss-newton
	% solve non-linear system iteratively
	[x, cflag] = picard(@(x) step(x), x0, sopt);
	obj.h  = x(1:n);
	obj.qx = x(n+1:2*n);
	obj.qy = x(2*n+1:3*n);

function x = step(x)
	x_ = x;
	[A,b] = obj.assemble_stationary(  x(1:n) ...
				 	, x(n+1:2*n) ...
					, x(2*n+1:3*n) ...
					);
	[A,b] = obj.apply_boundary_condition_stationary(A,b,x(1:n),x(n+1:2*n),x(2*n+1:3*n));
	% solve
%	x = gmres(A,b,n);
	x = A\b;
	figure()
	x_ = reshape(x(1:n),obj.mesh.n);
	subplot(2,2,1)
	imagesc(x_)
	x_ = reshape(x(n+1:2*n),obj.mesh.n);
	subplot(2,2,2)
	imagesc(x_)
	x_ = reshape(x(2*n+1:3*n),obj.mesh.n);
	subplot(2,2,3)
	imagesc(x_)
%	obj.mesh.plot(x(1:n));
%	plot([x_(1:n),x(1:n)])
end % step

end % solve_stationary

