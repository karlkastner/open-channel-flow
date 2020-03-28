% Sat 19 May 13:43:45 CEST 2018
%
%% infer the bed level
function obj = infer_bed_level2(obj)
	x    = obj.h;
	n    = obj.mesh.n;
	iter = 0;
	fun  = @(x) obj.objective_bed_level(x);
	while (1)
		iter = iter+1;
		[f0,A] = gradpde2d(fun,x,n);
		sum(f0)
		%x = x - A \ f0;
		%A = A + 1e-6*speye(size(A));
		%A = diag(diag(A));
		p = 1e-2;
		x = x + p*minres(A,f0,[],prod(n));
		figure(iter)
		clf();
%		obj.plot(diag(A));
		obj.plot(real(x))
		if (iter >= 4)
			break;
		end
	
	end % while
	
	obj.h = x;
end % infer_bed_level2

