% Sat 19 May 16:22:33 CEST 2018
% quasi time-stepping
function [resnorm, dqs, obj] = infer_bed_level2(obj)
	x    = obj.h;
	n    = obj.mesh.n;
	iter = 1;
	fun  = @(x) obj.objective_bed_level(x);
	maxiter = 1e3;
	sq  = 0;
	p   = 10;
	dqs = fun(x);
	resnorm = rms(dqs);
	p_min = 1;
	while (1)
		% TODO, quick boundary
		x = reshape(x,obj.mesh.n);
		x(:,1) = 15;
		x = flat(x);
		x    = max(0.1,x);
		
		% update
		%while (p>p_min)
			x_   = x + p*dqs;
			dqs  = fun(x_);
		%	dqs = obj.mesh.smooth(dqs);
			res_ = rms(dqs)
		%	if (res_ < resnorm(iter))
		%		break;
		%	end
		%	p = p/2;
		%end
		x = x_;
		iter      = iter+1
		resnorm(iter) = res_;
		if (p<=p_min)
			%warning('no convergence');
			%break;
		end
		%p = 2*p;
		if (iter > maxiter)
			break;
		end
		x = obj.mesh.smooth2(x,1,0.05);

		%[f0,A] = gradpde2d(fun,x,n);
		% sum(f0)
		%x = x - A \ f0;
		%A = A + 1e-6*speye(size(A));
		%A = diag(diag(A));
		%p = 1e-2;
		%x = x + p*minres(A,f0,[],prod(n));

	
	end % while
	
	figure()
	%iter)
	clf
	obj.plot(x);
	obj.h = x;
end % infer_bed_level2

