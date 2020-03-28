% Sat  4 Nov 19:45:25 CET 2017
%
%% stationary solution to the SWE
%
% TODO, make the solver more general and move to finite volume
function [x A Q] = swe_stationary(L,n,bcfun,zbfun,cd)
%	n=10
	reltol = 1/n;

	dx = L/(n-1);
	D  = derivative_matrix_1_1d(n,L);
	Z = sparse(n,n);
	g = 9.81;	

	x  = dx*(0:n-1)';
	zb = zbfun(x);
	dzbdx = cdiff(zb)/dx;

	rhs = zeros(2*n,1);
	A = rand(n,1);
	Q = rand(n,1);
	y =[A;Q];
	Amin = sqrt(eps);

	p = 1;

	M = 1/4*spdiags(ones(n,1)*[1,2,1],-1:1,n,n);
	M(1,1:2) = [1,1]/2;
	M(n,n-1:n) = [1,1]/2;

	for idx=1:2*n
		Aold = A;	
		% TODO, use matrix functions here
		Lhs = [Z,D;
	               (diag(sparse(-Q.^2./(A).^2)) + g*diag(sparse((A))))*D + g*diag(sparse(dzbdx)), ...
		       diag(sparse(2*Q./(A)))*D + cd*diag(sparse(abs(Q)./(A).^2))];

		% evaluate bc
		% TODO t0 not 0
		[lA lrhs]    = feval(bcfun{1},0);
		Lhs(1,:)     = 0;
		Lhs(1,1)     = lA(1,1) - lA(1,2)/dx;
		Lhs(1,2)     = lA(1,2)/dx;
		rhs(1)       = lrhs(1);
		Lhs(n+1,:)   = 0;
		Lhs(n+1,n+1) = lA(2,1) - lA(2,2)/dx;
		Lhs(n+1,n+2) = lA(2,2)/dx;
		rhs(n+1)     = lrhs(2);

		[rA rrhs]    = feval(bcfun{2},0);
		Lhs(n,:)     = 0;
		Lhs(n,n)     =  rA(1,1) + rA(1,2)/dx;
		Lhs(n,n-1)   = -rA(1,2)/dx;
		rhs(n)       = rrhs(1);
		Lhs(2*n,:)   = 0;
		Lhs(2*n,2*n)   =  rA(2,1) + rA(2,2)/dx;
		Lhs(2*n,2*n-1) = -rA(2,2)/dx;
		rhs(2*n)     = rrhs(2);

if (0)		
		Lhs(1,:)     = 0;
		Lhs(1,1)     = lA(1,1) - lA(1,3)/dx;
		Lhs(1,2)     = lA(1,3)/dx;
		Lhs(1,n+1)   = lA(1,4) - lA(1,6)/dx;
		Lhs(1,n+2)   = lA(1,6)/dx;
		rhs(1)       = lrhs(1); 	
		
		Lhs(n+1,:)   = 0;
		Lhs(n+1,1)   = lA(2,1) - lA(2,3)/dx;
		Lhs(n+1,2)   = lA(2,3)/dx;
		Lhs(n+1,n+1) = lA(2,4) - lA(2,6)/dx;
		Lhs(n+1,n+2) = lA(2,6)/dx;
		rhs(n+1)     = lrhs(2); 	
		
		Lhs(n,:)     = 0;
		Lhs(n,n)     = rA(1,1) + rA(1,3)/dx;
		Lhs(n,n-1)   = -rA(1,3)/dx;
		Lhs(n,2*n)   = rA(1,4) + rA(1,6)/dx;
		Lhs(n,2*n-1) = -rA(1,6)/dx;
		rhs(n)       = rrhs(1); 	
		
		Lhs(2*n,:)     = 0;
		Lhs(2*n,n)     =  rA(2,1) + rA(2,3)/dx;
		Lhs(2*n,n-1)   = -rA(2,3)/dx;
		Lhs(2*n,2*n)   = rA(2,4) + rA(2,6)/dx;
		Lhs(2*n,2*n-1) = -rA(2,6)/dx;
		rhs(2*n)       = rrhs(2); 	
end % if 0

		% smoothing with M, to asure stability
		y  = (1-p)*y + p*[M,Z;Z,M]*(Lhs\rhs);

		A = y(1:n);
		Q = y(n+1:2*n);

		delta       = abs(A-Aold);
		if (all(delta<reltol*abs(A)))
			fprintf(1,'Stationary flow solver converged in %d iterations with delta(A)/A*nx %f \n',idx,max(delta./abs(A))*n);
			break;
		end
	end % for idx
end % swe_stationary


