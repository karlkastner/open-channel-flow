'here'
AA = A(1:ns,1:ns);
cond(AA)
AA = A(ns+1:2*ns,1:ns);
cond(AA)
AA = A(1:ns,ns+1:end);
cond(AA)
AA = A(ns+1:end,ns+1:end);
cond(AA)
if (1)
		% enforce smoothmess with thikhonov regularization
		%D = spdiags(ones(ns,1)*[-1/2,0,1/2],-1:1,ns,ns);
		D = spdiags(ones(ns,1)*[-1,1,0],-1:1,ns,ns);
		D(1,1:2) = [-1,1];
		D(end,end-1:end) = [-1,1];
		D = 1./ns*D;
		DD = zeros(ns+ns_,ns+ns_);
		DD(1:ns,1:ns) = D'*D;
		lambda = 0*1e2;
		c = (A'*A + lambda*DD) \ (A'*b);
end

