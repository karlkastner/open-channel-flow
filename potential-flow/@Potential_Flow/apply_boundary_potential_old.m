	% boundary conditions
	% upper boundary and lower boundary (2n-2) vertices
	% no outflow, normal flow dphi/dn is zero
	%A                = spzeros(n(2),n(2));
	%A(1,1:2)         = 1/h(1)*[-1,1];
	%A(end,end-1:end) = 1/h(1)*[-1,1];
	%I1 = diag(sparse([0;ones(n(1)-2,1);0]));
	%AA = AA+kron(I1,A);

	m = prod(n)*[1,1];

	% TODO corners
	[h1,h2] = obj.mesh.hSN();
	% TODO quick fix
	h1(:,end+1) = h1(:,end);
	h2(end+1,:) = h2(end,:);

	% TODO, apply conditions with separate functions at each boundary

	% inflow across all boundaries is at first set to zero
	% this is later overwritten by dirichlet inflow conditions
	% as specified by the bc

	% flow across left boundary to zero (no inflow)
	id                          = (1:n(1));
	AA(id,:)                    = 0;
	AA(sub2ind(m,id,id))        = -1./h1(id);
	AA(sub2ind(m,id,id+n(1)))   =  1./h1(id);

	% flow across right boundary to zero (no inflow)
	id                          = ((n(2)-1)*n(1)+1:n(2)*n(1));
	AA(id,:)                    = 0;
	AA(sub2ind(m,id,id))        = -1./h1(id);
	AA(sub2ind(m,id,id-n(1)))   =  1./h1(id);

	% flow across upper boundary to zero (no inflow)
	id                          = (1:n(2)-2)*n(1)+1;
	AA(id,:)                    = 0;
	AA(sub2ind(m,id,id))        = -1./h2(id);
	AA(sub2ind(m,id,id+1))      = +1./h2(id);

	% flow across lower boundary to zero (no inflow)
	id                          = (2:n(2)-1)*n(1);
	%(2:n(2)-1)*n(1):n(1):(n(2)-1)*n(1);
	AA(id,:)                    = 0;
	AA(sub2ind(m,id,id))        = -1./h2(id);
	AA(sub2ind(m,id,id-1))      =  1./h2(id);
	
	% TODO tangential flow at outflow to zero ?

	% inflow boundary: set velocity potential to a prescribed value
	for idx=1:n(1)
		AA(idx,:) = 0;
		AA(idx,idx) = 1;
	end

	% outflow boundary: set velocity potential to a presecribed value
	for idx=1:n(1)
		id        = (n(2)-1)*n(1)+idx;
		AA(id,:)  = 0;
		AA(id,id) = 1;
	end

	% right hand side
	bb = kron(b{1},ones(n(2),1)) + kron(ones(n(1),1),b{2});

	obj.A = AA;
	obj.b = bb;

	obj.bc_dirichlet(varargin{:});

