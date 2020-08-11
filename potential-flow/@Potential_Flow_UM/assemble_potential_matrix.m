% Wed 13 May 20:19:41 +08 2020
function assemble_potential_matrix(obj)
	% assemble discretization matrix for weak formulation
	% int Delta u v d omega = int n u v d Gamma - int grad u grad v d Omega

	switch (obj.mode)
	case {'hermite'}
		A   = -obj.assemble_2d_dphi_dphi_hermite();
		rhs = zeros(3*obj.np+obj.nelem,1);
		%[A,rhs] = obj.apply_boundary_condition_hermite(A,rhs,xb,yb,p,val);
		[A,rhs] = obj.mesh.apply_boundary_condition(A,rhs,[],obj.mode);
	case {'lagrange'}
		A = -obj.mesh.assemble_2d_dphi_dphi_lagrange();
%		A   = -obj.assemble_2d_dphi_dphi();
		rhs = zeros(obj.mesh.np,1);
		[A,rhs] = obj.mesh.apply_boundary_condition(A,rhs,[],obj.mode);
	end
	obj.A = A;
	obj.b = rhs;
end
