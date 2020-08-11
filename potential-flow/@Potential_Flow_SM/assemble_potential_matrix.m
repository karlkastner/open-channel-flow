% Mon 12 Mar 11:54:08 CET 2018
% Karl Kastner, Berlin
%
%% assemble the discretisation matrix for potential flow
function [bcmark, obj] = assemble_potential_matrix(obj,cutfun)

	% not necessary to reassemble
	% [Dx, Dy, La] = 
	% obj.mesh.derivative_matrices();

	if (nargin() > 1 && ~isempty(cutfun))
		% cut from the difference matrices
		obj.mesh.cut_from_domain(cutfun);
	end

	% laplacian in the interior of the domain
	% to avoid spurious oscillations following identitiy is used
	% applyed
	% D(h D x) = h D2 x + Dx Dh 
	if (0)
		H = diag(sparse(obj.h));
		A =   obj.mesh.Dx*(H*obj.mesh.Dx) ...
		    + obj.mesh.Dy*(H*obj.mesh.Dy);
	else
		H = obj.h;
		if (~isempty(H))
		A =   diag(sparse(obj.mesh.Dx*H))*obj.mesh.Dx ...
		    + diag(sparse(H))*obj.mesh.D2x ...
		    + diag(sparse(obj.mesh.Dy*H))*obj.mesh.Dy ...
		    + diag(sparse(H))*obj.mesh.D2y;
		else
		A = obj.mesh.D2x + obj.mesh.D2y;
		end
	end

	% right hand side
	nn = obj.mesh.n;
	b  = zeros(prod(nn),1);

	% boundary condition
	[A,b,bcmark] = obj.mesh.apply_boundary_condition(A,b,obj.h);
	
	obj.A  = A;
	obj.b  = b;
end % assemble_potential_matrix

