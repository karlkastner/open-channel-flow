% 2018-03-15 10:08:09.889670709 +0100
%% assemble the discretisation matrix
function obj = assemble_discretisation_matrix_rectilinear(A)
		h = L./(n-1);
		x1 = L(1)*(0:n(1)-1)/(n(1)-1);
		x2 = L(2)*(0:n(2)-1)/(n(2)-2);
	
		% laplacian for the interior (n1-2)x(n2-2)
		for idx=1:2
			%I{idx} = speye(n(idx));
			I{idx} = diag(sparse([0;ones(n(idx)-2,1);0]));
			A{idx} = 1/h(idx)^2*spdiags(ones(n(idx),1)*[1,-2,1],-1:1,n(idx),n(idx));
			% not at boundary
			A{idx}(1,:)   = 0;
			A{idx}(end,:) = 0;
			b{idx} = zeros(n(idx),1);
		end
		% set up 2d-matrix
		AA = ( kron(A{1},I{2}) + kron(I{1},A{2}) );

		% TODO use derivative_matrix_2
		Dx = derivative_matrix_1_1d(n(1),L(1));
	%	D2x = derivative_matrix_2_1d(n(1),L(1));
		Dx = kron(Dx,speye(n(2)));
	%	D2x = kron(D2x,speye(n(2)));
		Dy = derivative_matrix_1_1d(n(2),L(2));
	%	D2y = derivative_matrix_2_1d(n(1),L(1));
		Dy = kron(speye(n(1)),Dy);
	%	D2y = kron(speye(n(1)),D2y);
end

