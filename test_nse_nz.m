function test_nse_xz()
	n = [10,20];
	L = [

end % test_nse

% stationary
% d/dx terms are zero
function [u, v, w] = nse_yz(n,L)
	dy = L(1)/(n(1)-1);
	dz = L(2)/(n(2)-1);
	
	% allocate memory
	u = flat(zeros(n(1)-1,n(2)-1));
	v = flat(zeros(n(1),n(2)-1));
	w = flat(zeros(n(1)-1,n(2)));

	uvw = [u;v;w];

	while (1)

	% how to satisfy continuity?

	% Dt = 0
	% Dx = 0
	% Dp/Dy = 0
	%diag_u = diag(sparse(u));
	diag_v  = diag(sparse(v));
	diag_vc = diag(sparse(mid(v.').'));
	diag_w  = diag(sparse(w));
	diag_wc = diag(sparse(mid(w)));
	% note that the difference matrices have different number of rows for w and z
	% note that difference matrice has to centres / cs has to be for two
	% TODO choose dimension 1
	Dy_cc = difference_matrix_2d(nu,[L(1)-dx,L(2)-dy]);
	Dz_cc = difference_matrix_2d(nu,[L(1)-dx,L(2)-dz]);
	A = [ diag_vc*Dy_cc + diag_wc*Dz_cc,   		      Z,                       Z,
                                      Z, diag_v*Dy + diag_wc*Dzc,                      Z,
                                      Z,                      Z, diag_vc*Dyc + diag_w*Dz];



	% TODO mu is not const
	A = A + mu*[ D2y+D2z,       Z, Z;
                           Z, D2y+D2z, Z;
	                   Z,       Z, D2y+D2z];

	% for u, only implicit bc has to be applied
	% bottom w = 0, top  dw/dz = 0
	% left   v = 0, right v = 0

	% TODO, boundary condition

	uvw = (1-p)*uvw + p*(A \ rhs);
	
	end		   
	 

end

