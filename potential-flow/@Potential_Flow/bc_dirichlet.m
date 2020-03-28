% Mon 12 Mar 11:54:08 CET 2018
%% apply Dirichlet boundary conditions
function obj = bc_dirichlet(obj,bc)
	A = obj.A;
	b = obj.b;
	X = flat(obj.mesh.X);
	Y = flat(obj.mesh.Y);
	m = length(X)*[1 1];
	for idx=1:length(bc)
		% get distance to boundary
		% TODO, only apply to boundary points
		% TODO, this fails for small step widths, if boundary is curved
		[xp,yp,p,d]        = Geometry.plumb_line(bc(idx).x0(1),bc(idx).y0(1),bc(idx).x0(2),bc(idx).y0(2),X,Y);
		% TODO make threshold dependent on step size
		d_min    = 0.1;
		id       = find(d<d_min);
		A(id,:) = 0;
		A(sub2ind(m,id,id)) = 1;
		b(id)   = bc(idx).f(X(id),Y(id));
	end
	obj.A = A;
	obj.b = b;
end

