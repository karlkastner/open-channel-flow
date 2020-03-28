% Sun 20 May 13:13:15 CEST 2018
%% apply boundary condition for stationary flow
function [A,b,obj] = apply_boundary_stationary(obj,A,b,h,qx,qy)
	nn = prod(obj.mesh.n);	

	% get all points on the boundary
	% and next inner point
	[id, jd] = obj.mesh.boundary_indices();
	nb = length(id);
	% TODO this is more complicated at the corners,
	% there derivatives with respect to two neighbours have to be
	% considered, if the mesh is curved, the diagonal point as well,
	% or when the conditions change along the boundary
	% singly removed points effect all 8 neighbours
	% -> no the 8 points are boundary points
	
	% get the conditions at the boundary
	[p, val] = obj.bcfun(obj.mesh.X(id),obj.mesh.Y(id),h(id),qx(id),qy(id));
	
	% apply them for each of the three variables h, qx, qy
	for idx=1:3
		val_ = val((1:nb)+nb*(idx-1));
		p_   = p((1:nb)+nb*(idx-1));
		% only set for those, where val is not NaN
		fdx = ~isnan(val_);
		id_ = nn*(idx-1) + id(fdx);
		if (~isempty(id))
			jd_ = nn*(idx-1) + jd(fdx);
			p_ = p_(fdx);
			dx  = hypot(obj.mesh.X(id(fdx))-obj.mesh.X(jd(fdx)),obj.mesh.Y(id(fdx))-obj.mesh.Y(jd(fdx)));
			A(sub2ind(3*[nn,nn],id_,id_)) = p_ - (1-p_)./dx;
			A(sub2ind(3*[nn,nn],id_,jd_)) = (1-p_)./dx;
			b(id_)                 = val_(fdx);
		end % if
	end % for
	
end % apply_boundary_condition
	
