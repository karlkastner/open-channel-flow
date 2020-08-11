% Sat 19 May 19:33:16 CEST 2018
function quiver(obj,u,v,varargin)
	if (nargin()<2 || isempty(u))
		u = obj.u;
		v = obj.v;
	end
%	x = flat(obj.mesh.X);
%	y = flat(obj.mesh.Y);
%	quiver(x,y,u,v,varargin{:});
	obj.mesh.quiver(u,v,varargin{:});
end
