% Sat 19 May 19:33:16 CEST 2018
% TODO move to mesh
function quiver(obj,u,v)
	x = flat(obj.mesh.X);
	y = flat(obj.mesh.Y);
	quiver(x,y,u,v);
end
