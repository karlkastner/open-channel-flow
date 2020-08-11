% Fri 16 Mar 11:18:45 CET 2018
%% surface plot
function h = plot(obj,varargin)
%	if (isvector(val))
%		val = reshape(val,obj.mesh.n);
%	end
	obj.mesh.plot(varargin{:});
	%surface(obj.mesh.X,obj.mesh.Y,val,'edgecolor','none');
	axis tight
	axis equal
	xlabel('x')
	ylabel('y');
end

