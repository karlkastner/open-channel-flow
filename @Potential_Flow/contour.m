% Fri 18 May 11:42:11 CEST 2018
%% contour plot of the potential flow solution
function h = contour(obj,val,varargin)
	if (isvector(val))
		val = reshape(val,obj.mesh.n);
	end
	% necessary as a workaround for bug: Warning: Contour not rendered for constant ZData
	ax = gca;
	%surface(obj.mesh.X,obj.mesh.Y,val,'edgecolor','none');
	%contour(obj.mesh.X,obj.mesh.Y,val,varargin{:});
	contourf(obj.mesh.X,obj.mesh.Y,val,varargin{:});
	axis tight;
	axis equal;
	xlabel('x');
	ylabel('y');
end

