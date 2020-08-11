% Thu 30 Aug 18:12:46 CEST 2018
% 
%% analytical solutions to various depth-averaged potential flow problems
%%
% notation : small    s,n : flow
%	     captital S,N : mesh
%
% TODO support quasi unstructured meshes
% TODO make phi protected (delete automatic variables)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
classdef Potential_Flow_Analytic < Potential_Flow
	properties
		fun
	end
	properties (Access = protected)
	end % methods
	methods
		% constructor
		function obj = Potential_Flow_Analytic()
			obj.mesh = SMesh();
		end % constructor

		function ubed = ubed(obj,varargin)
			ubed = obj.fun.ubed(varargin{:});
		end	
		function ubed = vbed(obj,varargin)
			ubed = obj.fun.vbed(varargin{:});
		end
%		function ubed = ubed(obj,varargin)
%			[ubed,vbed] = obj.velocity_near_bed(varargin{:});
%		end

%		function vbed = vbed(obj,varargin)
%			[ubed,vbed] = obj.velocity_near_bed(varargin{:});
%		end

		% TODO values should be stored in mesh
		function [h] = h(obj,x,y,h)
			if (nargin()>3)
				obj.h_ = h;
			else
				h = obj.h_;
			%	if (isscalar(h))
			%		obj.h_ = h*ones(prod(obj.mesh.n),1);
			%	else
			%		obj.h_ = flat(h); 
			%	end
			end
		end

		% velocity potential
		function phi = phi(obj,x,y)
			if (nargin() < 2 || isempty(x))
				x = obj.mesh.X;
			end
			if (nargin() < 3 || isempty(y))
				y = obj.mesh.Y;
			end
			phi = obj.fun.phi(x,y);
		end

		% velocity in cartesian coordinates
		function u = u(obj,x,y)
			if (nargin() < 2 || isempty(x))
				x = obj.mesh.X;
			end
			if (nargin() < 3||isempty(y))
				y = obj.mesh.Y;
			end
			u = obj.fun.u(x,y);
		end

		function [v] = v(obj,x,y)
			if (nargin() < 2 || isempty(x))
				x = obj.mesh.X;
			end
			if (nargin() < 3 || isempty(y))
				y = obj.mesh.Y;
			end
			v = obj.fun.v(x,y);
		end



		function R = streamline_radius_of_curvature(obj,x,y)
			if (nargin() < 2 || isempty(x))
				x = obj.mesh.X;
			end
			if (nargin() < 3 || isempty(y))
				y = obj.mesh.Y;
			end
			R = obj.fun.R(x,y);
		end

		function  [dfx_dx, dfx_dy, dfy_dx, dfy_dy] = grad(obj,selector)
			switch (selector)
			case {'u'}
				dfx_dx = -obj.mesh.D2x*obj.phi_;
				dfx_dy = -obj.mesh.Dxy*obj.phi_;
				dfy_dx =  dfx_dy;
				dfy_dy = -obj.mesh.D2y*obj.phi_;
			case {'q'}
				dfx_dx = obj.mesh.Dx*obj.qx;
				dfx_dy = obj.mesh.Dy*obj.qx;
				dfy_dx = obj.mesh.Dx*obj.qy;
				dfy_dy = obj.mesh.Dy*obj.qy;
			otherwise
				error('here')
			end
		end


		function [dfs_dx, dfs_dy] = grad_mag(obj,selector)
			[dfx_dx, dfx_dy, dfy_dx, dfy_dy] = obj.grad(selector);
			fmag     = obj.mag(selector);
			switch (selector)
			case {'u'}
				fx     = obj.u;
				fy     = obj.v;
			case {'q'}
				fx     = obj.qx;
				fy     = obj.qy;
			otherwise
				error('here');
			end
			dfs_dx = 1./fmag.*(fx.*dfx_dx + fy.*dfy_dx);
			dfs_dy = 1./fmag.*(fx.*dfx_dy + fy.*dfy_dy);
		end

		% flow direction = streamline tangent
		function [ds_dx, ds_dy, dn_dx, dn_dy] = dir(obj,selector)
			mag  = obj.mag(selector);
			switch (selector)
			case {'u'}
				ds_dx = obj.u./mag;
				ds_dy = obj.v./mag;
			case {'q'}
				ds_dx = obj.qx./mag;
				ds_dy = obj.qy./mag;
			otherwise
				error('here');
			end

			if (nargout()>2)
				% orthogonal = streamline perpendicular
				dn_dx = -ds_dy;
				dn_dy =  ds_dx;
			end
		end
		
		% rotate velocities from streamline sn coordinates to euclidean
		% xy coordinates
		function [df_dx, df_dy] = sn2xy(obj,df_ds,df_dn,selector)
			[ds_dx, ds_dy, dn_dx, dn_dy] = obj.dir(selector);
			df_dx = [ds_dx.*df_ds + dn_dx.*df_dn];
			df_dy = [ds_dy.*df_ds + dn_dy.*df_dn];
		end

		function [df_ds, df_dn] = xy2sn(obj,df_dx,df_dy,selector)
			[ds_dx, ds_dy, dn_dx, dn_dy] = obj.dir(selector);
			df_ds =  dn_dy.*df_dx - dn_dx.*df_dy;
			df_dn = -ds_dy.*df_dx + ds_dx.*df_dy;
		end
	end % methods
end % classdef Potential_Flow_Analytic

