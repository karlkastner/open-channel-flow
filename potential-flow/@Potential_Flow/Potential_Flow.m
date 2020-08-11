% Wed 20 May 11:34:33 +08 2020
%
%% numerical potential flow solver by various methods
%% analytic (series)
%% or FDM on non-orthogonal curvilinear grids
%% or FEM by unstructured meshes
%
% notation : small    s,n : flow
%	     captital S,N : mesh
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
classdef Potential_Flow < handle
	properties
		% mesh/grid of computational domain
		% either structured or unstructured
		mesh = [];
		% potential flow discretisation matrix
		% the bilinear operator A(phi,phi) = b(phi)
		A
		% right hand side
		% the linear operator
		b

		% path for servlet deployment
		deploypath = './';

		% if false, intepreted as streamfunction, not as potential	
		ispotential = true;

		% options for streamline computation
		%streamlineopt = struct('MaxStep',10);
		streamlineopt = struct('RelTol',1e-2, ...
				       'solver',@ode23s);

		% physical constants
		% TODO not here
		g = Constant.gravity();

		% chezy coefficient
		chezy = 60;

		% bed material grain size
		d_mm = 0.2;

		solver = 'gaussian-elimination';
	end % properties
	properties (Access = protected)
		% velocity potential
		phi_;
		% water depth
		h_
		u_
		v_
		% TODO also store Rc and umag
		ubed_
		vbed_
		R_
		interp_
	end % properties (protected)
	methods (Static)
		[mesh,fun]=test_case(id,n);
	end % methods (Static)
	methods
		function clear(obj)
			obj.phi_ = [];
			obj.h_ = [];
			obj.ubed_ = [];
			obj.vbed_ = [];
			obj.interp_ = [];
		end

		function val = interp(obj,field,x,y)
			if (~isfield(obj.interp_,field))
				% note : the value of the interpolation of scatteredInterpolant
				% can be reassigned, but this is slow
				obj.interp_.(field) = scatteredInterpolant(flat(obj.mesh.X),flat(obj.mesh.Y),obj.(field),'natural','linear');
			end
			val = obj.interp_.(field)(x,y);
		end % interp
i
		% TODO values should be stored in mesh
		function [h] = h(obj,h)
			if (nargin()<2)
				h = obj.h_;
			else
				if (isscalar(h))
					obj.h_ = h*ones(prod(obj.mesh.np),1);
				else
					obj.h_ = flat(h); 
				end
			end
		end

		% velocity potential
		function phi = phi(obj,varargin)
			if (nargin() < 2)
				phi = obj.phi_;
			else
				phi = obj.interp('phi',varargin{:});
			end
		end


		% velocity in cartesian coordinates
		function u = u(obj,varargin)
			if (nargin() < 2)
				if (obj.ispotential)
					u = -obj.mesh.Dx*obj.phi_;
				else
					u =  obj.mesh.Dy*obj.phi_;
				end
			else
				u = obj.interp('u',varargin{:});
			end
		end

		function v = v(obj,varargin)
			if (nargin() < 2)
				if (obj.ispotential)
					v = -obj.mesh.Dy*obj.phi_;
				else
					v = -obj.mesh.Dx*obj.phi_;
				end
			else
				v = obj.interp('v',varargin{:});
			end
		end

%		function umag = umag(obj,varargin)
%			umag = hypot(obj.u(varargin{:}), obj.v(varargin{:}));
%		end

		% velocity magnitude = velocity along streamline
		function mag = mag(obj,selector,varargin)
			switch (selector)
			case {'u'}
				mag = hypot(obj.u(varargin{:}),obj.v(varargin{:}));
			case {'q'}
				mag = hypot(obj.qx,obj.qy);
			otherwise
				error('here');
			end
		end

		% flow direction = streamline tangent
		function [ds, dn] = dir(obj,selector)
			mag  = obj.mag(selector);
			switch (selector)
			case {'u'}
				ds = [obj.u, obj.v]./mag;
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

		function [dfx_dx, dfx_dy, dfy_dx, dfy_dy] = grad(obj,selector,varargin)
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


		function ubed = ubed(obj,varargin)
			if (isempty(obj.ubed_))
				[obj.ubed_,obj.vbed_] = obj.velocity_near_bed();
			end
			if (nargin() < 2)
				ubed = obj.ubed_;
			else
				ubed = obj.interp('ubed',varargin{:});
			end
		end

		function vbed = vbed(obj,varargin)
			if (isempty(obj.vbed_))
				[obj.ubed_,obj.vbed_] = obj.velocity_near_bed();
			end
			if (nargin() < 2)
				vbed = obj.vbed_;
			else
				vbed = obj.interp('vbed',varargin{:});
			end
		end

		% TODO allow for bed
		function R = streamline_radius_of_curvature(obj,varargin)
			if (isempty(obj.R_))
				u = obj.u(varargin{:});
				v = obj.v(varargin{:});
				[du_dx, du_dy, dv_dx, dv_dy] = obj.grad('u',varargin{:});
				R = streamline_radius_of_curvature(u,du_dx,du_dy,v,dv_dx,dv_dy);
			else
				R = obj.interp('R',varargin{:});
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
end % class Potential_Flow

