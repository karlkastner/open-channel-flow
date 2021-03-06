% Sun  9 Feb 11:57:29 +08 2020
%
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
classdef Lateral_Diversion_Finite_Width_Gradual < Lateral_Diversion_Finite_Width
	properties
		% number of segments across outlet
		m  = 100;
		% coefficients for potentials and fourier expansion in the side branch
		cp
		cs
		cflag = 1;
		mode = 0;
		lmode = 0;
	end
	methods (Static)
		f = derive();
	end
	methods
		function obj = Lateral_Diversion_Finite_Width_Gradual(varargin)
			obj                 = obj@Lateral_Diversion_Finite_Width();
			%obj.funfilename     = 'functions.mat';
			for idx=1:2:length(varargin)
				field       = varargin{idx};
				obj.(field) = varargin{idx+1};
			end
			% TODO, stupid hack to avoid recursion
			%if (~issym(obj.alpha))
			%	obj.load_functions();
			%end
		end

		% Qin in non-normalized quantities
		function Qin = Qin(obj)
			Qin = obj.uin*(obj.ws/obj.gamma)*obj.h;
		end

		% Qs in natural units 
		%function Qs = Qs(obj)
		%	Qs = obj.uin*obj.alpha.*obj.ws*obj.h;
		%end

		% Q0 in non-normalized quantities
		% note that his is not Qin
		function Q0 = Q0(obj)
			Q0 = obj.uin*obj.w0*obj.h;
		end

		% main channel width in non-normalized quantities
		function w0 = w0(obj)
			w0  = obj.ws./obj.gamma;
		end

		% 
		% all following parameters are normalized
		%

		% mean velocity
		function u0 = u0(obj)
			u0 = 1 - 1/2*obj.alpha.*obj.gamma;
		end

		% initial coordinate of dividing streamline of depth averaged flow
		function ys_in = ys_in(obj)
			ys_in = obj.alpha;
		end

		% upstream point, where error is smaller than tol
		function xs_in = xs_in(obj,tol)
			xs_in  = obj.fun.xin(obj.alpha,tol);
		end

		function u  = u(obj,x,y)
			u   = obj.u0 + obj.evalk('u',x,y);
		end

		function v = v(obj,x,y)
			v = obj.evalk('v',x,y);
		end
		function du_dx = du_dx(obj,x,y)
			du_dx = obj.evalk('du_dx',x,y);
		end
		function du_dy = du_dy(obj,x,y)
			du_dy = obj.evalk('du_dy',x,y);
		end
		function dv_dx = dv_dx(obj,x,y)
			dv_dx = obj.evalk('dv_dx',x,y);
		end
		function dv_dy = dv_dy(obj,x,y)
			dv_dy = obj.evalk('dv_dy',x,y);
		end

		function d2u_dx2 = d2u_dx2(obj,x,y)
			d2u_dx2 = obj.evalk('d2u_dx2',x,y);
		end	
		function d2u_dxdy = d2u_dxdy(obj,x,y)
			d2u_dxdy = obj.evalk('d2u_dxdy',x,y);
		end	
		function d2u_dy2 = d2u_dy2(obj,x,y)
			d2u_dy2 = obj.evalk('d2u_dy2',x,y);
		end	

		function d2v_dx2 = d2v_dx2(obj,x,y)
			d2v_dx2 = obj.evalk('d2v_dx2',x,y);
		end	
		function d2v_dxdy = d2v_dxdy(obj,x,y)
			d2v_dxdy = obj.evalk('d2v_dxdy',x,y);
		end	
		function d2v_dy2 = d2v_dy2(obj,x,y)
			d2v_dy2 = obj.evalk('d2v_dy2',x,y);
		end	

		function [u,v] = uv(obj,x,y)
		switch (obj.mode)
		case {0}
			u   = obj.u(x,y); 
			v   = obj.v(x,y);
		case {1}
			[u,v] = obj.uv1(x,y);
		end % switch
		end

		function [ub,vb] = uvb(obj,x,y)
			u        = obj.u(x,y);
			v        = obj.v(x,y);
			R        = obj.R(x,y,u,v);
			% note that this is indeed beta, not h in the relation
			[ub,vb]  = bend_velocity_near_bed(u,v,obj.beta/obj.fs,R);
			if (0)
			uf       = obj.u_far(x,y);
			ubf      = obj.ubf(x,y);
			vbf      = obj.vbf(x,y);
			fdx      = (    ((x<-0.5) & (u<uf)) ...
			             | ((x>+0.5) & (u>uf)) );
			ub(fdx) = ubf(fdx);
			vb(fdx) = vbf(fdx);
			end
		end

		function ub = ubed(obj,x,y)
			u = obj.u(x,y);
			v = obj.v(x,y);
			R = obj.R(x,y,u,v);
			[ub,vb] = bend_velocity_near_bed(u,v,obj.beta/obj.fs,R);
		end

		function vb = vbed(obj,x,y)
			u = obj.u(x,y);
			v = obj.v(x,y);
			R = obj.R(x,y,u,v);
			[ub,vb] = bend_velocity_near_bed(u,v,obj.beta/obj.fs,R);
		end

		function ubf = ubf(obj,x,y)
			ubf = obj.fun.ubf(x,y,obj.alpha,obj.beta,obj.gamma);
		end

		function vbf = vbf(obj,x,y)
			vbf = obj.fun.vbf(x,y,obj.alpha,obj.beta,obj.gamma);
		end

		function R = R(obj,x,y,u,v,du_dx,du_dy,dv_dx,dv_dy)
			if (nargin()<6)
			if (nargin()<4)
				u = obj.u(x,y);
				v = obj.v(x,y);
			end
			du_dx = obj.du_dx(x,y);
			du_dy = obj.du_dy(x,y);
			dv_dx = obj.dv_dx(x,y);
			dv_dy = obj.dv_dy(x,y);
			end
			R = -streamline_radius_of_curvature( ...
						    u,du_dx,du_dy, ...
						    v,dv_dx,dv_dy);
		end
		function J = J(obj,x,y)
			du_dx = obj.du_dx(x,y);
			du_dy = obj.du_dy(x,y);
			dv_dx = obj.dv_dx(x,y);
			dv_dy = obj.dv_dy(x,y);
			J = [diag(sparse(du_dx)),diag(sparse(du_dy));
			     diag(sparse(dv_dx)),diag(sparse(dv_dy))];
		end
	end % methods
end % Lateral_Diversion_Finite_Width_Gradual

