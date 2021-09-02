% Wed 21 Mar 15:06:47 CET 2018
%
%% Transverse profile of the bed level and bed material grain size in
%% an equilibrium (infintely long) meander bend
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
classdef Equilibrium_Bend < handle
	properties
		phi;
		g     = Constant.gravity;
		C     = 60;
		Q
		Qs0   = 100;
		width
		Rc
		Hc
		% velocity distribution
		f_ =   1;
		D_max
		D_min
		% steps across channel
		nx      = 100;
		maxiter = 100;

		% TODO out struct
		% computed depth
		h
		% computed grain size
		D

		odeset;
		relaxation = struct('gsd',0.1);
		%u_bar = 1;
		T_C = 25;
	end % properties
	methods
		function obj = Equilibrium_Bend()
			obj.odeset.RelTol = 1e-3;
			obj.odeset.AbsTol = 1e-6;
		end % constructor

		% discretised across channel coordinate
		function r = r(obj)
			nx = obj.nx;
			Rc = obj.Rc;
			w  = obj.width;
			r = linspace(Rc-w/2,Rc+w/2,nx)';
		end
		function f = f(obj,h)
			if (~isscalar(obj.f_))
				r = obj.r();
				if (nargin() < 2)
					h = obj.h;
				end
				f = obj.f_(r,h);
				% normalize
				f  = f./mean(f);
			else
				f = ones(obj.nx,1);
			end
		end

%		function q = q(obj,varargin)

		function area = area(obj)
			area = mean(obj.h)*obj.width;
		end

		function u = u(obj,h)
			if (nargin() < 2)
				h = obj.h;
			end
			r  = obj.r;
			f  = obj.f(h); %r);
			Fr = f.*r./obj.Rc;
			qscale = mean(Fr.*h)*obj.width;
			u  = Fr.*obj.Q/qscale;
			%u  = Fr.*obj.u_bar;
		end

		function u_bar = u_bar(obj,h)
			u = obj.u(h);
			u_bar = sum(h.*u)./sum(h);
		end


%		function Q = discharge(obj,h)
%			r = obj.r;
%			q = h.*obj.u(h);
%			Q = (r(end)-r(1))*mean(h.*u);
%		end
	end % methods
end % Equilibrium_bend

