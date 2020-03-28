% Sun  5 Nov 16:37:03 CET 2017
% Mi 3. Feb 12:46:57 CET 2016
%
%% Class to solve the (cross sectionally averaged) shallow water equation
%% (st venant equation)
%
% shape, amplitude
% dissipation depends on t-step
% height of wave influences stability
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
classdef SWE < handle
	properties (Constant)
		BCVER = 2;

		% 1d, 2 unknowns
		m = 2;
	end
	properties
		g = 9.81;
		% switch betwee h-q and A-Q as unknowns
		QAflag = true;
		cd
		zb
		w
		cdfun = [];
		zbfun = [];
		wfun  = [];
		sourcefun
	end
	methods
		% what is the difference?
		J               = jacobian(q);
		[Lambda R Rinv] = fluxmateig(t, q, w_hat);
		dt              = dt_cfl(obj,q,dx);

		function obj = SWE(varargin)
			for idx=1:2:length(varargin)-1
				field = varargin{idx};
				val   = varargin{idx+1};
				obj.(field) = val;
			end	
		end
	
		function d = depth(obj,Y)
			nx = size(Y,1)/2;
			if (obj.QAflag)
				d = bsxfun(@times,Y(1:nx,:),1./obj.w);
			else
				d = Y(1:nx,:);
			end
		end

		function ret = area(obj,Y)
			nx = size(Y,1)/2;
			if (obj.QAflag)
				ret = Y(1:nx,:);
			else
				ret = bxfun(@times,Y(1:nx,:),obj.w);
			end
		end

		function d = velocity(obj,Y)
			nx = size(Y,1)/2;
			% if QAflag = this is Q/A, else q/h
			d = bsxfun(@times,Y(nx+1:2*nx,:),1./Y(1:nx,:));
		end

		function ret = discharge(obj,Y)
			nx = size(Y,1)/2;
			if (obj.QAflag)
				ret = Y(nx+1:2*nx,:);
			else
				ret = bxfun(@times,Y(nx+1:2*nx,:),obj.w);
			end
		end
		
		function obj = init(obj,x)
			obj.zb = obj.zbfun(x);
			obj.w  = obj.wfun(x);
			obj.cd = obj.cdfun(x);

			% ghost point fix
			% same as k-th
			k = 3;

			field_C = {'zb','w','cd'};
			for idx=1:length(field_C)
				y = obj.(field_C{idx});
	
				if (0)
					y(end-k+1:end)   = y(end-k);
					y(1:k)           = y(k+1);
				else
					% same as first
					y(2:2+k-1)     = y(1);
					y(end-k:end-1) = y(end);
				end
				obj.(field_C{idx}) = y;

			end

			obj.sourcefun   = @(t,x,q,w) obj.source_bed_level(t,x,q,w) ...
			                           + obj.source_width(t,x,q,w) ...
			                           + obj.source_friction(t,x,q,w);
		end % init
	end % methods non-static
	methods (Static)
		% TODO, this function is not swe specific
		u_dot           = dot(t, u, flux, bc);

		% TODO, no h0 needed
		[E Epot Ekin]   = energy(H,h0,V);
		f               = flux_lin(t, q, J0);

		H               = solve_analytic(T, X, h0, H0);
		[x A Q]         = solve_stationary(L,n,bcfun,zbfun,cd);

%		[A,B,rhs] = bc_inflow(t, q0, dt, q1fun);
%		[A,B,rhs] = bc_inflow_low_pass(t, q1, dt, q0fun, Tlp, cd, w0, dzb_dx);
%		[A,B,rhs] = bc_discharge_non_reflecting(t, q1, dt, q0fun);
%		[A,B,rhs] = bc_level(t, q, dt, hfun);
%		[A,B,rhs] = bc_level_sommerfeld(t, q, dt, hfun);
%		[A,B,rhs] = bc_incoming_non_reflecting(t, q, dt, hfun, w0);
%		[A,B,rhs] = bc_nonreflecting(t, q, dt);
%		[A,B,rhs] = bc_reflecting(t, q, dt);

	end % methods (Static)
	methods % non-static
	end
end % classdef

