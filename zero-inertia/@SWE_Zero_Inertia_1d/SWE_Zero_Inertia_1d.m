% Tue  8 Apr 11:20:35 CEST 2025
% Karl Kastner, Berlin
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
classdef SWE_Zero_Inertia_1d < RAD_Model
	properties
		g = 9.81;
		% TODO move to opt
		input_factor  = 1;
		output_factor = 1;
%		eps = 1e-4;
		infiltration_rate_linear     = 0;
		infiltration_rate_constant   = 0;
		precipitation_rate = @(t,dt) 0;
		% ground level
		zb;
	end
	methods
		function obj = SWE_Zero_Inertia_1d(varargin)
			obj = obj@RAD_Model();
			obj.opt.output.base_str = 'zi1d-';

			% linear resistance coefficient
			obj.p.cd1 = [];
			% quadratic resistance coefficient
			obj.p.cd2 = [];
% TODO make parameter
%			obj.p.infiltration_rate_linear     = 0;
%			obj.p.infiltration_rate_constant   = 0;
%			obj.p.precipitation_rate = @(t,dt) 0;
%			% ground level
%			obj.p.zb;
			obj.pmu.cd1 = 0;
			obj.pmu.cd2 = 0;
			obj.pss.cd1 = 0;
			obj.pss.cd2 = 0;
			obj.psl.cd1 = 0;
			obj.psl.cd2 = 0;
			obj.psdist.cd1 = 'normal';
			obj.psdist.cd2 = 'normal';

			obj.nvar = 1;
			obj.opt.analytic_jacobian = true;
			obj.opt.nlsolver.abstol   = sqrt(eps);
			obj.opt.nlsolver.maxiter  = 10;
			obj.opt.spatial_discretization = 'optimal';
			obj.opt.linear_solver.name = 'mldivide';
			obj.opt.linear_solver.maxiter = 10;
			obj.opt.linear_solver.tol = sqrt(eps);
			obj.opt.returnmat = true;
			obj.opt.heps      = sqrt(eps);
			obj.opt.time_integration.scheme = 'trapezoidal';
			obj.opt.twosided  = false;
			obj.opt.output.store_interface_values = false;
			obj.opt.Seps = 1e-2*sqrt(eps);
			obj.opt.step_split_source = 0;
			obj.hashfield_C{end+1} = 'infiltration_rate_linear';
			obj.hashfield_C{end+1} = 'infiltration_rate_constant';
			obj.hashfield_C{end+1} = 'zb';
			obj.hashfield_C{end+1} = 'opt.precipitation';
			obj.boundary_condition = {[0,0,1,0],[0,0,1,0]};
			obj.opt.analytic_flux_derivatives = true;
			obj.opt.Jeps = sqrt(eps);
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end % constructor
	end % methods
end % SWE_Zero_Inertia_1d

