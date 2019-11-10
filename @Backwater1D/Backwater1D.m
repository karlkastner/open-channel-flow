% Thu 19 Apr 10:13:40 CEST 2018
% Karl Kästner, Berlin
%
%% solve the gradually varied flow equation (backwater equation)
%% in one dimension
%%
%% c.f. Chow, Bresse
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
classdef Backwater1D < handle
	properties
		Q0
		Qt_
		chezy_
		width_
		zb_
		% domain length
		X
		widechannel = true;
		solver = @ode23s;
		sopt = struct();
		% momentum coefficient
 		beta = 1;
		%struct('InitialStep',0.1, ...
		%	      'MaxStep', 1e-3);
		% acceleration by gravity
		g    = Constant.gravity;
		rt
	end % properties

	methods
		function obj = Backwater1D(varargin)
	                for idx=1:2:length(varargin)
				switch(varargin{idx})
				case {'issym'}
					obj.issym = varargin{idx+1};
					if (obj.issym)
						syms g positive
						obj.g = g;
					end
				otherwise
	                            obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
				end
	                end %for idx
			% object has to be created here,
			% as otherwise matlab does not allow for parallel objects
			% object members have to be created in constructors,
			% otherwise side effects occur
			% does not work, because of recursion
			% obj.rt = River_Tide();
			if (isempty(obj.rt))
				obj.rt = River_Tide('backwater',obj);
			end
		end

		function Qt = Qt(obj,x)
			if (isa(obj.Qt_,'function_handle'))
				Qt = obj.Qt_(x);
			else
				Qt = obj.Qt_;
			end
		end

		% TODO, this needs only to be tested once and then set
		function chezy = chezy(obj,x,h)
			if(isa(obj.chezy_,'function_handle'))
				switch (nargin(obj.chezy_))
				    case {1}
					chezy = obj.chezy_(x);
				    case {2}
					chezy  = obj.chezy_(x,h);
				    otherwise
					error('Chezy coefficient function must take 1 (x) or 2 (x and h) arguments.');
				end
			else
				% constant value
				chezy = obj.chezy_;
			end
		end % chezy
			
		function width = width(obj,x)
			if (isa(obj.width_,'function_handle'))
				width = obj.width_(x);
			else
				% constant value
				width = obj.width_;
			end
		end % width
		
		function dw_dx = dw_dx(obj,x)
			if(isa(obj.width_,'function_handle'))
				if (abs(nargout(obj.width_))<2)
					L = (obj.X(2)-obj.X(1));
					% TODO no magic numbers
					n = 1024;
					dx = L/n;
					x1 = x-dx;
					x2 = x+dx;
					% TODO use grad
					dw_dx = (obj.width_(x2)-obj.width_(x1))./(x2-x1);
				else
					[~, dw_dx] = obj.width_(x);
				end
			else
				% constant value
				dw_dx = obj.width_(2);
			end
		end % dw_dx

		function zb = zb(obj,x)
			if(isa(obj.zb_,'function_handle'))
				zb = obj.zb_(x);
			else
				% constant value
				zb = obj.zb_;
			end
		end

		function dzb_dx = dzb_dx(obj,x)
			if(isa(obj.zb_,'function_handle'))
				if (abs(nargout(obj.zb_))<2)
					L = (obj.X(2)-obj.X(1));
					% TODO no magic numbers
					n = 1024;
					dx = L/n;
					x1 = x-dx;
					x2 = x+dx;
					% TODO use grad
					dzb_dx = (obj.zb_(x2)-obj.zb_(x1))./(x2-x1);
				else
					[~, dzb_dx] = obj.zb_(x);
				end
			else
				% constant value
				db_dx = obj.zb_(2);
			end
		end % dzb_dx

		% normalized to domain length,
		% to increase numerical stability for long domains
		function dh_dxi = dh_dxi(obj,xi,h)
			L = (obj.X(2)-obj.X(1));
			x = obj.X(1)+xi*L;
			dh_dxi = L*obj.dh_dx(x,h);
		end
	end % properties
end % class Backwater1D

