% Sat  1 Sep 11:13:38 CEST 2018
% Karl Kästner, Berlin
%%
%% wrapper to store precomputed streamlines of potential flows
%%
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
classdef Lateral_Diversion_Wide_Channel_Map < handle
	properties
		map
		filename = '';
		change = false;
		%erel = 2.5e-4;
		erel = 1e-3;
		sfilename = 'streamlines.mat';
	end
	methods
		function obj=Lateral_Diversion_Wide_Channel_Map()
			if (nargin()>0)
				obj.sfilename = sfilename;
			end
			obj.map = containers.Map();
		end

		function streamline_dimensional()
			% TODO when ws is scaled / set, xy0 must be scaled as well
			ws = 1;
			%h  = 1;
			u0 = 1;
			if (nargin() < 4)
				beta = 1;
				key    = [shape,' ',num2str(alpha)];
				ufield = 'u';
				vfield = 'v';
			else
				% when beta is given, determine streamline along bed
				key    = [shape,' ',num2str([alpha,beta])];
				ufield = 'ubed';
				vfield = 'vbed';
			end
			fs = 2/0.41^2;
			h  = beta*ws/fs
			Qs = alpha*u0*ws*h;
				%ws = h*fs/beta;
				%h = (0.41^2/2)*Qs/u0*beta;
				%h_div_ws = (0.41^2/2)*alpha*beta;

		end


		function obj = load(obj)
			what_ = what(class(obj));
			path_ = [what_.path,filesep(),obj.sfilename];

			if (exist(path_,'file'))
				load_ = load(path_,'map');
				obj.map = load_.map;
			else
				obj.map = containers.Map();
			end
			obj.change = false;
		end % load()
		function obj = save(obj)
			if (obj.change)
				what_ = what(class(obj))
				path_ = [what_.path,filesep(),obj.sfilename]
				map   = obj.map;
				save(path_,'map');
				obj.change = false;
			end
		end % save()
	end % methods()
end % Potential_Flow_Map

