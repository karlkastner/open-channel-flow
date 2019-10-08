% Tue 15 Aug 09:11:46 CEST 2017
%
%% velocity profile
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
classdef VP < handle
	properties
		% scalar double>0, effective sample size
		neff
		% scalar bolean, if true, acf is normalized by 1/n instead of 1/(n-k)
		acfbias
		% scalar double>0, distance between samples along cross section
		dw
		% scalar double>0, width of cross section
		width
		width_predict
		% indeces of inner region
		% [1xni] integer
		% [1xnn] bolean
		inner
		range
		% transverse profile
		t = struct(  'order', 0 ...
				);
		% vertical profile
		v = struct(  'order', 0 ...
			   , 'mode', 'profile' ...
				);
		% joint profile
		tv = struct();

		nf  = struct('pre',[0 0],'post',0);
	end % properties
	methods
		% constructor
		function obj = VP(varargin)
			obj = setfields(obj,varargin{:});
		end % constructor VP
	end % methods 
end % class VP

