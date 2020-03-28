% Sat 10 Mar 12:05:25 CET 2018
%% side weir, analytical solution to (critical) lateral outflow
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
classdef Side_Weir < handle
	properties
		g = 9.81;
		alpha = 1;
		Q
		width
		h
		method = 'exact';
		shape = 'rectangular';
	end % properties

	methods

	% sudden drop of water level at upstream end of weir
	function [delta_zs] = delta_zs(obj)
		% c.f. Paris 2012, eq 1 (typo also in 5)
		% 2*g*h1 + alpha*(Q0/A)^2 = 2*g*h2 + alpha*(Q0^2 - Qs^2)/A^2
		% dzs = -1/(2*g)*(Qs/A)^2;
		% there is a typo, it should read:
		% 2*g*h1 + alpha*(Q0/A)^2 = 2*g*h2 + alpha*(Q0 - Qs)^2/A^2
		% dzs = 1/g*(Q0*Qs-1/2*Qs^2)/A^2
		%dzs = -1/(2*g)*(Q0^2-Qs^2)/A^2;

		g  = obj.g;
		Qu = obj.Q.upstream;
		hu = obj.h.upstream;
		wu = obj.width.upstream;
		Qs = obj.Q.side;
		ws = obj.width.side;
		Au = hu*wu;

		switch (obj.method)
		case {'exact'}
			% exact integration
			delta_zs = hu/2*log( 1 + (Qu^2 - (Qu-Qs)^2)/(g*hu*Au^2 - Qu^2) );
		case {'approximal'}
			delta_zs = (Qu^2 - (Qu-Qs)^2)/(2*g*Au^2);
		case {'numerical'}
		otherwise
			error('here');
		end 
	end % delta_zs

	function derive_surface_elevation(obj)
		syms  Qs Qu ws hu wu x g z h zs
		obj.Q.side         = Qs;
		obj.Q.upstream     = Qu;
		obj.h.upstream     = h;
		obj.width.side     = ws;
		obj.width.upstream = wu;
		obj.g              = g;
		
		dzs_dx = obj.dzs_dx(x,zs);
		I      = int(dzs_dx,x);
		zs     = I - subs(I,x,0);
		zs     = simplify(zs,'ignoreanalyticconstraints',true);
		zs     = combine(zs,'log','IgnoreAnalyticConstraints',true)
		%zs = simplify(zs,'ignoreanalyticconstraints',true)
	end % derive_surface_elevation

	end % methods
end % class Side_Weir

