% Mon  2 Feb 10:48:01 CET 2026
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
%
% function f = drag_coefficient_colebrooke_white(ks,h,u)
%
% this is for turbulent flow close to the laminar transition
% to distinguish rough and smooth pipes
% the values do not approach the laminar limit f ~ 92/Re
%
function [f,f0,finf] = drag_coefficient_colebrooke_white(ks,h,u,isopen)
	nu = 1e-6;
	% the reference length for colebrook white is depth h
	% and not the hydraulic diamerter 4 h
	Re = u.*h./nu;
	if (isopen)
		Re = 4*Re;
		% TODO also scale ks/h?
	end
	es = ks./h;
	c2 = es/3.7;
	c3 = 2.51./Re;
	if (0)
	c0 = 0;
	c1 = 2/log(10);
	one_over_sqrt_f = c1.*lambertw(1./(c1.*c2).*exp((1./c1).*(c0 + (c2./c3)))) - c2./c3;
	else
		one_over_sqrt_f = zeros(size(c3));
		one_over_sqrt_f_inf = zeros(size(c3));
		one_over_sqrt_f_0 = zeros(size(c3));
		for idx=1:numel(c3)
			% laminar resistance as initial condition
			y0 = -2*log10(c2(idx));
			one_over_sqrt_f_0(idx) = y0;
			% at c3 = inf
			yinf = (2*wrightOmega(-log((2*c3(idx))/log(10))))/log(10);
			y0 = yinf;
			one_over_sqrt_f_inf(idx) = yinf;
			one_over_sqrt_f(idx) = lsqnonlin(@(y) y - -2*log10(c2(idx) + c3(idx).*y),y0);
		end
	end
	f = 1./one_over_sqrt_f.^2;
	finf = 1./one_over_sqrt_f_inf.^2;
	f0 = 1./one_over_sqrt_f_0.^2;
	%y0_ = 1./c3;
	f0(:,2) = c3.^2;
	%1./y0_.^2;
end

