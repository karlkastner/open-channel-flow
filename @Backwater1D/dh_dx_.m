% 2016-04-08 11:15:55.270978358 +0200
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
function dh_dx = dh_dx_(obj,h,Q0,Qt,Qmid,Qhr,C,W,dzb_dx)

	% acceleration by gravity
	g = obj.g;

	% momentum coefficient
	beta = obj.beta;
	
	if (~isempty(Qt))
		% TODO attenuation of the tide
		p = -obj.rt.friction_coefficient_dronkers(abs(Q0)./Qhr);
	else
		p = -pi*[0,0,1];
		Qt = 0; 
		Qhr = 0; 
	end

	if (~issym(h))
		h(h<0) = sqrt(eps);
	end

	% area
	A0   = h.*W;

	% perimeter
	if (obj.widechannel)
		% wide channel
		P = W;
	else	
		% rectangular channel
		P   = 2*h+W;
	end
	% hydraulic radius
	R   = A0./P;

	% velocity
	U0   = Q0./A0;

	% squared froude number
	% momentum coefficient beta compensates for int_cs(u^2)/(int_cs u)^2
	F2  = beta * U0.^2./(g*h);

	% friction slope
	if (~issym(h))
		% S_f = U.*abs(U)./(C^2.*R);
		% note: cross term Q1*Q2, produces Q1 and Q3
		Qt2 = sum(abs(Qt).^2,2);
		S_f = 1./C.^2./(pi*A0.^2.*R).*(p(1)*Qhr.^2 + p(2)*Q0.*Qhr + p(3).*(Q0.*abs(Q0) + 1/2*Qt2));
	else
		S_f = U0.^2./(C^2.*R);
	end

	% due to width-variation
	% 2018 11 02 note S_w = 0, see paper 3 
	%S_w = h/W*dw_dx;
	%S_w = 0;

	% change of flow depth along channel
	%dh_dx = (S_f - S_w - S_b)./(1 - F2);
	dh_dx = (S_f - dzb_dx)./(1 - F2);
	%dh_dx = (-dzb_dx + S_f)./(1 - F2);
end

