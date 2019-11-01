% 2016-04-08 11:15:55.270978358 +0200
% Karl Kastner, Berlin
%
%% change of depth along channel for the backwater equation
%
%  S0   : bed slope
%% beta : momentum coefficient
%% this is effectively an equation in h^3
%
% TODO improve for cases of high froude number and strong tidal flow
%      froude number correction only exact if tidal flow zero
%      and tidal friction only exact when froude number is low
% TODO interaction of tide-width is neglected
function dh_dx = dh_dx(obj,x,h)

	% acceleration by gravity
	g = obj.g;

	% momentum coefficient
	beta = obj.beta;

	Q0  = obj.Q0;
	Q1  = obj.Q1(x);
	Qhr = abs(Q1);
	p   = -obj.rt.friction_coefficient_dronkers(abs(Q0)./Qhr);

	S_b = obj.dzb_dx(x);
	
	C  = obj.chezy(x,h);

	W  = obj.width(x);
	dw_dx = obj.dw_dx(x);

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
		S_f = 1./C.^2/(pi*A0.^2.*R)*(p(1)*Qhr.^2 + p(2)*Q0.*Qhr + p(3)*(Q0.*abs(Q0) + 1/2*abs(Q1).^2));
	else
		S_f = U0.^2./(C^2.*R);
	end

	% due to width-variation
	% 2018 11 02 note S_w = 0, see paper 3 
	%S_w = h/W*dw_dx;
	S_w = 0;

	% change of flow depth along channel
	dh_dx = (S_f - S_w - S_b)./(1 - F2);
	%dh_dx = (S0 + S_f)./(1 - F2);
end

