% Tue 17 May 16:07:51 CEST 2016
% Karl Kastner, Berlin
%
%% determine bed level from surface elevation
%% (inverse backwater equation)
%% this is ill conditioned, as the surface is smooth for subcritical flow,
%% even if the bed is not smoth 
%%
%% C : chezy
%% W : width
%% Q : discharge
%% S : bed slope
%% y0 : surface elevation at outflow
function [x y_b] = inverse_backwater_curve(Q,C,W,ys,yb0,X)
	% dy/dx = 1./(1-F^2)*(S0 - S_f  - 2*beta*Q/(g*A^2)*I - Q^2/(g*A^2)*dbeta_dx)
	% beta = 1;
	% dbeta_dx = 0;
	%% lateral inflow
	% I = 0;
	
	% solve ODE
	opt.InitialStep = 0.1; % m
	opt.MaxStep = 100; % m
	% use 23 here, solution may have not have many continuous derivatives
	[x y_b] = ode23(@(x,yb) dyb_dx(x,yb,Q,C,ys,W), X, yb0, opt);
end % backwater

% S0   : bed slope
% beta : momentum coefficient
function dyb_dx = dyb_dx(x,yb,Q,C,ys,W)
	% acceleration by gravity
	g    = 9.81;
	% momentum coefficient
	beta = 1;
	if(isa(C,'function_handle'))
		C = feval(C,x);
	end
%	if(isa(ys,'function_handle'))
		[ys dys_dx] = feval(ys,x);
%	end
%	if(isa(dys_dx,'function_handle'))
%		dys_dx = feval(dys_dx,x);
%	end
	if(isa(W,'function_handle'))
		W = feval(W,x);
	end
	% flow depth
	h = ys - yb;
	h(h<0) = NaN;
	% area
	A   = h.*W;
	% hydraulic radius, wide channel approximation
	R   = h;
	% velocity
	U   = Q./A;
	% squared froude number
	F2  = beta * U.^2./(g*h);
	% friction slope
	S_f = sign(U).*U.^2./(C.^2.*R);
	% change of flow depth along channel
	%dh_dx = -(S0 - S_f)./(1 - F2);
	% TODO area change ???
	
	% derivative of bed slope along channel
	%dyb_dx = dh_dx.*(1-F2) + S_f;
	%dyb_dx = (S_f + dys_dx.*(1-F2))./(2 - F2);
	%dyb_dx = 1/(2 - F2)*((1-F2)*dys_dx - S_f);
	dyb_dx = -1./F2.*(S_f - dys_dx.*(1-F2));
end

