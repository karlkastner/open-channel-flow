% Sun  7 May 12:09:35 CEST 2017
% Mon  4 Sep 15:25:48 CEST 2017
% c.f. chow 10.4
%
%% analytical solution to the gradually varied flow equation (bresse method)
%
function [x, h, zs, u] = solve_analytic(obj,Q,C,W,S0,h0,nn)
	n = 3;
	m = 3;

	g = 9.81;

	% normal flow depth
	% chow 10-14
	yn = normal_flow_depth(Q,W,C,S0);

	% critical flow depth
	% chow 10-13
	yc = critical_flow_depth(Q,W);

	% relative depth
	y = h0 + (yn-h0)*(0:nn-1)'/nn;
	u = y./yn;

	% note that (yc/yn)^2 = C^2/g*S0

%	u = (0:nn-1)'/nn;

	% TODO make a class member
	me = 2;
	
	switch (me)
	case {0}
	for idx=1:length(u)

	x(idx,1)  = yn/S0*(u(idx) ...
		  - quad(@(u_) 1./(1-u_.^n),0,u(idx)) ...
		   + (yc/yn)^m*quad(@(u_) u_.^(n-m)./(1-u_.^n),0,u(idx)) );
		%1./(1-linspace(0,u(idx),10).^n)
		u_ = linspace(0,u(idx),100);
		u_.^(n-m)./(1-u_.^n)
		quad(@(u_) 1./(1-u_.^n),0,u(idx))
		quad(@(u_) u_.^(n-m)./(1-u_.^n),0,u(idx))
		pause
	end
	case {1}
		% rectangular channel
		% bresse 1860
		% allen 1968, eq 3
		% note probable typo in lamb
		x = -yn./(6*S0).*(fun1(y)-fun1(y(1)));
	case {2}
		% eq 2.11 Conventional Integral Solutions of the GVF Equation
		% eq 10.12 in chow
		%x = -yn/S0*(u - (1-(yc/yn)^3)*(fun2(u)-fun2(u(1))));

		x = gvf_x_chow(u,yc,yn,C,S0) - gvf_x_chow(u(1),yc,yn,C,S0);
	end
%	x = fun2(u) - fun2(u(1));
%	x(:,2) =  fun2(u) - fun2(u(1));
%	[x y]
	% u = y/yn;
%	y =yn.*u;
	h = y;

        % bed level
        zb = S0*x - h0;

        % flow depth
        zs = bsxfun(@plus,h,zb);

	function f = fun1(y)
		f = (6*y./yn - (1-sqrt(yc./yn)^3).*(log((y.^2+1*y*yn+yn.^2)./(y-yn).^2) ...
				- 2/3*atan(sqrt(3)*yn./(2*y+yn))));
	end
end

