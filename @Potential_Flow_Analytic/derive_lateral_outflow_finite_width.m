% Tue 21 Aug 11:48:34 CEST 2018
%%
%% derive coefficients for lateral outflow in the case of potential flow
%%
function derive_lateral_outflow_finite_width(obj)

% TODO cases have to be differentiated for integration depending on the exponent
% k w0 < 1 and k w0 > 1
% as well as k y < 1 and k y > 1

% note that as for an infinite channel
% du/dy  = dv/dx
% u(0,y) = u0, du/dy|_0 = 0
% v(0,0) = 1, du/dx|_0 = = 

syms v0 x y k ws w0 a b L real
assume(L>0);
assume(ws>0);
assume(w0>0);
assume(v0>0);
v0 = v0;

	sol = struct();
	% constant profile
	% f = 1 in -w/2, w/2, 0 otherwise
	% int_-inf^inf = int_-1/2 w ^ 1/2 w
	% along the side of the outlet:
	% v = {0, x<-w/2; 1, -w/2<x<w/2; 0, x>w/2}
	sol.const.f = v0;
	% symmetric linear v-shaped profile
	sol.v.f    = 2*v0*((1-2*abs(x)/ws));
	sol.quad.f = 30*v0*(x/ws-1/2)^2*(x/ws+1/2)^2;
	
	field_C = fieldnames(sol);
	for fdx=1:length(field_C)
	
	field = field_C{fdx};
	f = sol.(field).f;
	
	cx  = v0*int(f*exp(-1i*k*x),-ws/2,ws/2)
	fx  = 1/sym('2*pi')*cx*exp(-1i*k*x)
	
	% across the channel
	% v = {0, y=-w0; 1, y = 0}
	w0_lim = w0;
	fy  = a*exp(k*y) + b*exp(-k*y);
	if (0)
		dfy = diff(fy,y);
		cy   = solve({subs(dfy,y,-w0)-0,subs(cy,y,0)-1},{a,b});
	else
	%	cy   = solve({limit(subs(fy,y,-w0)-0,w0,w0_lim), ...
		cy   = solve([subs(fy,y,-w0)-0, ...
			      subs(fy,y,0)-1],[a,b]);
	end
	
	fy  = simplify(subs(subs(fy,a,cy.a),b,cy.b))
	
	% for each frequency component of the continuous spectrum
	dv_dk = real(fx*fy)
	
	% integrated
	% this integral has no analytical solution
	v = int(dv_dk,k)
	
	% solution for phi
	% v = dphi/dx
	dphi_dk = simplify(int(dv_dk,y))
	
	% solution for u
	du_dk = simplify(diff(dphi_dk,x))
	du_dx_dk = diff(du_dk,x)
	dv_dy_dk = diff(dv_dk,y)
	
	disp('velocity at zero for infintely wide channel');
	du_dk_ = simplify(taylor(taylor(du_dk,y,0,'order',1),x,0,'order',3))
	p=1; 
	u_ = limit(int((1-0*abs(k)/L)*taylor(du_dk_*w0^p,w0,0,'order',1)/w0^p,k,-L,L),L,inf)
	
	disp('derivative of velocity along main channel for infintely wide channel');
	% logically, this is u0 + 
	du_dx_dk_ = simplify(taylor(taylor(diff(du_dk,x),x,0,'order',1),y,0,'order',1))
	p=1; 
	du_dx_ = limit(int((1-0*abs(k)/L)*taylor(du_dx_dk_*w0^p,w0,0,'order',1)/w0^p,k,-L,L),L,inf)
	
	dv_dx_ = 0; % constant vel profile along outflow
	dv_dy_ = -du_dx_; % by continuity
	du_dy_ = 0; % = dv/dx, by continutiy
	v_     = v0;
	
	disp('curvature at zero');
	R = streamline_radius_of_curvature(u_,du_dx_,du_dy_,v_,dv_dx_,dv_dy_)
	
	% cesaro sum, allows for integration without oscillation
	f=simplify(real((1-k/L)*subs(subs(du_dx_dk,x,0),y,0))) 
	
	% rewrite to avoid inf/inf
	%syms cosh(x); cosh(x) = coth(x)*sinh(x); s=subs(c(2));
	du_dk = simplify(subs(expand(du_dk),'cosh(k*w0)','coth(k*w0)*sinh(k*w0)'));
	[n,d] = numden(dv_dk);
	s = exp(2*k*w0)*exp(2*k*y);
	dv_dk = ((expand(n/s))/(expand(d/s)));
	
	sol.(field).du_dk = du_dk;
	sol.(field).dv_dk = dv_dk;
	sol.(field).du_dx_dk = simplify(du_dx_dk);
	sol.(field).dv_dy_dk = simplify(dv_dy_dk);
	
	du_dk    = matlabFunction(du_dk);
	dv_dk    = matlabFunction(dv_dk);
	du_dx_dk = matlabFunction(du_dx_dk);
	dv_dy_dk = matlabFunction(dv_dy_dk);
	sol.(field).fun.u = @(u0,v0,w0,ws,x,y,K,tol) quadgk(@(k) (1-abs(k)/K).*du_dk(k,v0,w0,ws,x,y),-K,-sqrt(eps),'AbsTol',tol) + quadgk(@(k) (1-abs(k)/K).*du_dk(k,v0,w0,ws,x,y),sqrt(eps),K,'AbsTol',tol);
	sol.(field).fun.v = @(u0,v0,w0,ws,x,y,K,tol) quadgk(@(k) (1-abs(k)/K).*dv_dk(k,v0,w0,ws,x,y),-K,-sqrt(eps),'AbsTol',tol) + quadgk(@(k) (1-abs(k)/K).*dv_dk(k,v0,w0,ws,x,y),sqrt(eps),K,'AbsTol',tol);
	sol.(field).fun.du_dx = @(u0,v0,w0,ws,x,y,K,tol) quadgk(@(k) (1-abs(k)/(K)).*du_dx_dk(k,v0,w0,ws,x,y), -K,-sqrt(eps),'AbsTol',tol) + quadgk(@(k) (1-abs(k)/(K)).*du_dx_dk(k,v0,w0,ws,x,y), sqrt(eps),K,'AbsTol',tol);
	sol.(field).fun.dv_dy = @(u0,v0,w0,ws,x,y,K,tol) quadgk(@(k) (1-abs(k)/(K)).*dv_dy_dk(k,v0,w0,ws,x,y), -K,-sqrt(eps),'AbsTol',tol) + quadgk(@(k) (1-abs(k)/(K)).*dv_dy_dk(k,v0,w0,ws,x,y), sqrt(eps),K,'AbsTol',tol);
	
	end % fdx

	obj.lateral_outflow_function('finite',sol);	

end % derive_lateral_outflow_finite_width

