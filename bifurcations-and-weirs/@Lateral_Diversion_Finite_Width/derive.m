% Mon 20 Jan 11:23:39 +08 2020
function [sol] = derive(oname)

		syms a b g n x y k real
		obj2  = Lateral_Diversion_Finite_Width('alpha',a,'beta',b,'gamma',g,'n',n);
	
		[u,v]         = obj2.velocity(x,y);
		sol.sym.u = u;
		sol.sym.v = v;
		%sol.sym.y
		sol.sym.du_dx = diff(u,x) 
		sol.sym.du_dy = diff(u,y)
		sol.sym.dv_dx = diff(v,x)
		sol.sym.dv_dy = diff(v,y)
	% this is invalid, frequencies cannot be added!
	%	sol.sym.R = -streamline_radius_of_curvature( ...
	%					    u,diff(u,x),diff(u,y), ...
	%					    v,diff(v,x),diff(v,y));

		%h = obj2.h

		% note, R, ub, vb and Jb are non-linear and cannot composed
		% of the series components of there parts
		% sum_k ub != ub(sum k)

	%	[ub,vb] = bend_velocity_near_bed(sol.sym.u,sol.sym.v,h,sol.sym.R);
	%	sol.sym.ub = ub;
	%	sol.sym.vb = vb;


	sol.sym.J  = [	sol.sym.du_dx, sol.sym.du_dy;
			sol.sym.dv_dx, sol.sym.dv_dy ];
%	sol.sym.Jb = [diff(ub,x), diff(ub,y);
%          		   diff(vb,x), diff(vb,y)];

		var = {x,y,a,b,g,k};

		d2u_dx2 = diff(diff(u,x),x);
		d2u_dxdy = diff(diff(u,x),y);
		d2u_dy2 = diff(diff(u,y),y);
		d2v_dx2 = diff(diff(v,x),x);
		d2v_dxdy = diff(diff(v,x),y);
		d2v_dy2 = diff(diff(v,y),y);

		sol.fun.k.u     = matlabFunction(sol.sym.u,'var',var);
		sol.fun.k.v     = matlabFunction(sol.sym.v,'var',var);
	%	sol.fun.k.ub    = matlabFunction(sol.sym.ub,'var',var);
	%	sol.fun.k.vb    = matlabFunction(sol.sym.vb,'var',var);
	%	sol.fun.k.R     = matlabFunction(sol.sym.R,'var',var);

		sol.fun.k.du_dx = matlabFunction( sol.sym.du_dx,'var',var);	
		sol.fun.k.du_dy = matlabFunction( sol.sym.du_dy,'var',var);	
		sol.fun.k.dv_dx = matlabFunction( sol.sym.dv_dx,'var',var);	
		sol.fun.k.dv_dy = matlabFunction( sol.sym.dv_dy,'var',var);	
		
		sol.fun.k.d2u_dx2 = matlabFunction( d2u_dx2,'var',var);	
		sol.fun.k.d2u_dxdy = matlabFunction( d2u_dxdy,'var',var);	
		sol.fun.k.d2u_dy2 = matlabFunction( d2u_dy2,'var',var);	
		sol.fun.k.d2v_dx2 = matlabFunction( d2v_dx2,'var',var);	
		sol.fun.k.d2v_dxdy = matlabFunction( d2v_dxdy,'var',var);	
		sol.fun.k.d2v_dy2 = matlabFunction( d2v_dy2,'var',var);	

	
		sol.fun.k.J  = matlabFunction(sol.sym.J,'var',var);
		%sol.fun.k.Jb = matlabFunction(sol.sym.Jb,'var',var);


%		sol.fun.k.du_dx = matlabFunction(sol.sym.du_dx,'var',var);
%		sol.fun.k.du_dy = matlabFunction(sol.sym.du_dy,'var',var);
%		sol.fun.k.dv_dx = matlabFunction(sol.sym.dv_dx,'var',var);
%		sol.fun.k.dv_dx = matlabFunction(sol.sym.dv_dy,'var',var);

%		fname_C = {'u','v','ub','vb','R'};
%		for idx=1:length(fname_C)
%			sol.fun.(fname_C{idx}) = @(x,y) obj.evalk(fname_C{'idx'},x,y);
%		end
		
		% recommended input distance on left end
		syms x a tol xin;
		syms y positive;
		if (1)
			% this is to conservative for w0->inf
			xin=solve(int(a/pi*y/(x^2 + y^2),x,-inf,xin)-tol*y,xin);
			xin = subs(xin,y,a);
			xin = series(xin,tol,0,'order',2);
			sol.sym.xin = xin;
			sol.fun.xin = matlabFunction(xin,'var',{'a','tol'});
			syms x0;
		else
			syms x a g y0;
			l.alpha = a;
			l.gamma = g;
			v = l.v_far(x,y0)
			v = series(v,y0,0,'order',2)
			% matlab fails to see the simplification here, manually applied
			v = (a*pi*g.^2*y0*(1 - coth((pi*g*x)/2)^2))/4;
			I=int(v,x);
			I=(I-subs(I,x,-inf));
			s=solve(I-y0*tol,x);
			If = matlabFunction(I)
			sf = matlabFunction(s)
		end

%		syms a g x y real; l.alpha = a; l.gamma = g;
		uf     = obj2.u_far(x,y);
		vf     = obj2.v_far(x,y);
		duf_dx = simplify(diff(uf,x));
		duf_dy = diff(uf,y);
		dvf_dx = diff(vf,x);
		dvf_dy = diff(vf,y);
		Rf     = -streamline_radius_of_curvature( ...
						    uf,duf_dx,duf_dy, ...
						    vf,dvf_dx,dvf_dy);
		[ubf,vbf] = bend_velocity_near_bed(uf,vf,beta/obj.fs,Rf);

		% near field
		if (0)
			syms a b g x y; %l.alpha = a; l.beta=b; l.gamma=g;
			[u,v] = l.velocity(x,y);
			du_dx0 = subs(simplify(subs(diff(u,x),x,0)),y,0);
			dy_dy0 = simplify(subs(subs(diff(v,y),x,0),y,0));
			% test continuity
			simplify(du_dx0 + dy_dy0)
		else
			% matlab fails to sum the series, the series solution is:
			du_dx0 = a*g*coth(pi*g/4);
			dv_dy0 = -du_dx0;
			un = 1 + du_dx0*x;
			vn =     dv_dy0*y;
			series(un,g,0,'order',3)
			%fdu=matlabFunction(dun)
			r  =  du_dx0/limit(du_dx0,g,0), double(subs(r,g,1))
			Rn  = streamline_radius_of_curvature(un,du_dx0,0,v,0,dv_dy0)
			% assume(y>0); R=simplify(subs(subs(R,x,0),y,a/2))
			sol.fun.unear = matlabFunction(un);
			sol.fun.vnear = matlabFunction(vn);
			sol.fun.Rnear = matlabFunction(Rn);
		end

		var = {x,y,a,g};
		sol.fun.duf_dx = matlabFunction(duf_dx,'var',var);
		sol.fun.duf_dy = matlabFunction(duf_dy,'var',var);
		sol.fun.dvf_dx = matlabFunction(dvf_dx,'var',var);
		sol.fun.dvf_dy = matlabFunction(dvf_dy,'var',var);
		sol.fun.Rf = matlabFunction(Rf,'var',var);

		var = {x,y,a,b,g};
		sol.fun.ubf= matlabFunction(ubf,'var',var);
		sol.fun.vbf= matlabFunction(vbf,'var',var);

		if (0)
		% his is more accurate, but there is hardly any difference
		% dy/dx = dy/dt*dt/dx = v/u, in far field:
		% dy/dx = a*y/(pi*(x^2+y^2) - a*x)
		% syms x y(x) y0 a tol; syms a positive; y=dsolve(diff(y,x) - a*y/(pi*x^2-a*x),y(-inf)==y0,x), x0=solve((y-y0) - y*tol, x), series(x0,tol,0,'order',2) 
		% double(subs(subs(subs(y,y0,1),a,1),x,-100))
		end

		if (nargin()>0 && ~isempty(oname))	
			save(oname,'sol');
		end
end % function load functions

