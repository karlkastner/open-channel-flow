% Thu  9 Aug 14:24:52 CEST 2018
%
%% derive potential flow solution to lateral outlfow from an infinitely
%% wide main channel
%
function derive_lateral_outflow_analytic(obj)
	reset(symengine)

	% infinitessimal width, dirac delta function
	syms x y x_ y_ u0 v0 x0 Qs Q0 Qin ws r c s h real
	syms a b real;
	y0 = 0
	% real
	var = {x,y,u0,Qs,h,ws};
	v0  = Qs/(h*ws);
	fs  = 2/0.41^2;
	%assume(ws > 0)
	%assume(x ~= 0)
	assume(y<0)
	varg = {'ignoreanalyticconstraints',true};

	% potential of a point source
	%phi = 1/sym('pi')*log(sqrt(x^2 + y^2));

	%disp('test');
	%diff(u,x) + diff(v,y)

	sol           = struct();
%	u =  diff(phi,y);
%	v = -diff(phi,x);
%	sol.delta.phi = phi;
%	sol.delta.u   = u;
%	sol.delta.v   = v;
%	sol.delta.fun.phi = matlabFunction(phi,'var',var);
%	sol.delta.fun.u  = matlabFunction(u,'var',var);
%	sol.delta.fun.v  = matlabFunction(v,'var',var);

	% finite width
	% constant and v shaped velocity profile

	sol.delta.f = 'dirac';
	sol.const.f = 1;
	sol.v.f     = 2*(1-2*abs(x0)/ws);
	sol.quad.f  = 30*(x0/ws-1/2)^2*(x0/ws+1/2)^2;
	%ffun.cos   = (cos(2*pi*y0/ws)+1)
	% this is only for the side outflow, not including the main channel flow
	phi       = -1/sym('pi')*log(sqrt((x-x0)^2 + (y-y0)^2));

	field_C = fieldnames(sol)
	for idx=1:length(field_C)

		disp(['field: ',field_C{idx}]);
		f = sol.(field_C{idx}).f

		% integrating the velocity
if (0)
		u    = diff(phi,x);
		u    = simplify(u);
		v    = diff(phi,y);
		v    = simplify(v);
		u_   = u0 + simplify(int(f*u,x0,-ws/2,ws/2))
		v_   = simplify(int(f*v,x0,-ws/2,ws/2))
		%v__   = simplify(int(f*v,x0,-ws/2,0) ...
		%                 + int(f*v,x0,0,ws/2))
		phi_u = simplify(int(u,y))
end
		if (strcmp('delta',field_C{idx}))
			% dirac impulse, no integration required
			Phi = Qs*subs(phi,x0,0);
		else
		% direct integration of the velocity potential
		% int int f u dy0 dy = int f int u dy dy0 = int f phi dy0
		%Phi  = int(f*phi,x0)
		%Phi  = simplify(Phi);
		Phi  = v0*int_(f*phi,x0,-ws/2,ws/2)
		Phi  = simplify(Phi);
		end

		u    = u0 + diff(Phi,x);
		u    = simplify(u,varg{:})
		v    = diff(Phi,y);
		v    = simplify(v,varg{:})

		if (0)
		assume(x>0)
		v = simplify(v);
		syms('x','clear');
		assume(x,'real');
		end

if (0)
		disp('test');
		simplify(diff(u,x) + diff(v,y),varg{:})
		simplify(diff(u_,x) + diff(v_,y),varg{:})
		du = simplify(u-u_);
		simplify(combine(du,'log',varg{:}),varg{:})
		%simplify(u-u_,varg{:})
		%dv = v - v_;
		%simplify(combine(dv,'log',varg{:}))
		simplify(v-v_,varg{:})
end

		sol.(field_C{idx}).phi = Phi;
		sol.(field_C{idx}).u = u;
		sol.(field_C{idx}).v = v;

		sol.(field_C{idx}).fun.phi = matlabFunction(Phi,'var',var);
		sol.(field_C{idx}).fun.u   = matlabFunction(u,'var',var);
		sol.(field_C{idx}).fun.v   = matlabFunction(v,'var',var);

		% streamline curvature and near bed velocity

		R = -streamline_radius_of_curvature(u,diff(u,x),diff(u,y), ...
						    v,diff(v,x),diff(v,y));
		R = simplify(R);
		%R = -streamline_radius_of_curvature(v,diff(v,x),diff(v,y),u,diff(u,x),diff(u,y));
		sol.(field_C{idx}).R = R;
		sol.(field_C{idx}).fun.R = matlabFunction(R,'var',var);

		umag = sqrt(u.^2 + v.^2);
		Is = 11*umag./R;

		sol.(field_C{idx}).umag     = umag;
		sol.(field_C{idx}).fun.umag = matlabFunction(umag,'var',var);
		sol.(field_C{idx}).Is       = Is;
		sol.(field_C{idx}).fun.Is   = matlabFunction(Is,'var',var);

		% was plus
		[ub,vb] = bend_velocity_near_bed(u,v,h,R);

		sol.(field_C{idx}).ubed = ub;
		sol.(field_C{idx}).vbed = vb;

		varb            = var;
                %varb{end+1}     = h;
		sol.(field_C{idx}).fun.ubed = matlabFunction(ub,'var',varb);
		sol.(field_C{idx}).fun.vbed = matlabFunction(vb,'var',varb);

		% jacobian (for streamline computation)
		J = [diff(u,x), diff(u,y);
                     diff(v,x), diff(v,y)];
		sol.(field_C{idx}).J  = simplify(J);
		Jb = [diff(ub,x), diff(ub,y);
                      diff(vb,x), diff(vb,y)];
		sol.(field_C{idx}).Jb     = simplify(Jb);
		sol.(field_C{idx}).fun.J  = matlabFunction(J,'var',var);
		sol.(field_C{idx}).fun.Jb = matlabFunction(Jb,'var',varb);

		% stagnation point
		try
			y0 = solve(limit(u,y,0),x);
			y0 = y0(1);
	
			sol.(field_C{idx}).y0 = y0;
			sol.(field_C{idx}).fun.y0 = matlabFunction(y0,'var',var);
		catch e
			e
		end

		% near and far field linearisation
		try
		    ufar   = simplify(taylor(subs(subs(u,y,c*r),x,s*r),r,inf,'order',2));
		    ufar   = simplify(subs(subs(ufar,c,y/r),s,x/r));
		    vfar   = simplify(taylor(subs(subs(v,y,c*r),x,s*r),r,inf,'order',2));
		    vfar   = simplify(subs(subs(vfar,c,y/r),s,x/r));
		    %Rfar  = -streamline_radius_of_curvature(vfar,diff(vfar,x),diff(vfar,y),ufar,diff(ufar,x),diff(ufar,y))
		    Rfar  = -streamline_radius_of_curvature(ufar,diff(ufar,x),diff(ufar,y), ...
							    vfar,diff(vfar,x),diff(vfar,y))
		    sol.(field_C{idx}).ufar  = ufar;
		    sol.(field_C{idx}).vfar  = vfar;
		    sol.(field_C{idx}).fun.ufar  = matlabFunction(ufar,'var',var);
		    sol.(field_C{idx}).fun.vfar  = matlabFunction(vfar,'var',var);
		    sol.(field_C{idx}).Rfar  = Rfar;
		    sol.(field_C{idx}).fun.Rfar  = matlabFunction(Rfar,'var',var);
		catch
			e
		end

		% far field expansion
		try
		    xorder = 2;
		    yorder = 2;

		    unear  = taylor(taylor(u,y,0,'order',yorder),x,0,'order',xorder);
		    %vnear  = taylor((v0+simplify(limit(diff(v,x),x,0,'left'))*x),y,0,'order',xorder);
		    %vnear  = taylor((v0+simplify(limit(diff(v,),x,0,'left'))*x),y,0,'order',xorder);
		    %vnear  = taylor(taylor(v,y,0,'order',yorder),x,0,'order',xorder);
		    %v_ = v0;
		    assume(Qs>0);
		    assume(ws>0);
		    % order of differentiation does matter here
		    dv_dx = limit(limit(diff(v,x),x,0),y,0,'left');
		    dv_dy = limit(limit(diff(v,y),x,0),y,0,'left');
		    v00   = simplify(limit(limit(v,x,0),y,0,'left'));
		    vnear = v00 + dv_dx*x + dv_dy*y;
		    %syms('Qs','clear')
		    %syms('ws','clear')

		    %Rnear = -streamline_radius_of_curvature(vnear,diff(vnear,x),diff(vnear,y),unear,diff(unear,x),diff(unear,y))
		    Rnear = -streamline_radius_of_curvature(unear,diff(unear,x),diff(unear,y), ...
							    vnear,diff(vnear,x),diff(vnear,y))

		    sol.(field_C{idx}).unear = unear;
		    sol.(field_C{idx}).vnear = vnear;
		    sol.(field_C{idx}).Rnear = Rnear;

		    sol.(field_C{idx}).fun.unear = matlabFunction(unear,'var',var);
		    sol.(field_C{idx}).fun.vnear = matlabFunction(vnear,'var',var);
		    sol.(field_C{idx}).fun.Rnear = matlabFunction(Rnear,'var',var);
		catch e
			e
		end % try

	% normalization

	try
	Phi = simplify(subs(subs(subs(subs(sol.(field_C{idx}).phi,Qs,a*u0*h*ws),x,ws*x_),y,ws*y_),h,b*ws/fs)/u0,'ignoreanalyticconstraints',true);
	u  = simplify(subs(subs(subs(subs(sol.(field_C{idx}).u,Qs,a*u0*h*ws),x,ws*x_),y,ws*y_),h,b*ws/fs)/u0,'ignoreanalyticconstraints',true);
	v  = simplify(subs(subs(subs(subs(sol.(field_C{idx}).v,Qs,a*u0*h*ws),x,ws*x_),y,ws*y_),h,b*ws/fs)/u0,'ignoreanalyticconstraints',true);
	ub = simplify(subs(subs(subs(subs(sol.(field_C{idx}).ubed,Qs,a*u0*h*ws),x,ws*x_),y,ws*y_),h,b*ws/fs)/u0,'ignoreanalyticconstraints',true);
	vb = simplify(subs(subs(subs(subs(sol.(field_C{idx}).vbed,Qs,a*u0*h*ws),x,ws*x_),y,ws*y_),h,b*ws/fs)/u0,'ignoreanalyticconstraints',true);
	y0 = simplify(subs(sol.(field_C{idx}).y0,Qs,a*u0*ws*h)/ws);
	% test
%	n=1e4;
%	a_=randn(n,1);
%	b_=randn(n,1);
%	x__=randn(n,1);
%	y__=randn(n,1);
%	norm(u(a,b,rand(),rand(),x,y) - u(a,b,rand(),rand(),x,y))

	% matlab fails to simplify this, so manually cancel u0 and ws
	Phi= simplify( subs(subs(Phi,u0,1),ws,1),'ignoreanalyticconstraints',true);
	u  = simplify( subs(subs(u,u0,1),ws,1),'ignoreanalyticconstraints',true);
	v  = simplify( subs(subs(v,u0,1),ws,1),'ignoreanalyticconstraints',true);
	ub = simplify(subs(subs(ub,u0,1),ws,1),'ignoreanalyticconstraints',true);
	vb = simplify(subs(subs(vb,u0,1),ws,1),'ignoreanalyticconstraints',true);
	J  = [diff(u,x_), diff(u,y_);
              diff(v,x_), diff(v,y_)];
	Jb = [diff(ub,x_), diff(ub,y_);
              diff(vb,x_), diff(vb,y_)];
	var_ = {x_,y_,a,b};

	sol.(field_C{idx}).n.phi  = Phi;
	sol.(field_C{idx}).n.u    = u;
	sol.(field_C{idx}).n.v    = v;
	sol.(field_C{idx}).n.ubed = ub;
	sol.(field_C{idx}).n.vbed = vb;
	sol.(field_C{idx}).n.J    = J;
	sol.(field_C{idx}).n.Jb   = Jb;
	sol.(field_C{idx}).n.y0   = y0;

	sol.(field_C{idx}).nfun.phi  = matlabFunction(Phi,'var',var_);
	sol.(field_C{idx}).nfun.u    = matlabFunction(u,'var',var_);
	sol.(field_C{idx}).nfun.v    = matlabFunction(v,'var',var_);
	sol.(field_C{idx}).nfun.ubed = matlabFunction(ub,'var',var_);
	sol.(field_C{idx}).nfun.vbed = matlabFunction(vb,'var',var_);
	sol.(field_C{idx}).nfun.J    = matlabFunction(J,'var',var_);
	sol.(field_C{idx}).nfun.Jb   = matlabFunction(Jb,'var',var_);
	sol.(field_C{idx}).nfun.y0   = matlabFunction(y0,'var',var_);
	catch e
		e
	end % try

	end % for idx

	%obj.lateral_outflow_function('inf',sol);
	w     = what(class(obj));
	name  = 'lateral-outflow-functions.mat';
	path_ = [w.path,filesep(),name];
	save(path_,'sol');
end % function

% fails for abs
function Flr = int_(f,x,l,r)
	F = int(f,x);
	Flr =   (subs(F,x,r) - limit(F,x,(l+r)/2,'right')) ...
	      - (subs(F,x,l) - limit(F,x,(l+r)/2,'left'));
	%Flr = int(f,x,l,(r+l)/2) + int(f,x,(l+r)/2,r);
end


