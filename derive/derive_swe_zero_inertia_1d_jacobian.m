% Fri  4 Oct 14:57:13 CEST 2024
	syms zc zbl  bi zbr real
	syms C hc hl hr dx positive
	syms c1 c2 g positive

%	Sl = (h + z  - hl - zl)/dx;
%	Sr = (hr+ zr - h  - z )/dx;
	
	hi = (hl+hc)/2;

	Si = (hc + zc - hl - zbr)/dx;
	den = sqrt(c1^2 + 4*c2*g*abs(Si).*hi.^3);
	% velocity of interfaces
	ui  = sign(Si).*(c1 - den)./(2*c2*hi);
	qi    = 0.5*ui.*( (1-bi).*hl    + (1+bi).*hc);

	dui_dhl = diff(ui,hl)
	dui_dhc = diff(ui,hc)
	dqi_dhl = diff(qi,hl)
	dqi_dhc = diff(qi,hc)

dui_dhl_ =  - ui/(2*hi) - g*((1.5*Si*hi) - (sign(Si)^2*hi^2)/dx)/den
dui_dhc_ =  - ui/(2*hi) - g*((1.5*Si*hi) + (sign(Si)^2*hi^2)/dx)/den
simplify(dui_dhl - dui_dhl_)
simplify(dui_dhc - dui_dhc_)

% qi    = 0.5*ui.*( (1-bi).*hl    + (1+bi).*hc);
dqi_dhl_ = 0.5*dui_dhl_*(( (1-bi).*hl    + (1+bi).*hc)) + 0.5*ui.*(1-bi);
dqi_dhc_ = 0.5*dui_dhc_*(( (1-bi).*hl    + (1+bi).*hc)) + 0.5*ui.*(1+bi);
simplify(dqi_dhl_ - dqi_dhl)
simplify(dqi_dhc_ - dqi_dhc)

%dqi_dhl_ = - (dirac(Si*dx)*(qi/sign(Si)) - (sign(Si*dx)*(c1 - den)*(hc*(bi + 1) - hl*(bi - 1)))/(8*c2*hi^2) - (sign(Si*dx)*(c1 - den)*(bi - 1))/(4*c2*hi) - (sign(Si*dx)*(hc*(bi + 1) - hl*(bi - 1))*((6*c2*g*abs(Si*dx)*hi^2)/dx - (4*c2*g*sign(Si*dx)*hi^3)/dx))/(8*c2*den*hi)
%dqi_dhc_ = (dirac(Si*dx)*(qi/sign(Si)) - (sign(Si*dx)*(c1 - den)*(hc*(bi + 1) - hl*(bi - 1)))/(8*c2*hi^2) + (sign(Si*dx)*(c1 - den)*(bi + 1))/(4*c2*hi) - (sign(Si*dx)*(hc*(bi + 1) - hl*(bi - 1))*((6*c2*g*abs(Si*dx)*hi^2)/dx + (4*c2*g*sign(Si*dx)*hi^3)/dx))/(8*c2*den*hi)	


	% velocity, this is the coefficient
	%ul = -C*sign(Sl)*sqrt(0.5*(hl + h ).*abs(Sl)) 
	%ur = -C*sign(Sr)*sqrt(0.5*(h  + hr).*abs(Sr)) 

	

if (0)

Jl = diff(dh_dt,hl)
Jr = diff(dh_dt,hr)
Jc = diff(dh_dt,h )

	dul_dhl = diff(ul,hl)
	dur_dhr = diff(ur,hr)
	dul_dhc = diff(ul,h)
	dur_dhc = diff(ur,h)


	hli = 0.5*(hl + h);
	hri = 0.5*(hr + h);
	dul_dhl_ =  ul.*(0.25./hli - 0.5./(dx*Sl))
	dur_dhr_ =  ur.*(0.25./hri + 0.5./(dx*Sr))

	dul_dhc_ =  ul.*(0.25./hli + 0.5./(dx*Sl))
	dur_dhc_ =  ur.*(0.25./hri - 0.5./(dx*Sr))

	dh_dt = -(    0.5*ur.*( (1-br).*h   + (1+br).*hr) ...
	            - 0.5*ul.*( (1-bl).*hl  + (1+bl).*h ) ...
		 )/dx;


	%$dur_dhr = diff(ur,hr);
	Jl_  = (   0.5*ul.*(1-bl) ...
		      +0.5*dul_dhl_.*( (1-bl).*hl  + (1+bl).*h ) ...
		    )/dx;
	Jr_  = -(   0.5*ur.*(1+br) ...
		   +0.5*dur_dhr_*( (1-br).*h   + (1+br).*hr) ...
	        )/dx;

	Jc_ =  -( ...
		    0.5*ur.*(1-br) ...
		   +0.5*dur_dhc_.*((1-br).*h   + (1+br).*hr) ...
		   -0.5*ul.*(1+bl) ...
		   -0.5*dul_dhc_.*((1-bl).*hl  + (1+bl).*h ) ...
	       )/dx;
	Jc__ = -(Jl_ + Jr_);

du_dhl    = diff(ul,hl)
f.dul_dhl  = matlabFunction(dul_dhl);
f_.dul_dhl = matlabFunction(dul_dhl_);
f.dur_dhr  = matlabFunction(dur_dhr);
f_.dur_dhr = matlabFunction(dur_dhr_);
f.dur_dhc  = matlabFunction(dur_dhc);
f_.dur_dhc = matlabFunction(dur_dhc_);
f.Jl      = matlabFunction(Jl);
f_.Jl      = matlabFunction(Jl_);
f.Jr      = matlabFunction(Jr);
f_.Jr      = matlabFunction(Jr_);
f.Jc      = matlabFunction(Jc);
f_.Jc      = matlabFunction(Jc_);
f__.Jc      = matlabFunction(Jc__);
if (0)
f.dul_dhl(1,1,1,2,0,0)
f_.dul_dhl(1,1,1,2,0,0)  	
f.dul_dhl(1,1,2,1,0,0) 
f_.dul_dhl(1,1,2,1,0,0) 
'dur_dhr'
f.dur_dhr(1,1,1,2,0,0)
f_.dur_dhr(1,1,1,2,0,0)  	
f.dur_dhr(1,1,2,1,0,0) 
f_.dur_dhr(1,1,2,1,0,0) 
end
f.dur_dhc(1,1,1,2,0,0)
f_.dur_dhc(1,1,1,2,0,0)  	
f.dur_dhc(1,1,2,1,0,0) 
f_.dur_dhc(1,1,2,1,0,0) 
pause
%    (C,bl,br,dx,h,hl,hr,z,zl,zr)
f.Jc( 1, 0, 0, 1,1, 2, 3,0,0,0)
f_.Jc(1, 0, 0, 1,1, 2, 3,0,0,0)
f__.Jc(1, 0, 0, 1,1, 2, 3,0,0,0)
%    (C,bl,br,dx,h,hl,hr,z,zl,zr)
f.Jc( 1, 0, 0, 1,1, 2, 0.5,0,0,0)
f_.Jc(1, 0, 0, 1,1, 2, 0.5,0,0,0)
f__.Jc(1, 0, 0, 1,1, 2, 0.5,0,0,0)

f.Jc( 1, 0, 0, 1,2, 1, 0.5,0,0,0)
f_.Jc(1, 0, 0, 1,2, 1, 0.5,0,0,0)
f__.Jc(1, 0, 0, 1,2, 1, 0.5,0,0,0)
end
