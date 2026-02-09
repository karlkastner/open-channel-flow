% Tue  8 Oct 21:16:01 CEST 2024

syms C h hl hr hu hd dx1 dx2 positive
syms hul hur hdl hdr positive
syms bl br bu bd  z zl zr zu zd real
syms zul zur zdr zdl real

	% surface slope at left interfaces
	Slx   = (h + z  - hl - zl)/dx1;
	Srx   = (hr+ zr - h  - z)/dx1; % right(Slx)
	Suy   = (h + z  - hu - zu)/dx2;
	Sdy   = (hd+ zd - h  - z )/dx2; % down(Suy)
	
	% surface level above and below left interface midpoint
	zsl_u = 0.5*(hul + zul + hu + zu);
	zsl_d = 0.5*(hdl + zdl + hd + zd);

	% surface slope at left interface in y-direction
	Sly   = 0.5*(zsl_d - zsl_u)/dx2;

	zsr_u = 0.5*(hur + zur + hu + zu);
	zsr_d = 0.5*(hdr + zdr + hd + zd);

	% surface slope at left interface in y-direction
	Sry   = 0.5*(zsr_d - zsr_u)/dx2;

	% surface level left and right of above interface midpoint
	zsu_l = 0.5*(hul + zul + hl + zl);
	zsu_r = 0.5*(hur + zur + hr + zr);

	% slope in x-direction at up-interface
	Sux = 0.5*(zsu_r - zsu_l)/dx1;

	% surface level left and right of above interface midpoint
	zsd_l = 0.5*(hdl + zdl + hl + zl);
	zsd_r = 0.5*(hdr + zdr + hr + zr);

	% slope in x-direction at up-interface
	Sdx = 0.5*(zsd_r - zsd_l)/dx1;

	% slope magnitude at left interface
	Sl = hypot(Slx,Sly);
	Sr = hypot(Srx,Sry);
	Su = hypot(Sux,Suy);
	Sd = hypot(Sdx,Sdy);

	% velocity magnitude at left interface
	hli = 0.5*(hl+h);
	hri = 0.5*(hr+h);
	hui = 0.5*(hu+h);
	hdi = 0.5*(hd+h);
	Vl = C*sqrt(hli.*Sl);
	Vr = C*sqrt(hri.*Sr);
	Vu = C*sqrt(hui.*Su);
	Vd = C*sqrt(hdi.*Sd);

	% velocity in x-direction at left interface 
	vl = -Vl.*Slx./Sl;
	vr = -Vr.*Srx./Sr;
	vu = -Vu.*Suy./Su;
	vd = -Vd.*Sdy./Sd;

if (1)
	dvd_dhl = diff(vd,hl);
	dvd_dhl = vd*(diff(Vd,hl)/Vd + diff(Sdy,hl)/Sdy - diff(Sd,hl)/Sd);
	dvd_dhl = vd*(0.5*(diff(hdi,hl)/hdi + diff(Sd,hl)/hd)  + diff(Sdy,hl)/Sdy - diff(Sd,hl)/Sd);
	dvd_dhl = vd*(0.5* diff(hdi,hl)/hdi + diff(Sdy,hl)/Sdy - 0.5*diff(Sd,hl)/Sd);
	dvd_dhl = vd*( 0 + 0 - 0.5*diff(Sd,hl)/Sd);
	dvd_dhl = vd*( 0 + 0 - 0.5*(-Sdx/(4*Sd^2*dx1)));
	dvd_dhl = 0.125*vd*Sdx/(Sd^2*dx1);
	simplify(dvd_dhl-diff(vd,hl))
end
if (0)
	%dvl_dhul = diff(vl,hul);
	dvl_dhul = 0.25*vl*Sly/(Sl^2*dx2);
	simplify(dvl_dhul - diff(vl,hul))
	dvr_dhur = 0.25*vr*Sry/(Sr^2*dx2);
	simplify(dvr_dhur - diff(vr,hur))

	dvl_dhdl = -0.25*vl*Sly/(Sl^2*dx2);
	simplify(dvl_dhdl - diff(vl,hdl))
	dvr_dhdr = -0.25*vr*Sry/(Sr^2*dx2);
	simplify(dvr_dhdr - diff(vr,hdr))

	dvu_dhul = diff(vu,hul);
	dvu_dhul = vu*(0.25*Sux/(Su^2*dx1));
	simplify(dvu_dhul - diff(vu,hul))
	
	%dvu_dhur = -vu*(-0.5*diff(Sd,hur)/Sd);
if (0)
	dvu_dhur = vu*(diff(Vu,hur)/Vu + diff(Suy,hur)/Suy - diff(Su,hur)/Su)
	dvu_dhur = vu*(0.5*(diff(hui,hur)/hur + diff(Su,hur)/Su) + diff(Suy,hur)/Suy - diff(Su,hur)/Su)
	dvu_dhur = vu*(0.5* diff(hui,hur)/hur + diff(Suy,hur)/Suy - 0.5*diff(Su,hur)/Su)
	dvu_dhur = vu*( 0 + 0 - 0.5*diff(Su,hur)/Su)
end
	dvu_dhur = -0.25*vu*Sux/(Su^2*dx1);
	simplify(dvu_dhur - diff(vu,hur))

	dvd_dhdl = 0.25*vd*Sdx/(Sd^2*dx1);
	simplify(dvd_dhdl - diff(vd,hdl))
if (0)
	dvd_dhdr  = diff(vd,hdr);
	dvd_dhdr = vd*(diff(Vd,hdr)/Vd + diff(Sdy,hdr)/hdr - diff(Sd,hdr)/Sd);
	dvd_dhdr = vd*(0.5*(diff(hdi,hdr) + diff(Sd,hdr)/Sd) + diff(Sdy,hdr)/hdr - diff(Sd,hdr)/Sd);
	dvd_dhdr = vd*(0.5* diff(hdi,hdr) + diff(Sdy,hdr)/hdr - 0.5*diff(Sd,hdr)/Sd);
	dvd_dhdr = vd*(                 0 + 0 - 0.5*Sdx/(2*Sd^2*dx1));
end
	dvd_dhdr = -0.5*vd*Sdx/(2*Sd^2*dx1);
	simplify(dvd_dhdr - diff(vd,hdr))

if (0)
	dvd_dhd = diff(vd,hd);
	dvd_dhd = vd*(diff(Vd,hd)/Vd + diff(Sdy,hd)/Sdy - diff(Sd,hd)/Sd);
	dvd_dhd = vd*(0.5*(diff(hdi,hd)/hdi + diff(Sdy,hd)/Sdy) + diff(Sdy,hd)/Sdy - diff(Sd,hd)/Sd);
	dvd_dhd = vd*(0.5*diff(hdi,hd)/hdi + diff(Sdy,hd)/Sdy - 0.5*diff(Sd,hd)/Sd);
end
	dvd_dhd = vd*(0.25/hdi + 1/(dx2*Sdy) - 0.5*(Sdy/(dx2*Sd^2)));
	simplify(dvd_dhd - diff(vd,hd))

if (0)
	dvr_dhd = diff(vr,hd);
	dvr_dhd = vr*(diff(Vr,hd)/Vr + diff(Srx,hd)/Srx - diff(Sr,hd)/Sr)
	dvr_dhd = vr*(0.5*(diff(hri,hd)/hri + diff(Sr,hd)/Sr) + diff(Srx,hd)/Srx - diff(Sr,hd)/Sr)
	dvr_dhd = vr*(0.5*diff(hri,hd)/hri  + diff(Srx,hd)/Srx - 0.5*diff(Sr,hd)/Sr)
	dvr_dhd = vr*(0  + 0 - 0.5*(Sry/(2*Sr^2*dx2)))
end
	dvr_dhd = -vr*0.5*(Sry/(2*Sr^2*dx2))

	simplify(dvr_dhd - diff(vr,hd))
if (0)
	dvl_dhd = vl*(diff(Vl,hd)/Vl + diff(Slx,hd)/Slx - diff(Sl,hd)/Sl)
	dvl_dhd = vl*(0.5*(diff(hli,hd)/hli + diff(Sl,hd)/Sl) + diff(Slx,hd)/Slx - diff(Sl,hd)/Sl)
	dvl_dhd = vl*(0.5* diff(hli,hd)/hli + diff(Slx,hd)/Slx - 0.5*diff(Sl,hd)/Sl)
	dvl_dhd = vl*(0 + 0 - 0.5*Sly/(2*Sl^2*dx2))
end
	dvl_dhd = vl*(-0.5*Sly/(2*Sl^2*dx2));
	simplify(dvl_dhd - diff(vl,hd))

if (0)
	dvu_dhu = diff(vu,hu);
	dvu_dhu = vu*(diff(Vu,hu)/Vu + diff(Suy,hu)/Suy - diff(Su,hu)/Su)
	dvu_dhu = vu*(0.5*(diff(hui,hu)/hui + diff(Su,hu)/Su) + diff(Suy,hu)/Suy - diff(Su,hu)/Su)
	dvu_dhu = vu*(0.5*diff(hui,hu)/hui + diff(Suy,hu)/Suy - 0.5*diff(Su,hu)/Su)
	dvu_dhu = vu*(0.25/hui + -1/(dx2*Suy) - 0.5*(-Suy/(dx2*Su^2)))
end
	dvu_dhu = vu*(0.25/hui + -1/(dx2*Suy) + 0.5*(Suy/(dx2*Su^2)));
	simplify(dvu_dhu - diff(vu,hu))	
if (0)
	dvr_dhu = diff(vr,hu);
	dvr_dhu = vr*(diff(Vr,hu)/Vr + diff(Srx,hu)/Srx - diff(Sr,hu)/Sr)
	dvr_dhu = vr*(0.5*(diff(hri,hu)/hri + diff(Sr,hu)/Sr) + diff(Srx,hu)/Srx - diff(Sr,hu)/Sr);
	dvr_dhu = vr*(0.5* diff(hri,hu)/hri                   + diff(Srx,hu)/Srx - 0.5*diff(Sr,hu)/Sr);
	dvr_dhu = vr*(                    0                   +                0 - 0.5*((-Sry/(2*Sr^2*dx2))));
end
	dvr_dhu = vr*0.25*(Sry/(Sr^2*dx2));
	simplify(dvr_dhu - diff(vr,hu))

if (0)	
	dvl_dhu = diff(vl,hu);
	dvl_dhu = vl*(diff(Vl,hu)/Vl + diff(Slx,hu)/Slx - diff(Sl,hu)/Sl);
end
	dvl_dhu = vl*(0.5*(diff(hli,hu)/hli + diff(Sl,hu)/Sl) + diff(Slx,hu)/Slx - diff(Sl,hu)/Sl);
	dvl_dhu = vl*(0.5* diff(hli,hu)/hli + diff(Slx,hu)/Slx - 0.5*diff(Sl,hu)/Sl);
	dvl_dhu = vl*( 0 + 0 - 0.5*(Sly/(-2*Sl^2*dx2)));
	dvl_dhu = 0.25*vl*Sly/(Sl^2*dx2);
%	dvl_dhu = vl*(0.5*( diff(hli,hu)/hli + diff(Sl,hu)/Sl ) + diff(Slx,hu)/Slx - diff(Sl,hu)/Sl)
%	dvl_dhu = vl*(0.5*diff(hli,hu)/hli + diff(Slx,hu)/Slx - 0.5*diff(Sl,hu)/Sl)
	simplify(dvl_dhu - diff(vl,hu))
if (0)
	dvu_dhr = diff(vu,hr);
	dvu_dhr = vu*(diff(Vu,hr)/Vu + diff(Suy,hr)/Suy - diff(Su,hr)/Su);
	dvu_dhr = vu*( 0.5*(diff(hui,hr)/hui + diff(Su,hr)/Su) + diff(Suy,hr)/Suy - diff(Su,hr)/Su);
	dvu_dhr = vu*( 0.5* diff(hui,hr)/hui + diff(Suy,hr)/Suy - 0.5*diff(Su,hr)/Su);
end
	dvu_dhr = vu*(                     0 +                0 - 0.5*(Sux/(2*dx1*Su^2)));
	dvu_dhr = -0.25*vu*Sux/(Su^2*dx1);
	simplify(dvu_dhr - diff(vu,hr))

if (0)
	dvr_dhr = diff(vr,hr);
	dvr_dhr = vr*(diff(Vr,hr)/Vr + diff(Srx,hr)/Srx - diff(Sr,hr)/Sr)
	dvr_dhr = vr*(0.5*(diff(hri,hr)/hri + diff(Sr,hr)/Sr) + diff(Srx,hr)/Srx - diff(Sr,hr)/Sr)
	dvr_dhr = vr*(0.5*(diff(hri,hr)/hri) + diff(Srx,hr)/Srx - 0.5*diff(Sr,hr)/Sr)
	dvr_dhr = vr*(0.5*(diff(hri,hr)/hri) + 1/dx1/Srx - 0.5*(Srx/Sr^2/dx1))
end	
	dvr_dhr = vr*(0.25/hri + 1/dx1/Srx - 0.5*(Srx/Sr^2/dx1))
%	dvr_dhr = vr*(0.5*(0.5/hri) + 1/dx1/Srx - 0.5*(Srx/Sr/dx1))
	simplify(dvr_dhr - diff(vr,hr))

	dvd_dhr = diff(vd,hr);
	dvd_dhr = vd*(diff(Vd,hr)/Vd + diff(Sdy,hr)/Sdy - diff(Sdy,hr)/Sd);
	dvd_dhr = vd*(0.5*(diff(hdi,hr)/hdi + diff(Sd,hr)/Sd) + diff(Sdy,hr)/Sdy - diff(Sd,hr)/Sd);
	dvd_dhr = vd*(0.5*diff(hdi,hr)/hdi + diff(Sdy,hr)/Sdy - 0.5*diff(Sd,hr)/Sd);
	dvd_dhr = vd*(                   0 +                0 - 0.5*(Sdx/(2*dx1*Sd^2)))
	dvd_dhr = -0.25*vd*Sdx/(dx1*Sd^2)
	simplify(dvd_dhr - diff(vd,hr))



if (0)
	dVl_dhl_ = 1/2*Vl*(1/hli*diff(hli,hl) + diff(Sl,hl)/Sl)
	dVl_dhl = diff(Vl,hl)
	dvl_dhl = vl*(diff(Vl,hl)./Vl + diff(Slx,hl)./Slx - diff(Sl,hl)/Sl)
	dvl_dhl = vl*(0.5*(diff(hli,hl)/hli + diff(Sl,hl)/Sl) + diff(Slx,hl)./Slx - diff(Sl,hl)/Sl)
	dvl_dhl = vl*(0.5/hli*diff(hli,hl) + diff(Slx,hl)./Slx - 1/2*diff(Sl,hl)/Sl)
end
	dvl_dhl = vl*(0.25/hli + (-1/dx1)./Slx + Slx/(2*Sl^2*dx1));
	simplify(dvl_dhl-diff(vl,hl))
if (0)
%	J_.dvl_dhl = vl*(dVl_dhl_./Vl + (-1/dx1)./Slx + Slx/(Sl^2*dx1))
	dvu_dhl = diff(-Vu.*Suy./Su,hl)
	dvu_dhl = vu*(diff(Vu,hl)/Vu + diff(Suy,hl)/Suy - diff(Su,hl)/Su)
	dvu_dhl = vu*(0.5*(diff(hui,hl)/hui + diff(Su,hl)/Sl ) + diff(Suy,hl)/Suy - diff(Su,hl)/Su)
	dvu_dhl = vu*(0.5*diff(hui,hl)/hui + diff(Suy,hl)/Suy - 0.5*diff(Su,hl)/Su)
	dvu_dhl = vu*(0 + 0 - 0.5*-Sux/(2*dx1*Su^2) );
end
	dvu_dhl = 0.25*vu*Sux/(4*dx1*Su^2);
	simplify(dvu_dhl-diff(vu,hl))
end
	dvd_dhd  =  vd*(0.25/hdi + (1/dx2)./Sdy - 0.5*Sdy/(Sd^2*dx2));
	dvl_dhl  =  vl*(0.25/hli - (1/dx1)./Slx + 0.5*Slx/(Sl^2*dx1));
	dvr_dhr  =  vr*(0.25/hri + (1/dx1)./Srx - 0.5*Srx/(Sr^2*dx1));
	dvu_dhu  =  vu*(0.25/hui - (1/dx2)./Suy + 0.5*Suy/(Su^2*dx2));
%	simplify(dvd_dhd - diff(vd,hd))
%	simplify(dvl_dhl - diff(vl,hl))
	dvd_dhdl =  0.125*vd*Sdx/(Sd^2*dx1);
	dvd_dhdr = -0.125*vd*Sdx/(Sd^2*dx1);
	dvd_dhl  =  0.125*vd*Sdx/(Sd^2*dx1);
	dvd_dhr  = -0.125*vd*Sdx/(dx1*Sd^2)
	dvl_dhd  = -0.125*vl*Sly/(Sl^2*dx2);
	dvl_dhdl = -0.125*vl*Sly/(Sl^2*dx2);
	dvl_dhu  =  0.125*vl*Sly/(Sl^2*dx2);
	dvl_dhul =  0.125*vl*Sly/(Sl^2*dx2);
	dvr_dhd  = -0.125*vr*Sry/(Sr^2*dx2);
	dvr_dhdr = -0.125*vr*Sry/(Sr^2*dx2);
	dvr_dhu  =  0.125*vr*Sry/(Sr^2*dx2);
	dvr_dhur =  0.125*vr*Sry/(Sr^2*dx2);
	dvu_dhl  =  0.125*vu*Sux/(dx1*Su^2);
	dvu_dhr  = -0.125*vu*Sux/(Su^2*dx1);
	dvu_dhul =  0.125*vu*Sux/(Su^2*dx1);
	dvu_dhur = -0.125*vu*Sux/(Su^2*dx1);

%	simplify(dvl_dhl-diff(vl,hl))
%	simplify(dvu_dhl-diff(vu,hl))

	dh_dt =	-(...
		   ( ...
		       0.5*      vr.*( (1-br).*h  + (1+br).*hr) ...
	             - 0.5*      vl.*( (1-bl).*hl + (1+bl).*h) ...
		   )/dx1 ...
		 + ( ...
		       0.5*     vd.*( (1-bd).*h          + (1+bd).*hd) ...
	             - 0.5*     vu   .*( (1-bu).*hu    + (1+bu).*h) ...
		   )/dx2  ...
		)


	% clockwise starting with top
	clear J_ J
	% up
	J_.hu =  0.5*dvu_dhu.*( (1-bu).*hu  + (1+bu).*h)/dx2 ...
		 +0.5*vu*(1-bu)/dx2 ...
		-0.5*dvr_dhu.*( (1-br).*h  + (1+br).*hr)/dx1 ... 
		+0.5*dvl_dhu.*( (1-bl).*hl + (1+bl).*h)/dx1;            

	% up right
	J_.hur = -0.5*dvr_dhur* ((1-br).*h   + (1+br).*hr)/dx1 ...
		 +0.5*dvu_dhur*( (1-bu).*hu + (1+bu).*h)/dx2;

	% right
	J_.hr   = -0.5*dvr_dhr.*( (1-br).*h  + (1+br).*hr)/dx1 ...
		  - 0.5*vr*(1+br)/dx1 ...
		  -0.5*dvd_dhr.*( (1-bd).*h  + (1+bd).*hd)/dx2 ...
		  +0.5*dvu_dhr.*( (1-bu).*hu + (1+bu).*h )/dx2;
	% d-r	
	J_.hdr = -0.5*dvr_dhdr.*( (1-br).*h + (1+br).*hr)/dx1 ...
		 -0.5*dvd_dhdr.*( (1-bd).*h + (1+bd).*hd)/dx2;
	% d
	J_.hd = -0.5*dvd_dhd.*( (1-bd).*h  + (1+bd).*hd)/dx2 ...
		- 0.5*vd*(1+bd)/dx2 ...
		-0.5*dvr_dhd.*( (1-br).*h  + (1+br).*hr)/dx1 ...
		+0.5*dvl_dhd.*( (1-bl).*hl + (1+bl).*h)/dx1;

	% down left
	J_.hdl   = +0.5*dvl_dhdl.*((1-bl).*hl + (1+bl).*h)/dx1 ...
		   -0.5*dvd_dhdl.*((1-bd).*h  + (1+bd).*hd)/dx2;

	% left
	J_.hl    = +0.5*dvl_dhl.*( (1-bl).*hl + (1+bl).*h)/dx1 ...
		   +0.5*vl*(1-bl)/dx1 ...
		   -0.5*dvd_dhl.*( (1-bd).*h  + (1+bd).*hd)/dx2 ...
	           +0.5*dvu_dhl.*( (1-bu).*hu + (1+bu).*h)/dx2;

	% top left
	J_.hul   = +0.5*dvl_dhul*((1-bl).*hl + (1+bl).*h)/dx1 ...
		   +0.5*dvu_dhul*((1-bu).*hu + (1+bu).*h)/dx2;

	% ul u vr
	% l  h  r
	% dl d dr

vars = {C,dx1,dx2,h,hd,hdl,hdr,hl,hr,hu,hul,hur,z,zd,zdl,zdr,zl,zr,zu,zul,zur,bd,bl,br,bu};

field_C = {'hl','hr','hu','hd','hul','hur','hdl','hdr'}
for idx=1:length(field_C)
	f = field_C{idx};
	J.(f) = diff(dh_dt,f);
end


field_C2 = {'Slx','Sly','Sl','Sux','Suy','Su','vl','vu', ...
		'dvd_dhd', 'dvl_dhl', 'dvr_dhr', 'dvu_dhu', 'dvd_dhdl', 'dvd_dhdr', 'dvd_dhl', 'dvd_dhr', 'dvl_dhd', 'dvl_dhdl', 'dvl_dhu', 'dvl_dhul', 'dvr_dhd', 'dvr_dhdr', 'dvr_dhu', 'dvr_dhur', 'dvu_dhl', 'dvu_dhr', 'dvu_dhul', 'dvu_dhur'}
%field_C = {'vl','vu'}
for idx=1:length(field_C2)
	f = field_C2{idx};
	J.(f) = eval(f);
	J_.(f) = eval(f);
end

field_C = {'hl','hr','hu','hd','hul','hur','hdl','hdr',field_C2{:}}
for idx=1:length(field_C)
	f = field_C{idx};
	fun.(f) = matlabFunction(J.(f),'Vars',vars);
	fun_.(f) = matlabFunction(J_.(f),'Vars',vars);
end
%J.dul_dhl  = diff(ul,hl)
%field_C = {'dul_dhl'}
%J.hul = diff(dh_dt,hul)
%J.hdl = diff(dh_dt,hdl)
%J.hl  = diff(dh_dt,hl)
%J.hr  = diff(dh_dt,hl)
%J.h   = diff(dh_dt,h)

%simplify(J.hul - J.hul_) 
%simplify(J.hdl - J.hdl_) 

rng=1;
x={};
for idx=1:30;
	x{idx} = rand();
end
if (1)
C   = 1; % rand()
dx1 = 1; %rand();
dx2 = 1;% rand();

h   = rand();
hd  = rand();
hl  = rand();
hdl = rand();
hdr = rand();
hr  = rand();
hu  = rand();
hul = rand();
hur = rand();

e = sqrt(eps);
if (0)
hul = e;
hu  = e;
hur = e;
hl  = 1+e;
h   = 1;
hr  = 1+e;
hdl = 2;
hd  = 2;
hdr = 2;
end
if (0)
hul = 2*e;
hu  = 1+e;
hur = 2+e;
hl  = e;
h   = 1;
hr  = 2;
hdl  = 2*e;
hd   = 1+e;
hdr  = 2+e;
end

z   = 0*rand();
zd  = 0*rand();
zdl = 0*rand();
zdr = 0*rand();
zl  = 0*rand();
zr  = 0*rand();
zu  = 0*rand();
zul = 0*rand();
zur = 0*rand();
bd = 0;
bl = 0;
br = 0;
bu = 0;
end
hh = [hul,hu,hur;
     hl,h,hr;
     hdl,hd,hdr];
zz = [zul, zu, zur;
      zl, z, zr;
      zdl,zd,zdr];  



x= {C,dx1,dx2,h,hd,hdl,hdr,hl,hr,hu,hul,hur,z,zd,zdl,zdr,zl,zr,zu,zul,zur,bd,bl,br,bu};

for idx=1:length(field_C)
	f = field_C{idx};
	n = nargin(fun.(f));
	m = nargin(fun_.(f));
	val_.(f) = fun_.(f)(x{1:n});
	val.(f) = fun.(f)(x{1:m});
	err = val_.(f) - val.(f); 
	printf('%s %g\n',f,err);
end

% Test

n = [3,3];
L = n.*[dx1,dx2];
[Jn,out] = swe_zero_inertia_2d_jacobian((hh),zz,C,L,n,'central',false)
Jn(5,:)
[val.hu val.hur val.hr val.hdr val.hd val.hdl val.hl val.hul]
%     f                             f 
out
val

% left to right
%    2.8284    0.0000    2.7747   -0.0000    2.8284    0.0000    0.2094   -0.0000   -8.6409
%    2.8284         0    3.1623         0    2.8284         0    0.0000         0
% up-down
%   -1.0607    0.0000    1.0000   -0.0000    3.0619   -0.0000    1.0000    0.0000   -4.0012
%   -0.3536    0.0000    1.0000   -0.0000    1.8371   -0.0000    1.0000    0.0000

