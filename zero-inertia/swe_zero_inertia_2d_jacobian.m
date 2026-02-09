% Tue  8 Oct 21:18:13 CEST 2024
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

function [J,out] = swe_zero_inertia_2d_jacobian(h,z,C,L,n,mode,returnmat)
	e = sqrt(eps);
	dx  = L./n;
	
	h = reshape(h,n);
	z = reshape(z,n);

	hu  = up(h,1);
	hd  = down(h,1);
	hl  = left(h,1);
	hr  = right(h,1);
	zu  = up(z,1);
	zd  = down(z,1);
	zl  = left(z,1);
	zr  = right(z,1);
	hdl = left(down(h,1),1);
	hdr = right(down(h,1),1);
	hul = left(up(h,1),1);
	hur = right(up(h,1),1);
	zdl = left(down(z,1),1);
	zdr = right(down(z,1),1);
	zul = left(up(z,1),1);
	zur = right(up(z,1),1);

	% surface slope in x-direction at left interface
	Slx   = (h + z  - hl - zl)/dx(1);
	Srx = right(Slx,1);
	% surface slope in y-direction at up-interface
	Suy   = (h + z  - hu - zu)/dx(2);
	Sdy = down(Suy,1);

	% surface level above and below left interface midpoint
	zsl_u = 0.5*(hul + zul + hu + zu);
	zsl_d = 0.5*(hdl + zdl + hd + zd);

	% surface slope in y-direction at left interface
	Sly   = 0.5*(zsl_d - zsl_u)/dx(2);
	Sry   = right(Sly);


	% surface level left and right of up interface midpoint
	zsu_l = 0.5*(hul + zul + hl + zl);
	zsu_r = 0.5*(hur + zur + hr + zr);

	% slope in x-direction at up-interface
	Sux   = 0.5*(zsu_r - zsu_l)/dx(1);
	Sdx = down(Sux,1);

	% slope magnitude at left interface
	Sl  = sqrt(e + Slx.^2 + Sly.^2);
	Sr = right(Sl,1);
	% slope magnitude at up interface
	Su  = sqrt(e + Sux.^2 + Suy.^2);
	Sd = down(Su,1);

	% depth at interfaces
	hli = 0.5*(hl+h);
	hri = 0.5*(hr+h);
	hui = 0.5*(hu+h);
	hdi = 0.5*(hd+h);

	% velocity magnitude at left interface
	Vl  = C*sqrt(hli.*Sl);
	% velocity-magnitude at up-interface
	Vu  = C*sqrt(hui.*Su);

	% velocity in x-direction at left interface 
	vl  = -Vl.*Slx./Sl;
	vr  = right(vl,1);
	% velocity in y-direction at up-interface 
	vu  = -Vu.*Suy./Su;
	vd  = down(vu,1);

if (0)
close all
subplot(2,2,1)
imagesc(reshape(vl,n))
subplot(2,2,2)
imagesc(reshape(vu,n))

subplot(2,2,3)
imagesc(reshape(Slx,n))
subplot(2,2,4)
imagesc(reshape(Suy,n))
pause
end


	% celerity, advection coefficient at interface
	al = 1.5*vl;
	au = 1.5*vu;

	% diffusion coefficient at left interface
	dl = -0.25*(h+left(h,1)).*Vl./Sl; % was v instead of V
	du = -0.25*(h+  up(h,1)).*Vu./Su;

	% mesh peclet number
	s = -1;
	Pel = -s*(0.5*dx(1))*al./dl;
	Peu = -s*(0.5*dx(2))*au./du;

	% weights of difference terms at left interface
	switch (mode)
	case {'central'}
		bl = zeros(n);
		bu = zeros(n);
	case{'upwind'}
		bl = 2*(al<0)-1;
		bu = 2*(au<0)-1;
	case{'optimal'}
		% magical coefficient
		bl = coth(Pel) - 1./Pel;
		bu = coth(Peu) - 1./Peu;
	case {'optimal-linearized'}
		% linearized magical coefficient
		bl = max(-1,min(1,Pel/3));
		bu = max(-1,min(1,Peu/3));
	otherwise
		error('here')
	end
	fdx = isnan(bl);
	bl(fdx) = 2*(al(fdx)<0)-1;
	br = right(bl,1);
	fdx = isnan(bu);
	bu(fdx) = 2*(au(fdx)<0)-1;
	bd = down(bu,1);


	% derivatives
	%dvd_dhd  =  vd.*(0.25./hdi + 1./(dx(2).*Sdy) - 0.5.*(Sdy./(dx(2).*Sd.^2)));
	dvd_dhd  =  vd.*(0.25./hdi - 0.5.*(Sdy./(dx(2).*Sd.^2))) - down(Vu./(dx(2)*Su),1);
	dvl_dhl  =  vl.*(0.25./hli + 0.5*Slx./(Sl.^2.*dx(1))) + Vl./(dx(1)*Sl);
	%dvr_dhr  =  vr.*(0.25./hri + 1./dx(1)./Srx - 0.5.*(Srx./Sr.^2./dx(1)));
	dvr_dhr  =  vr.*(0.25./hri - 0.5.*(Srx./Sr.^2./dx(1))) - right(Vl./(dx(1)*Sl),1);
	dvu_dhu  =  vu.*(0.25./hui + 0.5.*Suy./(dx(2).*Su.^2)) + Vu./(dx(2)*Su);
	dvd_dhdl =  0.125.*vd.*Sdx./(Sd.^2.*dx(1));
	dvd_dhdr = -0.125.*vd.*Sdx./(Sd.^2.*dx(1));
	dvd_dhl  =  0.125.*vd.*Sdx./(Sd.^2.*dx(1));
	dvd_dhr  = -0.125.*vd.*Sdx./(Sd.^2.*dx(1));
	dvl_dhd  = -0.125.*vl.*Sly./(Sl.^2.*dx(2));
	dvl_dhdl = -0.125.*vl.*Sly./(Sl.^2.*dx(2));
	dvl_dhu  =  0.125.*vl.*Sly./(Sl.^2.*dx(2));
	dvl_dhul =  0.125.*vl.*Sly./(Sl.^2.*dx(2));
	dvr_dhd  = -0.125.*vr.*Sry./(Sr.^2.*dx(2));
	dvr_dhdr = -0.125.*vr.*Sry./(Sr.^2.*dx(2));
	dvr_dhu  =  0.125.*vr.*Sry./(Sr.^2.*dx(2));
	dvr_dhur =  0.125.*vr.*Sry./(Sr.^2.*dx(2));
	dvu_dhl  =  0.125.*vu.*Sux./(Su.^2.*dx(1));
	dvu_dhr  = -0.125.*vu.*Sux./(Su.^2.*dx(1));
	dvu_dhul =  0.125.*vu.*Sux./(Su.^2.*dx(1));
	dvu_dhur = -0.125.*vu.*Sux./(Su.^2.*dx(1));

	nn = prod(n);
	J  = zeros(nn,9);
	if (issym(h))
		J = sym(J);
	end
	wc = 1;
	wlr = 1;
	wud = 1;
	% up
	J(:,1) = wud*flat(+0.5*dvu_dhu.*( (1-bu).*hu  + (1+bu).*h)/dx(2) ...
		      +0.5*vu.*(1-bu)/dx(2) ...
		      -0.5*dvr_dhu.*( (1-br).*h  + (1+br).*hr)/dx(1) ... 
		      +0.5*dvl_dhu.*( (1-bl).*hl + (1+bl).*h)/dx(1));
	out.hu = J(3+2,1);

	% up right
	J(:,2) = wc*flat(-0.5*dvr_dhur.* ((1-br).*h   + (1+br).*hr)/dx(1) ...
		      +0.5*dvu_dhur.*( (1-bu).*hu + (1+bu).*h)/dx(2));
	out.hur = J(3+2,2);

%
%	dh_dt = -( ( ...
%		       0.5*right(ul,1).*( (1-br).*h         + (1+br).*right(h,1)) ...
%	             - 0.5*      ul   .*( (1-bl).*left(h,1) + (1+bl).*h) ...
%		   )/dx(1) ...
%		 + ( ...
%		       0.5*down(uu,1).*( (1-bd).*h          + (1+bd).*down(h,1)) ...
%	             - 0.5*     uu   .*( (1-bu).*up(h,1)    + (1+bu).*h) ...
%		   )/dx(2) ...
%		 );

	% right dhr
	J(:,3)   = wlr*flat(-0.5*dvr_dhr.*( (1-br).*h  + (1+br).*hr)/dx(1) ...
		        -0.5*vr.*(1+br)/dx(1) ...
		        -0.5*dvd_dhr.*( (1-bd).*h  + (1+bd).*hd)/dx(2) ...
		        +0.5*dvu_dhr.*( (1-bu).*hu + (1+bu).*h )/dx(2));
	out.hr = J(3+2,3);

	% dhdr	
	J(:,4) = wc*flat(-0.5*dvr_dhdr.*( (1-br).*h + (1+br).*hr)/dx(1) ...
		      -0.5*dvd_dhdr.*( (1-bd).*h + (1+bd).*hd)/dx(2));
	out.hdr = J(3+2,4);
	% down
	J(:,5) = wud*flat(-0.5*dvd_dhd.*( (1-bd).*h  + (1+bd).*hd)/dx(2) ...
		      -0.5*vd.*(1+bd)/dx(2) ...
		      -0.5*dvr_dhd.*( (1-br).*h  + (1+br).*hr)/dx(1) ...
		      +0.5*dvl_dhd.*( (1-bl).*hl + (1+bl).*h)/dx(1) );
	out.hd = J(3+2,5);

	% down left
	J(:,6)   = wc*flat(+0.5*dvl_dhdl.*((1-bl).*hl + (1+bl).*h)/dx(1) ...
		        -0.5*dvd_dhdl.*((1-bd).*h  + (1+bd).*hd)/dx(2));
	out.hdl = J(3+2,6);

	% left
	J(:,7)    = wlr*flat(+ 0.5*vl.*(1-bl)/dx(1) ...
			 +0.5*dvl_dhl.*( (1-bl).*hl + (1+bl).*h)/dx(1) ...
		         -0.5*dvd_dhl.*( (1-bd).*h  + (1+bd).*hd)/dx(2) ...
	                 +0.5*dvu_dhl.*( (1-bu).*hu + (1+bu).*h)/dx(2));
	out.hl = J(3+2,7);

	% top left
	J(:,8)   = wc*flat(+0.5*dvl_dhul.*((1-bl).*hl + (1+bl).*h)/dx(1) ...
		        +0.5*dvu_dhul.*((1-bu).*hu + (1+bu).*h)/dx(2));
	out.hul = J(3+2,8);

	% u ur r dr d dl l ul
	J(:,9) = flat( 	-  down(      reshape(J(:,1),n),1) ...
			-  down( left(reshape(J(:,2),n),1),1) ...
			-  left(      reshape(J(:,3),n),1) ...
			-    up( left(reshape(J(:,4),n),1),1) ...
			-    up(      reshape(J(:,5),n),1) ...
			-    up(right(reshape(J(:,6),n),1),1) ...
			- right(      reshape(J(:,7),n),1) ...
			-  down(right(reshape(J(:,8),n),1),1) ...
		     );
	% J(:,9) = -sum(J,2);

	if (returnmat)
		J = diags2mat(J,n);
	end

	out.Sl = Sl;
	out.Slx = Slx;
	out.Sly = Sly;
	out.Su = Su;
	out.Sux = Sux;
	out.Suy = Suy;
	out.vu = vu;
	out.vl = vl;
	out.bl = bl;
	out.bu = bu;
	out.dvd_dhd = dvd_dhd(2,2);
	out.dvl_dhl = dvl_dhl(2,2);
	out.dvr_dhr = dvr_dhr(2,2);
	out.dvu_dhu = dvu_dhu(2,2);
	out.dvd_dhdl = dvd_dhdl(2,2);
	out.dvd_dhdr = dvd_dhdr(2,2);
	out.dvd_dhl = dvd_dhl(2,2);
	out.dvd_dhr = dvd_dhr(2,2);
	out.dvl_dhd = dvl_dhd(2,2);
	out.dvl_dhdl = dvl_dhdl(2,2);
	out.dvl_dhu = dvl_dhu(2,2);
	out.dvl_dhul = dvl_dhul(2,2);
	out.dvr_dhd = dvr_dhd(2,2);
	out.dvr_dhdr = dvr_dhdr(2,2);
	out.dvr_dhu = dvr_dhu(2,2);
	out.dvr_dhur = dvr_dhur(2,2);
	out.dvu_dhl = dvu_dhl(2,2);
	out.dvu_dhr = dvu_dhr(2,2);
	out.dvu_dhul = dvu_dhul(2,2);
	out.dvu_dhur = dvu_dhur(2,2);
end % swe_zero_inertia_2d_jacobian

