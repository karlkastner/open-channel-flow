% Mon 19 May 15:19:26 CEST 2025
	colorbar
% Fri  4 Jul 14:38:57 CEST 2025

% test in 1D
% 2D left
% 2D right

unit = 'm'
%tscale = 1e-2; % stretch 2 hour events over 7 day return interval, so reduce rate of change
tscale =1;
switch(unit)
case {'mm'}
	%factor = 0.001;
case {'m'}
	%factor = 1;
end
input_factor = 1e-3;
input_factor = 0.1;
tscale = 10;

L    = 9;
dx   = 1;
n    = L/dx;
% water depth
hmax_m = 1;
hmax = hmax_m; %/factor;
hmin = 0.01*hmax;
hmin = 1;
sd   = 0.1*L;

% mean bed slope
dzb_dx = 1e-3; %1e-2*[-1,0,+1];
dzb_dy = 1e-3; %1e-2*[-1,0,+1];
abstol = 1e-7;

% chezy coefficient
C   = 10;
%C = 15;
lcd = 1e-0;
%lcd = 0;
% end time % seconds
T = 100/tscale;
dt = 10/tscale;
dte = 1/tscale;

% 1 - dt*K*h*2/dx^2 > 0

%discretization_C = {'central','upwind','optimal','optimal-linearized'};
%discretization_C = {'upwind'} % ,'optimal','optimal-linearized'};
discretization_C = {'central'}
solver_C = {'explicit'} %,'implicit'}
%solver_C = {'implicit'}
linear_solver = 'ilu-ks';
linear_solver = 'direct';
maxiter = 200;

% spatial axis
x  = (0:n-1)'*L/(n-1);
xe = (-1:n)'*L/(n-1);

% spatial step
dx = x(2)-x(1);

% initial water level, a gaussian bump
mu = L/2;
p = 1
% h0 = hmin + hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
h0   = hmin + 0.5*hmax*(1 + cos(2*pi*x/L));
%h0 = hmax*ones(n,1);
%h0 = hmin + hmax*(h0>0.5*hmax);
%h0 = circshift(h0,-0.4*n);

% mean coordinate, depth averaged
xc0 = sum(h0.*x)./sum(h0);

% for each bed slope
for zdx=1:length(dzb_dx)

% for each solver
for kdx=1:length(solver_C)

% for each discretization
for jdx=1:length(discretization_C)
	% bed level
	zb   = dzb_dx(zdx)*(xe-L/2);
	zb_y   = dzb_dy(zdx)*(xe-L/2);
	h0   = flat(h0);
	returnmat = true;
	
	zi = SWE_Zero_Inertia_1d();
	zi.discretization = discretization_C{jdx};
	zi.returnmat = returnmat;
	zi.zb = zb;
	zi.C  = C;
	zi.L  = L;
	zi.n  = n;
	zi.gn_abstol = abstol;
	zi.gn_maxiter = maxiter;
	zi.linear_solver_name = linear_solver;
	zi.linear_solver_tol = abstol;
	zi.linear_solver_maxiter = maxiter;
	zi.integration_scheme = 'explicit'; %trapezoidal';
	zi.lcd = lcd;
	zi.input_factor = input_factor;
	zi.output_factor = tscale;
	dh_dt_1d = cvec(zi.dh_dt(0,h0));
	J0 = zi.jacobian(h0(:));
	J0 = reshape(J0,n,n);
	

	zb_2d = cvec(zb) + rvec(zb_y);%repmat(cvec(zb),1,n+2);
	h_2d  = repmat(cvec(h0),1,n);
	scale = 1/max(abs(dh_dt_1d));
	
	% for each dimension
	for ddx = 1
	if (3 == ddx)
		h0 = hmin ...
			+ hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p)) ...
			     .*(1-(1-(normpdf(x',mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
		h0   = hmin + 0.25*hmax*(1 + cos(2*pi*x/L)) + 0.25*hmax*(1 + cos(2*pi*x'/L));
		h_2d = h0;
	end

	zi2d = SWE_Zero_Inertia_2d();
	zi2d.discretization = discretization_C{jdx};
	zi2d.returnmat = returnmat;
	zi2d.zb = zb_2d;
	zi2d.C  = C;
	zi2d.lcd = lcd;
	zi2d.L  = L*[1,1];
	zi2d.n  = n*[1,1];
	zi2d.gn_abstol  = abstol;
	zi2d.gn_maxiter = maxiter;
	zi2d.linear_solver_name = linear_solver;
	zi2d.linear_solver_tol = abstol;
	zi2d.linear_solver_maxiter = maxiter;
	zi2d.integration_scheme = 'trapezoidal';
	zi2d.lcd = lcd;
	zi2d.input_factor = input_factor;
	zi2d.output_factor = tscale;

	[qi,qj] = zi2d.interface_values(h_2d);

	figure(1e3);
	subplot(3,3,ddx)
	imagesc(h_2d)
	title('h');
	subplot(3,3,3+ddx)
	imagesc(qi)
	title('qx');
	subplot(3,3,6+ddx)
	imagesc(qj)
	title('qy');


	dh_dt = zi2d.dh_dt(0,h_2d);
	dh_dt = reshape(dh_dt,n*[1,1]);
	if (1==ddx)
		res = dh_dt - cvec(dh_dt_1d);
	else
		res = dh_dt - rvec(dh_dt_1d);
	end
	max(abs(res(:)))/max(abs(dh_dt_1d))
	%(res(:))/rms(dh_dt)

	zi2d.analytic_jacobian = false;
	zi2d.eps = 1e-6;
	zi2d.twosided = true;	
	J = zi2d.jacobian(h_2d);

	% compute the jacobian
	J_ = zeros(n*n,n*n);
	e_ = 1e-2;
	%e = 1e-2*max(h_2d(:));
	e = 1e-7*max(h_2d(:));
	for idx=1:(n*n)
		%e = e_*h_2d(idx);
		%e  = e_*h_2d(idx);
		hl = h_2d;
		hl(idx) = hl(idx)-e;
		hr = h_2d;
		hr(idx) = hr(idx)+e;
		%d = (zi2d.dh_dt(0,hr) - zi2d.dh_dt(0,hl))/(2*e); 
		dhr = zi2d.dh_dt(0,hr);
 		dhl = zi2d.dh_dt(0,hl);
		d = (dhr - dhl)/(2*e);
if (idx == round(n*(n-1)/2+n/2)) 
figure(2000)
if (0)
subplot(3,4,1+4*(ddx-1))
%imagesc(reshape(d,n,n)); axis square
imagesc(reshape(dhl,n,n))
colorbar
subplot(3,4,2+4*(ddx-1))
imagesc(reshape(dhr,n,n))
%imagesc(reshape(J(:,idx),n,n)); axis square
colorbar
end
subplot(3,3,1+3*(ddx-1))
imagesc(reshape(d,n,n))
%imagesc(reshape(J(:,idx)-d,n,n)); axis square
colorbar
title('numerical');
subplot(3,3,2+3*(ddx-1))
imagesc(reshape(J(:,idx),n,n));
colorbar
subplot(3,3,3+3*(ddx-1))
imagesc(reshape(J(:,idx)-d(:),n,n));
colorbar
title('Analytic');
%idx
%pause(1)
end
		J_(:,idx) = d; 
	end % for idx
	%J_ = J_';

	id = reshape((1:(n*n)^2),n*n,n*n);
	for i1 = 1:3
	for i2 = 1:3
		Jc  = circshift(J, (i2-2) + n*(i1-2),2);
		J_c = circshift(J_,(i2-2) + n*(i1-2),2);
		d   = diag(Jc);
		d_  = diag(J_c);
		d = reshape(d,n,n);
		d_ = reshape(d_,n,n);
		J0_ = diag(circshift(J0,i2-2,2));
	figure(100)
	subplot(3,3,ddx)
	%spy(J)
	imagesc(J)
	axis square
	title('analytic')
	colorbar
	subplot(3,3,3+ddx)
	%spy(J_)
	imagesc(J_)
	axis square
	title('numerical')
	colorbar
	subplot(3,3,6+ddx)
	%spy(J_)
	imagesc(J-J_)
	axis square
	title('difference')
	colorbar
	figure(ddx);
	subplot(3,9,i2 + 3*(i1-1))
		cla
		%d = reshape(J(diag(id)),n,n);
		imagesc(reshape(d,n,n));
		colorbar
	%	subplot(3,9,9+sdx);
	subplot(3,9,i2 + 3*(i1-1)+9)
		cla
		d_ = reshape(J_(diag(id)),n,n);
		imagesc(reshape(d_,n,n));
		colorbar
	%	subplot(3,9,18+sdx);
	subplot(3,9,i2 + 3*(i1-1)+18)
		cla
		imagesc(reshape(d-d_,n,n));
		colorbar
	if (1)
	figure(100+ddx);
		%subplot(3,9,sdx);
		subplot(3,9,i2 + 3*(i1-1))
		cla
		if (ddx==1)
			plot(d(:,2));
			hold on
			plot(d_(:,2),':');
		else
			plot(d(2,:));
			hold on
			plot(d_(2,:),':');
		end
		plot(J0_,'--');
		%imagesc(reshape(d,n,n));
		%subplot(3,9,9+sdx);
		%subplot(3,9,i2 + 3*(i1-1)+9)
		%cla
		%if (ddx==1)
		%	plot(d_(:,1));
		%else
		%	plot(d_(1,:));
		%end
		%hold on
		%plot(J0_);
	end
	end % i2
	end % i1

	figure(1e3+ddx)
	%J=zi2d.jacobian(h_2d,true);
	mJ = max(max(abs(J)));
 Ji = zi2d.jacobian_i(h_2d);
 Jj = zi2d.jacobian_j(h_2d);
 subplot(2,2,1);
 imagesc(Ji);
 caxis(mJ*[-1,1]);
 colorbar;
 subplot(2,2,2);
 imagesc(Jj);
caxis(mJ*[-1,1]);
 colorbar;
 subplot(2,3,4);
 Jij = (Ji+Jj);
 imagesc(Jij);
 caxis(mJ*[-1,1]);
colorbar;
 subplot(2,3,5);
 imagesc(J);
 colorbar;
 caxis(mJ*[-1,1]);
 subplot(2,3,6);
 imagesc(J-Jij);
 colorbar;
 caxis(mJ*[-1,1]);




	zb_2d = zb_2d';
	h_2d  = h_2d';

	

	end % for ddx each dimension
end % jdx
end % kdx
end % zdx

