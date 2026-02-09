% 2025-05-15 19:04:11.173966755 +0200
% test in 1D
% 2D left
% 2D right

unit = 'm'
tscale = 1e-2; % stretch 2 hour events over 7 day return interval, so reduce rate of change

switch(unit)
case {'mm'}
	factor = 0.001;
case {'m'}
	factor = 1;
end

L    = 2*99;
%n    = 33;
n = L;
% water depth
hmax_m = 1e-2;
hmax = hmax_m/factor;
hmin = 0.01*hmax;
sd    = L/10;

% mean bed slope
dzb_dx = 1e-2; % 1e-2*[-1,0,+1];
dzb_dx = 5e-3;
abstol = 1e-7;

% chezy coefficient
C   = 1;
%C = 15;
lcd = 1e-1;
%lcd = 0;
% end time % seconds
T = 100/tscale;
dt = 10/tscale;
dte = 1/tscale;

% 1 - dt*K*h*2/dx^2 > 0

%discretization_C = {'central','upwind','optimal','optimal-linearized'};
discretization_C = {'upwind'} % ,'optimal','optimal-linearized'};
%discretization_C = {'central'}
discretization_C = {'optimal'}
solver_C = {'explicit'} %,'implicit'}
%solver_C = {'implicit'}
linear_solver = 'direct';
integration_scheme = 'aid';
%linear_solver = 'ilu-ks';
%integration_scheme = 'trapezoidal';
maxiter = 200;

% spatial axis
x  = (1/2+(0:n-1)')*L/(n-1);
xe = (-1:n)'*L/(n-1);

% spatial step
dx = x(2)-x(1);

% initial water level, a gaussian bump
mu = L/2;
mu = mean(x);
p = 1
h0 = hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
%h0 = hmin + 
%h0 = hmin + hmax*(h0>0.5*hmax);
%h0 = circshift(h0,-0.4*n);

% mean coordinate, depth averaged
xc0 = sum(h0.*x)./sum(h0);

% for each bed slope
for zdx=1:length(dzb_dx)

for kdx=1:length(solver_C)

for jdx=1:length(discretization_C)
	% bed level
	zb   = dzb_dx(zdx)*(xe-L/2);
	h0   = flat(h0);
	h_1d = hmin + h0;
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
	zi.linear_solver_name = 'direct'; %linear_solver;
	zi.linear_solver_tol = abstol;
	zi.linear_solver_maxiter = maxiter;
	zi.integration_scheme = 'explicit';
	zi.integration_scheme = 'trapezoidal';
	zi.lcd = lcd;
	zi.input_factor = factor;
	zi.output_factor = tscale;
	dh_dt_1d = cvec(zi.dh_dt(0,h_1d));
	s = 4;
	%T = [0:(s*20)]*1e6/s;
	%T = [0:(s*20)]*1e6/s;
	%T = [0:(s*5)]*1e6/s*4;
	%T = [0:(s*3)]*1e6/s*7;
	T = [0:(s*4)]*1e6/s*5;
	hh_1d = zi.solve(h_1d,T);

	zb_2d = repmat(cvec(zb),1,n+2);
	h_2d  = repmat(cvec(h0),1,n)+hmin;
	scale = 1/max(abs(dh_dt_1d));

	figure(10);
	subplot(2,4,1)
	cla
	plot(x-mu,h0);
	hold on

	figure(1)
	subplot(2,4,1);
	cla
	plot(x-mu,scale*dh_dt_1d);
	hold on

	figure(3);
	subplot(2,4,1)
	cla
	plot(x-mu,hh_1d(:,end));
	hold on

	for ddx = 1:3
	printf('ddx %g\n',ddx);
	zi2d = SWE_Zero_Inertia_2d();
	zi2d.discretization = discretization_C{jdx};
	zi2d.returnmat = returnmat;
	zi2d.C  = C;
	zi2d.L  = L*[1,1];
	zi2d.n  = n*[1,1];
	zi2d.gn_abstol  = abstol;
	zi2d.gn_maxiter = maxiter;
	zi2d.linear_solver_name = linear_solver;
	zi2d.linear_solver_tol = abstol;
	zi2d.linear_solver_maxiter = maxiter;
	zi2d.analytic_jacobian  = false;
	%zi2d.integration_scheme = 'explicit'; %trapezoidal';
	%zi2d.integration_scheme = 'trapezoidal';
	zi2d.integration_scheme = integration_scheme;



	if (3 == ddx)
		zb_2d_   = (zb + zb')/sqrt(2);
		xx    = repmat(cvec(x),1);
		yy    = xx';
		s     = 1/sqrt(2);
		h_2d_0 = hmin + hmax*(1-(1-(normpdf(s*(xx-mu)+s*(yy-mu),0,sd)/normpdf(0,0,sd)).^p).^(1/p));
		%h_2d_ = imrotate(h_2d_,-45,'crop','bilinear');
	else
		zb_2d_ = zb_2d;
		h_2d_0 = h_2d;
	end	

	zi2d.zb = zb_2d_;
	zi2d.lcd = lcd;
	zi2d.input_factor = factor;
	zi2d.output_factor = tscale;
	dh_dt = zi2d.dh_dt(0,h_2d_0);
	dh_dt = reshape(dh_dt,n*[1,1]);

	if (0)
	rms(zi2d.dh_dt_i(0,h_2d_)+zi2d.dh_dt_j(0,h_2d_)-zi2d.dh_dt(0,h_2d_))/rms(zi2d.dh_dt(0,h_2d_))
	
	% test the jacobian
	Ji  = jacobian_unstructured(@(h) zi2d.dh_dt_i(0,h),h_2d_(:));
	Ji_ = zi2d.jacobian_i(0,h_2d_);
	Jj  = jacobian_unstructured(@(h) zi2d.dh_dt_j(0,h),h_2d_(:));
	Jj_ = zi2d.jacobian_j(0,h_2d_);

	figure(1e5+ddx);
	subplot(2,3,1)
	imagesc(Ji);
	colorbar
	subplot(2,3,2)
	imagesc(Ji_);
	colorbar
	subplot(2,3,3)
	imagesc(Ji_-Ji);
	colorbar
	subplot(2,3,4)
	imagesc(Jj);
	colorbar
	subplot(2,3,5)
	imagesc(Jj_);
	colorbar
	subplot(2,3,6)
	imagesc(Jj_-Jj);
	colorbar
	figure(1e6+ddx);
	subplot(2,3,1)
	imagesc(reshape(diag(Ji),n,n));
	colorbar
	subplot(2,3,2)
	imagesc(reshape(diag(Ji_),n,n));
	colorbar
	subplot(2,3,3)
	imagesc(reshape(diag(Ji_-Ji),n,n));
	colorbar
	subplot(2,3,4)
	imagesc(reshape(diag(Jj),n,n));
	colorbar
	subplot(2,3,5)
	imagesc(reshape(diag(Jj_),n,n));
	colorbar
	subplot(2,3,6)
	imagesc(reshape(diag(Jj_-Jj),n,n));
	colorbar
	end

	
	% [hh,exitflag,rmsr,rmsg,res,iter,out] = solve(obj,h0,T) 
	tic()
	[hh_2d, exitflag, rmsr, rmsg, res, iter] = zi2d.solve(h_2d_0(:),T);
	rt(ddx) = toc()
	figure(1e4)
	subplot(2,3,ddx)
	plot(iter)	
	switch (ddx)
	case {1}
		res = dh_dt - cvec(dh_dt_1d);
		e_rel = max(abs(hh_2d(1:n,:)-hh_1d),[],'all')/max(h0)
	case{2}
		res = dh_dt - rvec(dh_dt_1d);
		e_rel = max(abs(hh_2d(1:n:n^2,:)-hh_1d),[],'all')/max(h0)
	otherwise
		res = dh_dt - cvec(dh_dt_1d);
		e_rel = max(abs(hh_2d(1:n,:)-hh_1d),[],'all')/max(h0)
		
	end
	max(abs(res(:)))/max(abs(dh_dt_1d))
		%(res(:))/rms(dh_dt)
	figure(10);
	subplot(2,4,1)
	switch (ddx)
	case {1}
		plot(x-mu,scale*h_2d_0(:,1));
	case {2}
		plot(x-mu,scale*h_2d_0(1,:));
	case {3}
		%plot((x-mu)*sqrt(2),scale*diag(h_2d_));
		plot(sqrt(2)*(x-mu),scale*diag(h_2d_0));
	end
	title('Initial condition h(0)')
	subplot(2,4,1+ddx);
	imagesc(h_2d_0)
	axis equal

	figure(1);
	subplot(2,4,1)
	switch (ddx)
	case {1}
		plot(x-mu,scale*dh_dt(:,1));
	case {2}
		plot(x-mu,scale*dh_dt(1,:));
	case {3}
		dh_dt_ = reshape(dh_dt,n,n);
		plot((x-mu)*sqrt(2),scale*diag(dh_dt_));
	end

	subplot(2,4,1+ddx);
	imagesc(dh_dt);
	caxis([-1,1]*max(abs(dh_dt_1d)));

	subplot(2,4,5+ddx)
	imagesc(res);
	caxis([-1,1]*max(abs(dh_dt_1d)));

	figure(3);
	subplot(2,4,1)
	h_ = reshape(hh_2d(:,end),n,n);
	switch (ddx)
	case {1}
		plot(x-mu,scale*h_(:,1)); %dh_dt(:,1));
	case {2}
		plot(x-mu,scale*h_(1,:)); %dh_dt(1,:));
	case {3}
		%dh_dt_ = reshape(dh_dt,n,n);
		plot(sqrt(2)*(x-mu),scale*diag(h_));
		%h_ = imrotate(h_,-45,'crop','bilinear');
		%plot((x-mu),scale*h_(:,1),'--'); %diag(h_));
		%dh_dt_));
	end
	legend('x','y','diagona');
	subplot(2,4,1+ddx);
	imagesc(reshape(hh_2d(:,end),n,n))
	axis equal

	zb_2d = zb_2d';
	h_2d  = h_2d';

	h_C{ddx} = hh_2d;
	end % for ddx each dimension
end % jdx
end % kdx
end % zdx

zb = (zb + zb')/sqrt(2);
zi2d.zb = zb;

figure(2)
clf
h_2d = hmin + cvec(h0)*rvec(h0)/max(h0);
dh_dt = zi2d.dh_dt(0,h_2d);
dh_dt = reshape(dh_dt,n,n);
hh_2d = zi2d.solve(h_2d(:),T);

subplot(2,2,1)
imagesc(dh_dt);
title('dh/dt|_(t=0)');

subplot(2,2,2)
surface(reshape(hh_2d(:,end),n,n),'edgecolor','none');
view(-30,45)
title('h(T)');

subplot(2,2,3)
x = (1:n)-0.5;
x = x-mean(x);
plot(x,hh_1d(:,[1,end]))
d = diag(reshape(hh_2d(:,end),n,n));
d(:,2) = h_C{1}(1:n,end);
d(:,3) = h_C{2}(1:n:n^2,end);
%,':')
%,'--')
hold on
plot(x,d(:,2:end),'--')
plot(sqrt(2)*x,d(:,1))
max(hh_2d(:,end),[],'all')


