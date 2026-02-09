
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

L    = 100;
n    = L;
% water depth
hmax_m = 1e-2;
hmax = hmax_m/factor;
hmin = 0.01*hmax;
sd    = 10;

% mean bed slope
dzb_dx = 1e-2; % 1e-2*[-1,0,+1];
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
h0 = hmin + hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
%h0 = hmin + hmax*(h0>0.5*hmax);
h0 = circshift(h0,-0.4*n);

% mean coordinate, depth averaged
xc0 = sum(h0.*x)./sum(h0);

% for each bed slope
for zdx=1:length(dzb_dx)

for kdx=1:length(solver_C)

for jdx=1:length(discretization_C)
	% bed level
	zb   = dzb_dx(zdx)*(xe-L/2);
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
	zi.input_factor = factor;
	zi.output_factor = tscale;
	dh_dt_1d = cvec(zi.dh_dt(0,h0));

	zb_2d = repmat(cvec(zb),1,n+2);
	h_2d  = repmat(cvec(h0),1,n);
	scale = 1/max(abs(dh_dt_1d));
	
	figure(1)
	subplot(2,3,1);
	cla
	plot(scale*dh_dt_1d);
	hold on

	for ddx = 1:2

	zi2d = SWE_Zero_Inertia_2d();
	zi2d.discretization = discretization_C{jdx};
	zi2d.returnmat = returnmat;
	zi2d.zb = zb_2d;
	zi2d.C  = C;
	zi2d.L  = L*[1,1];
	zi2d.n  = n*[1,1];
	zi2d.gn_abstol  = abstol;
	zi2d.gn_maxiter = maxiter;
	zi2d.linear_solver_name = linear_solver;
	zi2d.linear_solver_tol = abstol;
	zi2d.linear_solver_maxiter = maxiter;
	zi2d.integration_scheme = 'trapezoidal';
	zi2d.lcd = lcd;
	zi2d.input_factor = factor;
	zi2d.output_factor = tscale;
	dh_dt = zi2d.dh_dt(0,h_2d);
	dh_dt = reshape(dh_dt,n*[1,1]);
	if (1==ddx)
		res = dh_dt - cvec(dh_dt_1d);
	else
		res = dh_dt - rvec(dh_dt_1d);
	end
	max(abs(res(:)))/max(abs(dh_dt_1d))
	%(res(:))/rms(dh_dt)

	subplot(2,3,1)
	if (1==ddx)
	plot(scale*dh_dt(:,1));
	else
	plot(scale*dh_dt(1,:));
	end

	subplot(2,3,1+ddx);
	imagesc(dh_dt);
	caxis([-1,1]*max(abs(dh_dt_1d)));

	subplot(2,3,4+ddx)
	imagesc(res);
	caxis([-1,1]*max(abs(dh_dt_1d)));

	zb_2d = zb_2d';
	h_2d  = h_2d';

	end % for ddx each dimension
end % jdx
end % kdx
end % zdx

figure(2)
h_2d = cvec(h0)*rvec(h0)/max(h0);
dh_dt = zi2d.dh_dt(0,h_2d);
dh_dt = reshape(dh_dt,n,n);
subplot(2,2,1)
imagesc(dh_dt);


