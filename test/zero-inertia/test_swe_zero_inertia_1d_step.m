% Fri  4 Oct 14:15:48 CEST 2024
% explicit solver works fine, so the problem is with the jacobian
%
% TODO time steps for zero interia has to be lower than for linearized flow
% 


unit = 'm'
tscale = 1e-2; % stretch 2 hour events over 7 day return interval, so reduce rate of change

switch(unit)
case{'mm'}
	factor = 0.001;
case {'m'}
	factor = 1;
end
%factor = 1;

L    = 100;
n    = L*10;
% 10mm water depth
hmax = 1e-2/factor;
hmin = 0.01*hmax;
sd    = 10;


% mean bed slope
dzb_dx = 1e-2*[-1,0,+1];
abstol  = 1e-7;

% chezy coefficient
C = 10;
%C = 15;
lcd = 1e-6; %sqrt(eps);
%lcd = 0;
% end time % seconds
T = 100/tscale;
dt = 10/tscale;
dte = 1/tscale;

% 1 - dt*K*h*2/dx^2 > 0

discretization_C = {'central','upwind','optimal','optimal-linearized'};
solver_C = {'explicit','implicit'}
solver_C = {'implicit'}
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

% mean coordinate, depth averaged
xc0 = sum(h0.*x)./sum(h0);

% for each bed slope
for zdx=1:length(dzb_dx)
%z = z.*cvec(tukeywin(n,0.1));

% bed level
%zb = -dzb_dx(zdx)*(L/(2*pi))*sin(2*pi*xe/L);
zb  = dzb_dx(zdx)*(xe-L/2);
%zb = zb + 0*dzb_dx(zdx)*randn(n,1);
%z = circshift(z,-10);
%h0 = hmax*(h0>0.5*hmax);
%h0 = meanfilt1(h0,3);

for kdx=1:length(solver_C)
figure(kdx+(zdx-1)*10)
for jdx=1:length(discretization_C)
	zb = flat(zb);
	h0 = flat(h0);
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
%	zi.integration_scheme = 'trapezoidal';
	zi.integration_scheme = 'midpoint';
	zi.lcd = lcd;
	zi.input_factor = factor;
	zi.output_factor = tscale;

switch (solver_C{kdx})
case{'explicit'}
	%dt_max = 0.1*min(dx^2./(2*C*h0))
	%nt     = ceil(T/dt_max)
	%nt     = 100;
	%dte     = T/nt;
	nt = round(T/dte);
	h = h0;
	for idx=1:nt
		%[dh_dt,bl,Pel,u,d,S] = swe_zero_inertia_1d_dhdt(h,z,C,L,n,discretization_C{jdx});
		dh_dt = zi.dh_dt(h); %h0);
		h = h + dte*dh_dt;
	end % for idx
	subplot(3,4,4+jdx);
	cla
	plot(x,[h0,h])
	%plot(x,bl)
	hold on
	%yyaxis right
	%plot(x,Pel/3)
	%ylim([-1,1]*1.05)
	%end % if 1 == idx
	
case{'implicit'}
%	zi.integration_scheme = 'euler-implicit';

	
	t = 0:dt:T;
	t = [0,T];
	%[h,exitflag,rmsr,rmsg,res,iter,out] = zi.step(dt,h0);
	%[h,exitflag,rmsr,rmsg,res,iter,out] = zi.step(dt,h0);
	[h] = zi.solve(h0,t);
	h = h(:,end);
%	[h,rmsr,res,iter,out] = swe_zero_inertia_1d_step(dt,h0,zb,C,L,n,discretization_C{jdx},linear_solver,abstol,maxiter,returnmat);

	subplot(3,4,jdx+4);
	%semilogy([cvec(rmsr),cvec(rmsg)],'.-')
	xlabel('GN-Step');
	ylabel('rms(res)');
	%title(exitflag)

	[hi,he,ui,bi] = zi.interface_values(h0,true);
%hi,he,ui,bi,Si,ai,di,Pei
	dh_dt = zi.dh_dt(h0);
	subplot(3,4,jdx+8);
	plot(dh_dt);
	yyaxis right
	plot(bi);
	
end % switch solver

% depth averaged mid coordinate after time step
xc(jdx) = sum(x.*h)./sum(h);

subplot(3,4,jdx);
cla();
%yyaxis left
plot(x,[h0,h]);
xlabel('Position x');
ylabel(['Water depth h/', unit]);
lh=legend('0','T');
title(lh,'Time');
%yyaxis right
%plot(x,z)
title([solver_C{kdx}, ' ', discretization_C{jdx}, ' dzb/dx=', num2str(dzb_dx(zdx))])



hh(:,jdx,zdx) = h;
end % for jdx
dxc=xc-xc0
[max(abs(diff(h0))),max(abs(diff(zb))),dzb_dx(zdx)]
end % for kdx
end % for zdx

if (0)
figure(1e3)
for idx=1:4;
 subplot(2,2,idx);
 plot(real([hh(:,idx,1)-flipud(hh(:,idx,3))]));
 end
end
