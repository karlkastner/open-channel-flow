% Mon  7 Oct 10:43:01 CEST 2024
% n.b.: explizit solver still unstable when dz/dx large (maybe bc of time step limitation by advection)
% RQ : nonlinear (solved by me?)
% RQ : crosswind diffusion
% RQ : role of source terms

L    = 100;
n    = L;
hmax = 0.01; %1e-2;
sd   = 5;
dz_dx_ = 1e-2*[-1,1e-5,1];

K = 1.2e3;
T = 0.1;

for zdx=1:length(dz_dx_)
dz_dx = dz_dx_(zdx);
x  = (0:n-1)'*L/n;
%z = z.*cvec(tukeywin(n,0.1));
z = -dz_dx*(L/(2*pi))*sin(2*pi*x/L);
z  = dz_dx*(x-L/2);
z = 0.*z + 0*dz_dx*randn(n,1);
%z = circshift(z,-10);
mu = L/2;
p = 1;
hmin = 0.01*hmax;
h0 = hmin + hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
%h0 = hmax*(h0>0.5*hmax);
%h0 = meanfilt1(h0,3);
xc0 = sum(h0.*x)./sum(h0);
	
dx = x(2)-x(1);

% 1 - dt*K*h*2/dx^2 > 0

method_C = {'central','upwind','optimal','optimal-linearized'};
solver_C = {'explicit','implicit'}
for kdx=1:length(solver_C)
figure(kdx+(zdx-1)*10)
for jdx=1:length(method_C)
switch (solver_C{kdx})
case{'explicit'}
dt_max = min(dx^2./(2*K*h0))
nt     = ceil(T/dt_max)
dt     = T/nt;
h = h0;
for idx=1:nt
	[dh_dt,bl,Pel] = swe_simplified_1d_dhdt(h,z,K,L,n,method_C{jdx});
	h = h + dt*dh_dt;
if (1==idx)
subplot(2,4,4+jdx);
cla
plot(x,bl)
hold on
%yyaxis right
plot(x,Pel/3)
ylim([-1,1]*1.05)
end % if 1 == idx
end % for idx

case{'implicit'}
	%method = 'ilu-ks';
	method = 'direct';
	maxiter = 200;
	z = flat(z);
	h0 = flat(h0);
	h = h0;
	dt = T;
	returnmat         = true;
	resfun            = @(h) euler_implicit_residual(@(h) swe_simplified_1d_dhdt(h,z,K,L,n,method_C{jdx}), dt, h, h0);
	m = method_C{jdx};
	%m = 'upwind';
	jfun              = @(h) euler_implicit_jacobian(@(h) swe_simplified_1d_jacobian(h,z,K,L,n,m,returnmat), dt, h, returnmat);
	[h,rmse,res,iter,out] = gauss_newton(resfun,h,[],maxiter,method,jfun);
subplot(2,4,jdx+4);
	semilogy(rmse,'.-')
end % switch solver

xc(jdx) = sum(x.*h)./sum(h);

subplot(2,4,jdx);
cla();
%yyaxis left
plot(x,[h0,h]);
%yyaxis right
%plot(x,z)
title([solver_C{kdx}, ' ', method_C{jdx}, ' ', num2str(dz_dx)])

end % for jdx
dxc=xc-xc0
[max(abs(diff(h0))),max(abs(diff(z))),dz_dx]
end % for kdx
end

