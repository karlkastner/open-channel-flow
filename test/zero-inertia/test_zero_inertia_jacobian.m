% Fri 25 Apr 19:02:15 CEST 2025

ifactor = 0.001;
ofactor = 1; % for time
ofactor = 1e-2;
zero_inertia = SWE_Zero_Inertia_1d();
%rad.aux.zero_inertia;
dx=1;
L  = 100;
n  = L/dx;
nx = n;
x = (0:n-1)'*L/n;
xe = (-1:n)'*L/n;

zb = zeros(n+2,1);
zb = 1e-2*xe;
zero_inertia.zb = zb;
zero_inertia.L  = L;
zero_inertia.n  = n;
zero_inertia.lcd = 1e-2;
zero_inertia.C   = 0.7;
zero_inertia.returnmat   = true;
zero_inertia.input_factor = ifactor;
zero_inertia.output_factor = ofactor;


zero_inertia.discretization = 'central';
h = 1;

if (0)
zero_inertia.boundary_condition = {[0,0,1,0];
                          [0,0,1,0]};
else
zero_inertia.boundary_condition = {[0,1,0,0];
                          [0,1,0,0]};
end

% 2025-04-22 17:05:16.998366415 +0200
close all;
y0 = (1+0.5*sin(2*pi*x/L)).*h;
%mean(hh(:,1));
%y0 = mean(hh(:,1))+x.^2/param.L.^2;
%y0 = mean(hh(:,1)) + 0.5*(1+sin(2*pi*x/param.L));
%y0 = mean(hh(:,1)) + x/param.L;
%y0 = h*(1 + 0.1*x/L);
%0.5*(1+sin(2*pi*x/param.L));
%.*mean(hh(:,1));


s = 1;

y=sin(4*pi*x/L);
e=1e-7;
J = zero_inertia.jacobian(y0);
Jy = s*J*y;
Jy(:,2) = (zero_inertia.dh_dt(y0+e/2*y) - zero_inertia.dh_dt(y0-e/2*y))/e;

J_ = [];
dh_dt0 = zero_inertia.dh_dt(y0);
for idx=1:nx
	ei = zeros(nx,1);
	ei(idx)=1;
	%J_(:,idx) = (zero_inertia.dh_dt(y0+e*ei)-dh_dt0)/e;
	J_(:,idx) = (zero_inertia.dh_dt(y0+0.5*e*ei) - zero_inertia.dh_dt(y0-0.5*e*ei))/e;
end
%J_ = J_';


clf;
subplot(2,2,1)
plot(x,Jy(:,1));
hold on;
plot(x,Jy(:,2))     
plot(x,J_*y);
%plot(Jx(:,1),Jx(:,2),'.');
%median(Jx(:,2)./Jx(:,1))
legend('analytic','numerical')
title('Jy')

subplot(2,3,4)
plot([s*diag(J,-1),diag(J_,-1)])
legend('analytic','numerical');
title('J_-1')

subplot(2,3,5);
plot([s*[diag(J,0)],[diag(J_,0)]])
title('J_0')

subplot(2,3,6);
plot([[s*diag(J,+1)],diag(J_,+1)])
title('J_+1')

% [hi,he,ui,bi] = zero_inertia.interface_values(y0);

