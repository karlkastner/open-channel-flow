% Sun  6 Oct 10:08:52 CEST 2024
L    = 30;
n    = L;
hmax = 0; %1e-2;
hmax = 0;
sd    = 5;
dzb_dx_ = 1e-2*[-1,0,1];
dzb_dx_ = 0;
%dzb_dx_ = -1e-2;
%dzb_dx_ = 0;
method_C = {'central','upwind','optimal','optimal-linearized'};
method_C = {'central'};

for zdx=1:length(dzb_dx_)
dzb_dx = dzb_dx_(zdx);
x  = (0:n-1)'*L/n;
xe  = (-1:n)'*L/n;
%z = z.*cvec(tukeywin(n,0.1));
%zb = -dz_dx*(L/(2*pi))*sin(2*pi*xe/L);
zb  = dzb_dx*(xe-L/2);
%zb = z + 0*dz_dx*randn(n,1);
%z = circshift(z,-10);
mu = L/2;
p = 1;
hmin = 1;
%h0_ = hmin + hmax*(1./(1+((x-mu)/sd).^2));
%h0 = hmin+hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
h0 = hmin+0.5*hmax*(1-cos(2*pi*x/L));
%h0 = hmax*(h0>0.5*hmax);
%h0 = meanfilt1(h0,3);
xc0 = sum(h0.*x)./sum(h0);
	
dx = x(2)-x(1);

%K = 1.2e3;
C = 10;
lcd = 1;
T = 1;

% 1 - dt*K*h*2/dx^2 > 0


zi = SWE_Zero_Inertia_1d;
zi.n = n;
zi.C = C;
zi.lcd = lcd;
zi.L = L;
zi.zb = zb;
zi.p = 1;
for jdx=1:length(method_C)
dt_max = 0.1*min(dx^2./(2*C*h0))
nt     = ceil(T/dt_max)
dt     = T/nt;
	zi.discretization = method_C{jdx};

	h = h0;
	%returnmat = true;
	%m = method_C{jdx};
	J = zi.jacobian(h);
	%J = euler_implicit_jacobian(@(h) swe_zero_inertia_1d_jacobian(h,z,C,L,n,m,returnmat), dt, h, returnmat);

	%@(h) euler_implicit_residual(@(h) swe_zero_inertia_1d_dhdt(h,z,C,L,n,method_C{jdx}), dt, h, h0);

if (1)
	J_ = sparse([],[],[],n,n);
else
	J_ = zeros(n);

	for idx=1:n
		e = sqrt(eps)*h0(idx);
		hl = h;
		hl(idx) = hl(idx)-e;
		hr = h;
		hr(idx) = hr(idx)+e;
	%resfun = @(h) zi.dh_dt(0,h);
		dh_dt_l = zi.dh_dt(0,hl);
		dh_dt_r = zi.dh_dt(0,hr);
		d = (dh_dt_r-dh_dt_l)/(2*e);
	if (idx == n/2)
		figure(1e3)
		subplot(3,3,1+3*(zdx-1))
		plot(dh_dt_l)
		subplot(3,3,2+3*(zdx-1))
		plot(dh_dt_r)
		subplot(3,3,3+3*(zdx-1))
		plot([J(:,idx),d]);
		legend('analytic','numeric')
		title(rms(J(:,idx)-d))
	end
 
if (0)
		J_(idx,idx) = d(idx);
		J_(idx-1,idx) = d(idx-1);
		J_(idx+1,idx) = d(idx+1);
end
		J_(:,idx) = d; 
	end % for idx
end
%J_ = J_';

figure(10)
subplot(1,3,zdx);%kdx+(zdx-1)*10)
plot(x,zi.dh_dt(0,h),'.-');

figure(jdx);
subplot(3,3,zdx);%kdx+(zdx-1)*10)
plot([diag(J),diag(J_)]-0*1)
%ylim([-1,1]*1)
title(['diagonal slope=',num2str(dzb_dx), ' ',method_C{jdx}]);

subplot(3,3,zdx+3);%kdx+(zdx-1)*10)
plot([diag(J,-1),diag(J_,-1)])
%ylim([-1,1]*1)
title(['left slope=',num2str(dzb_dx),' ',method_C{jdx}]);

subplot(3,3,zdx+6);%kdx+(zdx-1)*10)
plot([diag(J,+1),diag(J_,+1)])
%ylim([-1,1]*1)
title(['right slope=',num2str(dzb_dx),' ',method_C{jdx}]);

end % for jdx
end % for zdx

zi.analytic_jacobian=false;
zi.twosided = true;
zi.eps = 1e-6;
Jn = zi.jacobian(h);
 figure(100);
 subplot(2,2,1);
 imagesc(J);
 colorbar
 subplot(2,2,2);
 imagesc(Jn);
 colorbar
 subplot(2,2,3);
 imagesc(Jn-J)
 colorbar

