% Sun  6 Oct 10:08:52 CEST 2024
L    = 100;
n    = L;
hmax = 0.01; %1e-2;
sd    = 5;
dz_dx_ = 1e-2*[-1,0,1];
for zdx=1:3
dz_dx = dz_dx_(zdx);
x  = (0:n-1)'*L/n;
%z = z.*cvec(tukeywin(n,0.1));
z = -dz_dx*(L/(2*pi))*sin(2*pi*x/L);
z  = dz_dx*(x-L/2);
%z = z + 0*dz_dx*randn(n,1);
%z = circshift(z,-10);
mu = L/2;
p = 1;
hmin = 0.01*hmax;
h0_ = hmin + hmax*(1./(1+((x-mu)/sd).^2));
h0 = hmin+hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
%h0 = hmax*(h0>0.5*hmax);
%h0 = meanfilt1(h0,3);
xc0 = sum(h0.*x)./sum(h0);
	
dx = x(2)-x(1);

K = 1.2e3;
%C = 1e3;
dt = 1;

% 1 - dt*K*h*2/dx^2 > 0

%method_C = {'central','upwind','optimal','optimal-linearized'};
%method_C = {'upwind'}
method_C = {'central'}
%,'optimal-linearized'};
%method_C ={'optimal'}

for jdx=1:length(method_C)
%dt_max = 0.1*min(dx^2./(2*C*h0))
%nt     = ceil(T/dt_max)
%dt     = T/nt;

	h = h0;
	returnmat = true;
	m = method_C{jdx};

	J      = euler_implicit_jacobian(@(h) swe_simplified_1d_jacobian(h,z,K,L,n,m,returnmat), dt, h, returnmat);
	resfun = @(h) euler_implicit_residual(@(h) swe_simplified_1d_dhdt(h,z,K,L,n,method_C{jdx}), dt, h, h0);

	J_ = zeros(n);
	for idx=2:n-1
		e = sqrt(eps)*h0(idx);
		hl = h; hl(idx) = hl(idx)-e;
		hr = h; hr(idx) = hr(idx)+e;
		d = (resfun(hr) - resfun(hl))/(2*e); 
		J_(idx,idx) = d(idx);
		J_(idx-1,idx) = d(idx-1);
		J_(idx+1,idx) = d(idx+1);
	end

yl = 25;
figure(1);
subplot(3,3,zdx);%kdx+(zdx-1)*10)
plot([diag(J),diag(J_)]-1)
ylim([-1,1]*yl)
title(['diagonal ',num2str(dz_dx)]);

subplot(3,3,zdx+3);%kdx+(zdx-1)*10)
plot([diag(J,-1),diag(J_,-1)])
ylim([-1,1]*yl)
title(['left ',num2str(dz_dx)]);

subplot(3,3,zdx+6);%kdx+(zdx-1)*10)
plot([diag(J,+1),diag(J_,+1)])
ylim([-1,1]*yl)
title(['right ',num2str(dz_dx)]);

end % for jdx
end % for zdx
