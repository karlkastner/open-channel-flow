% Sun  6 Oct 10:08:52 CEST 2024
L    = 25*[1,1];
n    = L;
hmax = 0.01; %1e-2;
hmin = 0.01*hmax;
sd    = 5*[1,1];
dz_dx_ = 1e-2*[-1,0,1];
dz_dy_ = 1e-2*[-1,0,1]*0;
dz_dx_ = 0*[1,1];
dy_dx_ = 0*[1,1];
sdx_ = [5,1/sqrt(eps)];
sdy_ = [1e3*L,1/sqrt(eps),5];

for zdx=1:2 % :length(dz_dx_)
sdx = sdx_(zdx);
sdy = sdy_(zdx);
dz_dx = dz_dx_(zdx);
dz_dy = dz_dy_(zdx);
x  = (0:n(1)-1)'*L(1)/n(1);
y  = x';
%z = z.*cvec(tukeywin(n,0.1));
%z = -dz_dx*(L/(2*pi))*sin(2*pi*x/L);
z  = dz_dx*(x - L(1)/2) + dz_dy*(y - L(2)/2);
%z  = z + 0*dz_dx*randn(n,1);
%z = circshift(z,-10);
mu = L/2;
p = 1;
%h0_ = hmin + hmax*(1./(1+((x-mu)/sd).^2));
mu = L /2;
hx = normpdf(x,mu(1),sdx)/normpdf(0,0,sd(1));
hy = normpdf(y,mu(2),sdy)/normpdf(0,0,sd(2));%/normpdf(0,0,s)^2;
h0 = hmin + hmax*(hx+hy);
h0 = hmin+hmax*(cvec(hx)*ones(1,n(2)) + ones(n(1),1)*rvec(hy));
h0 = h0';
%h0 = hmin+hmax*(1-(1-(normpdf(x,mu,sd)/normpdf(0,0,sd)).^p).^(1/p));
%h0 = hmax*(h0>0.5*hmax);
%h0 = meanfilt1(h0,3);
xc0 = sum(h0.*x)./sum(h0);
	
dx = x(2)-x(1);

K = 1.2e3;
C = 1e3;
T = 1;
dt = 1;
% 1 - dt*K*h*2/dx^2 > 0

method_C = {'central'}
%,'upwind','optimal','optimal-linearized'};
%method_C = {'upwind'}
%method_C = {'central'}
%,'optimal-linearized'};
%method_C ={'optimal'}

for jdx=1:length(method_C)
	method = method_C{jdx};
	h = h0;
	returnmat = true;
	m = method_C{jdx};

	%gnsolver = 'ilu-ks';
	z  = flat(z);
	h0 = flat(h0);
	h  = h0;
	returnmat         = true;
	
	resfun            = @(h) euler_implicit_residual(@(h) swe_simplified_2d_dhdt(h,z,K,L,n,method,returnmat), dt, h, h0);
	jfun              = @(h) euler_implicit_jacobian(@(h) swe_simplified_2d_jacobian(h,z,K,L,n,method,returnmat), dt, h, returnmat);
	J = jfun(h(:));
%	resfun            = @(h) euler_implicit_residual(@(h) swe_zero_inertia_1d_dhdt(h,z,C,L,n,method_C{jdx}), dt, h, h0);
%	J                 = euler_implicit_jacobian(@(h) swe_zero_inertia_1d_jacobian(h,z,C,L,n,m,returnmat), dt, h, returnmat);

	J_ = jacobian_unstructured(resfun,h0);


if (1) %pdx==1)
	zx = zeros(size(hx));
	hx0 = hx;
	resfun1            = @(h) euler_implicit_residual(@(h) swe_simplified_1d_dhdt(h,zx,K,L(1),n(1),method,returnmat), dt, h, hx0);
	jfun1              = @(h) euler_implicit_jacobian(@(h) swe_simplified_1d_jacobian(h,zx,K,L(1),n(1),method,returnmat), dt, h, returnmat);
	J1                 = jfun1(hx);
	J2                 = jfun1(hy);

figure(999);
J1_ = J1(1:n(1),1:n(1));
k = -1;
subplot(2,3,1)
cla
plot(diag(J1,-1))
hold on
plot(diag(J1_,-1),'--');
subplot(2,3,2)
cla
plot(diag(J1,0))
hold on
plot(diag(J1_,0),'--');
subplot(2,3,3)
cla
plot(diag(J1,+1));
hold on
plot(diag(J1_,+1),'--');
subplot(2,3,4)
plot([diag(J2,-1),NaN*diag(J1_,-1)]);
subplot(2,3,5)
plot([diag(J2,0),NaN*diag(J1_,0)]);
subplot(2,3,6)
plot([diag(J2,+1),NaN*diag(J1_,+1)]);
%clf
%plot([hx(:),hy(:)])
pause
end

dJ  = diag(J);
dJ = reshape(dJ,n);
dJ_ = diag(J_);
dJ_ = reshape(dJ_,n);


figure(2+jdx*10);
subplot(2,3,3*(zdx-1)+1);
imagesc(reshape(h0,n))
title('initial');

subplot(2,3,3*(zdx-1)+2);
imagesc(dJ);
title('analytic');

subplot(2,3,3*(zdx-1)+3);
imagesc(dJ_);
title('numerical');

figure(1+jdx*10);
subplot(2,2,2*(zdx-1)+1);
plot([dJ(:,round(end/2)),dJ_(:,round(end/2))]);
%ylim([-1,1]*1)

subplot(2,2,2*(zdx-1)+2);
plot([dJ(round(end/2),:)',dJ_(round(end/2),:)']);
%ylim([-1,1]*1)

%plot([diag(J),diag(J_)]-1)
%ylim([-1,1]*1)

if (0)
subplot(3,5,zdx+2);%kdx+(zdx-1)*10)
plot([diag(J,-1),diag(J_,-1)])
ylim([-1,1]*1)

subplot(3,5,zdx+10);%kdx+(zdx-1)*10)
plot([diag(J,+1),diag(J_,+1)])
ylim([-1,1]*1)
end

end % for jdx
end % for zdx
